from pyplink import PyPlink
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Read the sample information
psam_df = pd.read_csv('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.psam', sep='\s+')
# Read the variant information
pvar_df = pd.read_csv('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.pvar', sep='\s+')

if not os.path.exists('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_matrix.pkl') or not os.path.exists('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_df.pkl'):
    # Load the genomic data using PyPLINK
    ctx = PyPlink('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos')

    # Initialize a list to store genotype data
    genotypes = []
    # Extract genotype data for each SNP
    for snp, geno in ctx.iter_geno():
        genotypes.append(geno)

    # Convert to DataFrame and transpose so that the rows are the samples
    #geno_df = pd.DataFrame.from_records(genotypes, columns=['snp', 'genotype']).set_index('snp').T
    geno_df = pd.DataFrame(genotypes, index=pvar_df['ID']).T
    geno_df.columns = pvar_df['ID']
    geno_df.index = psam_df['#IID']

    # Only retian numeric part of the genotype
    geno_matrix = geno_df.applymap(lambda x: x if isinstance(x, int) else np.nan)
    geno_df.to_pickle('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_df.pkl')
    geno_matrix.to_pickle('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_matrix.pkl')
else:    
    geno_df = pd.read_pickle('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_df.pkl')
    geno_matrix = pd.read_pickle('/home/jereiard/data/gch1-workdir/data/pca/pruned_data_filtered_aos.geno_matrix.pkl')

# Read the eigenvector file
eigenvec_df = pd.read_csv('/home/jereiard/data/gch1-workdir/data/pca/gch1_filtered_pca_result_aos.eigenvec', sep='\s+')
#eigenvec_df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 3)]

# Read the eigenvalues file
eigenval_df = pd.read_csv('/home/jereiard/data/gch1-workdir/data/pca/gch1_filtered_pca_result_aos.eigenval', sep='\s+', header=None)
eigenval_df.columns = ['eigenval']

# Extract PC1 scores
pc1_scores = eigenvec_df['PC1']

# Convert PC1 scores to numpy array for calculation
pc1_scores_array = np.array(pc1_scores)

# Calculate the correlation between each SNP and PC1 scores
# Handle missing values and impute if necessary
geno_matrix = geno_df.apply(pd.to_numeric, errors='coerce')

correlations = geno_matrix.apply(lambda var: np.corrcoef(var, pc1_scores_array)[0, 1], axis=0)

# Get the top 10 variants with the highest absolute correlation
top_10_variants = correlations.abs().sort_values(ascending=False).head(10).index
print("10 variables with the highest explanatory power in PC1:", top_10_variants)

top_n_variants = correlations.abs().sort_values(ascending=False).head(5000).index
geno_matrix_top_n = geno_matrix[top_n_variants]

psam_df['PHENO1'] = psam_df['PHENO1'].replace({1: 'Early onset', 2: 'Late onset > 50yrs'})
merged_df = pd.merge(eigenvec_df, psam_df[['#IID', 'PHENO1']], on='#IID')

# Extract PC1 and PC2 scores
pc1_scores = merged_df['PC1']
pc2_scores = merged_df['PC2']

# Calculate the proportion of variance explained by PC1 and PC2
variance_explained_pc1 = eigenval_df['eigenval'][0] / eigenval_df['eigenval'].sum() * 100
variance_explained_pc2 = eigenval_df['eigenval'][1] / eigenval_df['eigenval'].sum() * 100

# Plot PC1 and PC2 with group coloring
plt.figure(figsize=(10, 7))

for group in merged_df['PHENO1'].unique():
    group_df = merged_df[merged_df['PHENO1'] == group]
    plt.scatter(group_df['PC1'], group_df['PC2'], label=group, edgecolor='k')

plt.xlabel(f'PC1 ({variance_explained_pc1:.2f}%)')
plt.ylabel(f'PC2 ({variance_explained_pc2:.2f}%)')
plt.title('PCA Plot by Group')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('/home/jereiard/data/gch1-workdir/data/pca/pca_aos.jpg', dpi=192)
#plt.show()

# Perform PCA on the top 1000 variants
pca = PCA(n_components=2)
pca_result = pca.fit_transform(geno_matrix_top_n.fillna(0))  # Replace NaN with 0 or an appropriate imputation value

# Create a DataFrame with the PCA results
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['#IID'] = geno_matrix_top_n.index

# Merge group information with PCA results
#psam_df['PHENO1'] = psam_df['PHENO1'].replace({1: 'Group1', 2: 'Group2'})
merged_df = pd.merge(pca_df, psam_df[['#IID', 'PHENO1']], on='#IID')

# Calculate the proportion of variance explained by PC1 and PC2 (if needed)
variance_explained_pc1 = pca.explained_variance_ratio_[0] * 100
variance_explained_pc2 = pca.explained_variance_ratio_[1] * 100

# Plot PC1 and PC2 with group coloring
plt.figure(figsize=(10, 7))

for group in merged_df['PHENO1'].unique():
    group_df = merged_df[merged_df['PHENO1'] == group]
    plt.scatter(group_df['PC1'], group_df['PC2'], label=group, edgecolor='k')

plt.xlabel(f'PC1 ({variance_explained_pc1:.2f}%)')
plt.ylabel(f'PC2 ({variance_explained_pc2:.2f}%)')
plt.title('PCA Plot by Group (Top 1000 Variants)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('/home/jereiard/data/gch1-workdir/data/pca/pca_aos_top1000.jpg', dpi=192)
#plt.show()