# Use the os package to interact with the environment
import os
# Bring in Pandas for Dataframe functionality
import pandas as pd
# Numpy for basics
import numpy as np
# Use StringIO for working with file contents
from io import StringIO
# Enable IPython to display matplotlib graphs
import matplotlib.pyplot as plt
# Enable interaction with the FireCloud API
from firecloud import api as fapi
# Import the iPython HTML rendering for displaying links to Google Cloud Console
from IPython.display import display, HTML
# Import urllib modules for building URLs to Google Cloud Console
import urllib.parse
# BigQuery for querying data
from google.cloud import bigquery
#Import Sys
import sys as sys

import argparse, psutil, tqdm, math
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from pck.generics import SingletonFactory
from pck.log import logger as lg
from pck.gchutil import *
from configure import Settings 
from PIL import Image

def separate_sample_case_control():
    # Retrieve settings
    sf = SingletonFactory()
    s=sf.get_service_as('settings', Settings)

    # Check if previously processed
    if os.path.exists(os.path.join(s.PROJECT_DATA_DIR, 'case-control-variants.csv.gz')):
        return    

    # 1. Load data
    lg.info('Loading data...')
    df_master = pd.read_csv(os.path.join(s.CLINICAL_DATA_DIR,'master_key_release7_final.csv'))

    # 2. Case group filtering
    case_group = df_master[(df_master['pruned'] == 0) & (df_master['baseline_GP2_phenotype_for_qc'] == 'PD')]

    # 3. Control group filtering
    control_group = df_master[(df_master['pruned'] == 0) & (df_master['baseline_GP2_phenotype_for_qc'] == 'Control')]

    # 4. Extract necessary columns and add group names
    case_group = case_group[['GP2ID', 'GP2sampleID']].copy()
    case_group['Group'] = 'Case'
    control_group = control_group[['GP2ID', 'GP2sampleID']].copy()
    control_group['Group'] = 'Control'

    # 5. Merge Case and control group
    merged_group = pd.concat([case_group, control_group])

    # 6. Export to CSV
    lg.info('Saving step 1 results...')
    merged_group.to_csv(os.path.join(s.PROJECT_DATA_DIR, 'case-control-samples.csv.gz'), compression='gzip', index=False)

def separate_variant_case_control(classification='VUS'):      
    # Retrieve settings
    sf = SingletonFactory()
    s=sf.get_service_as('settings', Settings)

    # Check if previously processed
    if os.path.exists(os.path.join(s.PROJECT_DATA_DIR, 'case-control-variants.csv.gz')):
        return   
    
    # 1. Load data
    lg.info('Loading data...')
    df_samples = pd.read_csv(os.path.join(s.PROJECT_DATA_DIR, 'case-control-samples.csv.gz'))
    df_variants = pd.read_csv(os.path.join(s.USER_DATA_DIR, 'all_variants.csv.gz'))
    df_ancestry = pd.read_csv(os.path.join(s.USER_DATA_DIR, 'selected_merged.tsv.gz'), delimiter='\t', dtype=str)

    # 2. Assign Variant Key
    df_ancestry['@VARIANT_KEY'] = df_ancestry.apply(lambda row: f"{row['#CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", axis=1)
    df_variants_filtered = df_variants[df_variants['@Class'] == classification]

    # 3. Mapping IDs and Groups for Samples with Variants
    sample_dict = df_samples.set_index('GP2sampleID')['Group'].to_dict()

    output_data = []
    total_ancestry_rows = len(df_ancestry.index)
    with tqdm.tqdm(total=total_ancestry_rows, desc='Processing variants') as pbar:
        for idx, variant in df_ancestry.iterrows():
            variant_key = variant['@VARIANT_KEY']
            if variant_key in df_variants_filtered['@VARIANT_KEY'].values:
                for sample_id in df_ancestry.columns[57:]:  # Column index at which GP2sampleID starts
                    genotype = variant[sample_id]
                    if genotype in {'0/1', '1/0', '1/1'}:  # If the sample has a variant
                        group = sample_dict.get(sample_id, None)
                        if group:
                            row = {
                                'GP2ID': df_samples[df_samples['GP2sampleID'] == sample_id]['GP2ID'].values[0],
                                'GP2sampleID': sample_id,
                                'Group': group,
                                'Case': group == 'Case',
                                'Control': group == 'Control',
                                '@VARIANT_KEY': variant_key,
                                '@Class': classification,
                            }
                            row.update(variant[:57].to_dict())  # Variants 컬럼 추가
                            output_data.append(row)
            pbar.update(1)

    # 4. 결과 DataFrame 생성 및 저장
    lg.info('Saving step 2 results...')
    df_output = pd.DataFrame(output_data)
    df_output.to_csv(os.path.join(s.PROJECT_DATA_DIR, 'case-control-variants.csv.gz'), compression='gzip', index=False)

def plot_results_separate(output_png_filename, item_per_plot):
    # Retrieve settings
    sf = SingletonFactory()
    s=sf.get_service_as('settings', Settings)

    # Load data
    df_variants = pd.read_csv(os.path.join(s.PROJECT_DATA_DIR, 'case-control-variants.csv.gz'))

    # Calculate variant frequency
    variant_counts = df_variants.groupby(['@VARIANT_KEY', 'Group']).size().unstack(fill_value=0).reset_index()
    variant_counts['Total'] = variant_counts['Case'] + variant_counts['Control']
    variant_counts = variant_counts.sort_values(by='Total', ascending=False).drop(columns='Total')
    variant_counts_melted = variant_counts.melt(id_vars='@VARIANT_KEY', value_vars=['Case', 'Control'], var_name='Group', value_name='Frequency')

    # Export melted dataframe
    variant_counts_melted.to_csv(os.path.join(s.PROJECT_DATA_DIR, 'variants_counts_melted.csv'), index=False)

    # Divide the variant key by n to plot the graph (split into multiple graphs for readability)
    variants_per_plot = item_per_plot
    num_variants = len(variant_counts['@VARIANT_KEY'].unique())
    num_plots = math.ceil(num_variants / variants_per_plot)
    max_y = variant_counts_melted['Frequency'].max() + 50

    # Create a subplot
    fig, axes = plt.subplots(num_plots, 1, figsize=(18, num_plots * 10), sharey=True)
    axes = axes.flatten()  # Convert to list

    for i in range(num_plots):
        start_idx = i * variants_per_plot
        end_idx = start_idx + variants_per_plot
        subset = variant_counts_melted[variant_counts_melted['@VARIANT_KEY'].isin(variant_counts['@VARIANT_KEY'][start_idx:end_idx])]

        ax = axes[i]
        sns.barplot(x='@VARIANT_KEY', y='Frequency', hue='Group', data=subset, ax=ax)
        
        # Fixing the tick labels
        ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
        
        # Adding data labels
        for p in ax.patches:
            height = p.get_height()
            
            if height == 0:
                continue

            ax.text(p.get_x() + 0.05 + p.get_width() / 2., height+5,
                    f'{int(height)}', ha='center', va='bottom', fontsize=9, rotation=90)
        
        # Adding vertical grid
        unique_variants = subset['@VARIANT_KEY'].unique()            
        for x in range(0, len(unique_variants)+1):
            # Adjust x position for the vertical line
            ax.axline((x-0.5, 0), (x-0.5, 1), color='grey', alpha=0.3, linewidth=0.5, linestyle='-.')

        ax.set_xlabel('')
        ax.set_title(f'Variants {start_idx + 1} to {end_idx if end_idx < num_variants else num_variants}')
        
        # Save each subplot
        fig_subplot = plt.figure(figsize=(18, 10))
        ax_subplot = fig_subplot.add_subplot(111)
        sns.barplot(x='@VARIANT_KEY', y='Frequency', hue='Group', data=subset, ax=ax_subplot)
        ax_subplot.xaxis.set_major_locator(FixedLocator(ax_subplot.get_xticks()))
        ax_subplot.set_xticklabels(ax_subplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        for p in ax_subplot.patches:
            height = p.get_height()            
            if height == 0:
                continue
            ax_subplot.text(p.get_x() + 0.05 + p.get_width() / 2., height+5,
                            f'{int(height)}', ha='center', va='bottom', fontsize=9, rotation=90)
        for x in range(0, len(unique_variants)+1):
            # Adjust x position for the vertical line
            ax_subplot.axline((x-0.5, 0), (x-0.5, 1), color='grey', alpha=0.3, linewidth=0.5, linestyle='-.')
        ax_subplot.set_ylim(0, max_y)
        ax_subplot.set_xlabel('')
        ax_subplot.set_title(f'Variants {start_idx + 1} to {end_idx if end_idx < num_variants else num_variants}')

        plt.tight_layout()
        # Save cropped image
        base_name = os.path.basename(output_png_filename)
        # Split the filename and extension
        name, extension = os.path.splitext(base_name)
        fig_subplot.savefig(os.path.join(s.PROJECT_DATA_DIR,f'{name}_{i+1}{extension}'), dpi=192)
        plt.close(fig_subplot)

    plt.tight_layout()
    lg.info('Saving step 3 results...')
    plt.savefig(os.path.join(s.PROJECT_DATA_DIR, f'{output_png_filename}'), dpi=192)
    lg.info(f'The result was saved as {output_png_filename} in {s.PROJECT_DATA_DIR}.')

    return num_plots

def main():
    lg.info('GP2 v7.0')
    lg.info(f'Path to GP2 v7.0 Clinical Data: {GP2_CLINICAL_RELEASE_PATH}')
    lg.info(f'Path to GP2 v7.0 Raw Genotype Data: {GP2_RAW_GENO_PATH}')
    lg.info(f'Path to GP2 v7.0 Imputed Genotype Data: {GP2_IMPUTED_GENO_PATH}')

    if not os.path.exists(s.PROJECT_DIR):
        os.mkdir(s.PROJECT_DIR)
    if not os.path.exists(s.TOOL_DIR):
        os.mkdir(s.TOOL_DIR)
    if not os.path.exists(s.DATA_ROOT):
        os.mkdir(s.DATA_ROOT)
    if not os.path.exists(s.CLINICAL_DATA_DIR):
        os.mkdir(s.CLINICAL_DATA_DIR)        
    if not os.path.exists(s.USER_DATA_DIR):
        os.mkdir(s.USER_DATA_DIR)        
    if not os.path.exists(s.PROJECT_DATA_DIR):
        os.mkdir(s.PROJECT_DATA_DIR)
    
    # Preparation: Copy files required from GCS
    if not os.path.exists(os.path.join(s.CLINICAL_DATA_DIR, 'master_key_release7_final.csv')):
        shell_do(f'gsutil -u {s.BILLING_PROJECT_ID} -m cp {s.GP2_CLINICAL_RELEASE_PATH}/master_key_release7_final.csv {s.CLINICAL_DATA_DIR}')
    if not os.path.exists(os.path.join(s.USER_DATA_DIR,'selected_merged.tsv.gz')):
        shell_do(f'gsutil cp "gs://fc-66ffe54f-b734-406f-8f37-ec2292017a9e/user_data/selected_merged.tsv.gz" {os.path.join(s.USER_DATA_DIR, "selected_merged.tsv.gz")}')
    if not os.path.exists(os.path.join(s.USER_DATA_DIR,'all_variants.csv.gz')):
        shell_do(f'gsutil cp "gs://fc-66ffe54f-b734-406f-8f37-ec2292017a9e/user_data/all_variants.csv.gz" {os.path.join(s.USER_DATA_DIR, "all_variants.csv.gz")}')

    # Step 1: Categorize samples into Case and Control groups
    lg.info('Step 1: Categorize samples into Case and Control groups')
    separate_sample_case_control()
    
    # Step 2: Assign variants into Case and Control samples
    lg.info('Step 2: Assign variants into Case and Control samples')
    separate_variant_case_control('VUS')

    # Step 3: Visualize data
    lg.info('Step 3: Visualize data')
    plot_results_separate('VUS_frequency_by_case_control.jpg', 66)

if __name__ == '__main__':
    sf = SingletonFactory()
    sf.add_service('settings', Settings())
    parser = argparse.ArgumentParser(description = 'GCH1 Project on GP2')    

    #group = parser.add_mutually_exclusive_group(required = True)    
    #group.add_argument("-tsv", "--vcf-to-tsv", help = "Converts ancestry VCFs to tsv", action = 'store_true')
    #group.add_argument("-e", "--ancestry", type = str, help = "Enter one of the following: AAC, AFR, AJ, AMR, CAH, CAS, EAS, EUR, FIN, MDE, or SAS.", choices = ['AAC', 'AFR', 'AJ', 'AMR', 'CAH', 'CAS', 'EAS', 'EUR', 'FIN', 'MDE', 'SAS'])
    #group.add_argument("-wgs", "--using-wgs", help = "Collect variants from WGS data", action = 'store_true')
    parser.add_argument("-id", "--billing-project-id", type = str, help = "GP2 Terra Billing Project ID")
    parser.add_argument("-ns", "--workspace-namespace", type = str, help = "GP2 Terra Namespace of Workspace")
    parser.add_argument("-ws", "--workspace-name", type = str, help = "GP2 Terra Workspace Name")
    parser.add_argument("-gcs", "--workspace-bucket", type = str, help = "GP2 Terra Workspace bucket (Google Cloud Storage)")
    args = parser.parse_args()

    BILLING_PROJECT_ID = args.billing_project_id if args.billing_project_id is not None else os.environ['GOOGLE_PROJECT']
    WORKSPACE_NAMESPACE = args.workspace_namespace if args.workspace_namespace is not None else os.environ['WORKSPACE_NAMESPACE']
    WORKSPACE_NAME = args.workspace_name if args.workspace_name is not None else os.environ['WORKSPACE_NAME']
    WORKSPACE_BUCKET = args.workspace_bucket if args.workspace_bucket is not None else os.environ['WORKSPACE_BUCKET'] 
    GP2_RELEASE_PATH = "gs://gp2tier2/release7_30042024"
    GP2_CLINICAL_RELEASE_PATH = f"{GP2_RELEASE_PATH}/clinical_data"
    GP2_META_RELEASE_PATH = f'{GP2_RELEASE_PATH}/meta_data'
    GP2_SUMSTAT_RELEASE_PATH = f'{GP2_RELEASE_PATH}/summary_statistics'
    GP2_RAW_GENO_PATH = f"{GP2_RELEASE_PATH}/raw_genotypes"
    GP2_IMPUTED_GENO_PATH = f"{GP2_RELEASE_PATH}/imputed_genotypes"
    GP2_WGS_VCF_PATH = "gs://gp2tier2/release6_21122023/wgs/var_calling/deepvariant"
    AMP_RELEASE_PATH = "gs://amp-pd-data/releases/2022_v3release_1115"
    AMP_CLINICAL_RELEASE_PATH = f"{AMP_RELEASE_PATH}/clinical"
    AMP_WGS_RELEASE_PATH = "gs://amp-pd-genomics/releases/2022_v3release_1115/wgs-WB-DWGS"
    AMP_WGS_RELEASE_PLINK_PATH = os.path.join(AMP_WGS_RELEASE_PATH, 'plink')
    AMP_WGS_RELEASE_PLINK_PFILES = os.path.join(AMP_WGS_RELEASE_PLINK_PATH, 'pfiles')    
    #WORKSPACE_ATTRIBUTES = fapi.get_workspace(WORKSPACE_NAMESPACE, WORKSPACE_NAME).json().get('workspace',{}).get('attributes',{})

    s=sf.get_service_as('settings', Settings)
    s.HOME_DIR = os.path.expanduser("~")
    s.PROJECT_DIR = os.path.join(s.HOME_DIR, 'data', "gch1-workdir")
    s.TOOL_DIR = os.path.join(s.PROJECT_DIR, 'tools')    
    s.DATA_ROOT = os.path.join(s.PROJECT_DIR, 'data')
    s.CLINICAL_DATA_DIR = os.path.join(s.DATA_ROOT, 'clinical_data')    
    s.PROJECT_DATA_DIR = os.path.join(s.DATA_ROOT, 'project_generated_files')
    s.USER_DATA_DIR = os.path.join(s.DATA_ROOT, 'user_data')
    s.BILLING_PROJECT_ID = BILLING_PROJECT_ID
    s.WORKSPACE_NAMESPACE = WORKSPACE_NAMESPACE
    s.WORKSPACE_NAME = WORKSPACE_NAME
    s.WORKSPACE_BUCKET = WORKSPACE_BUCKET
    s.GP2_WGS_VCF_PATH = GP2_WGS_VCF_PATH
    s.GP2_RELEASE_PATH = GP2_RELEASE_PATH
    s.GP2_CLINICAL_RELEASE_PATH = GP2_CLINICAL_RELEASE_PATH
    s.GP2_RAW_GENO_PATH = GP2_RAW_GENO_PATH
    s.GP2_IMPUTED_GENO_PATH = GP2_IMPUTED_GENO_PATH
    s.AMP_RELEASE_PATH = AMP_RELEASE_PATH
    s.AMP_CLINICAL_RELEASE_PATH = AMP_CLINICAL_RELEASE_PATH
    s.AMP_WGS_RELEASE_PATH = AMP_WGS_RELEASE_PATH
    s.AMP_WGS_RELEASE_PLINK_PATH = AMP_WGS_RELEASE_PLINK_PATH
    s.AMP_WGS_RELEASE_PLINK_PFILES = AMP_WGS_RELEASE_PLINK_PFILES     
    s.GP2_BUCKET = 'gp2tier2'
    s.AMP_BUCKET = 'amp-pd-genomics'
    s.TOTAL_MEM_SIZE = 2 ** (int((psutil.virtual_memory().total / (1024 ** 3))).bit_length() - 1)
      
    main()