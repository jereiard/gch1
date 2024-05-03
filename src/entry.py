import os
import subprocess
import argparse
import shlex
import uuid
import psutil
import multiprocessing
import pandas as pd
from pathlib import Path
from pck.gchutil import *
from pck.log import logger as lg

class C:
  def __init__(self, f, before=None, after=None):
    self.f=f
    self.before=before
    self.after=after
    
class B:
    def __init__(self, f): self.f=f

class Pipe  (B): __ror__=lambda self, x: self.f(x)
# class Map   (B): __ror__=lambda self, x: map   (self.f, x)
# class Filter(B): __ror__=lambda self, x: filter(self.f, x)
class Pipe2  (C): __ror__=lambda self, x: self.f(x, self.before, self.after)

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

class Singleton:
    _instance=None

    def __init__(self, value):
        self.value=value

    def get_value(self):
        return Struct(**self.value)
    def set_value(self, struct):
        self.value=struct.__dict__

class SingletonFactory:
    _instance=None

    @staticmethod
    def get_singleton_instance(value):
        if SingletonFactory._instance is None:
            SingletonFactory._instance=Singleton(value)
        return SingletonFactory._instance
    
singletonFactory=SingletonFactory()
settings=singletonFactory.get_singleton_instance({})

def reheader_tsv_py(input, reheader, dummy=None):
    s=settings.get_value()
    body=os.path.join(s.dataDir, uuid.uuid1().hex)
    header=os.path.join(s.dataDir, uuid.uuid1().hex)
    output=append_id(input, "reheader")
    workDir=os.getcwd()
    
    try:
        lg.info("Changing column headers...")
        if not is_null_or_small(output, 10240):
            return output
        os.chdir(s.dataDir)  

        os.system(f"tail -n +2 {input} > {body}")
        os.system(f"head -n 1 {input} > {header}")
        with open(header, "rt") as f:
            line = f.readline()

        for item in reheader:
            line = line.replace(item[0], item[1])
        
        with open(header, "wt") as f:
            f.write(line)

        os.system(f"cat {header} {body} > {output}")

        os.remove(header)
        os.remove(body)

        if is_null_or_small(output, 10240):
            raise Exception(f"{output} file(s) may be corrupted.")
        
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def reheader_tsv(input, reheader, dummy=None):
    s=settings.get_value()
    body=os.path.join(s.dataDir, uuid.uuid1().hex)
    header=os.path.join(s.dataDir, uuid.uuid1().hex)
    output=append_id(input, "reheader")
    workDir=os.getcwd()
    
    try:
        lg.info("Changing column headers...")
        if not is_null_or_small(output, 10240):
            return output
        os.chdir(s.dataDir)  

        os.system(f"tail -n +2 {input} > {body}")
        os.system(f"head -n 1 {input} > {header}")
        for item in reheader:
            os.system(f'sed -i -e "s/{item[0]}/{item[1]}/g" {header}')
        os.system(f"cat {header} {body} > {output}")

        os.remove(header)
        os.remove(body)

        if is_null_or_small(output, 10240):
            raise Exception(f"{output} file(s) may be corrupted.")
        
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def to_tsv(input, fields, dest_dir=None):
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    effOnePerLineBin=os.path.join(snpeffDir, "snpEff", "scripts/vcfEffOnePerLine.pl")
    snpsiftBin=os.path.join(snpeffDir, "snpEff", "SnpSift.jar")
    output=change_ext(append_id(input, "extractField"), ".tsv")
    workDir=os.getcwd()

    try:
        lg.info("Converting VCF to TSV...")
        if not is_null_or_small(output, 10240):
            return output
        os.chdir(s.dataDir)   

        if Path(input).suffix != ".gz":
            cmd = f'cat {input} | {effOnePerLineBin} | java -Xmx{s.totalMemSizeGB}g -jar {snpsiftBin} extractFields -v -e "." - {fields} > {output}'
        else:
            cmd = f'zcat {input} | {effOnePerLineBin} | java -Xmx{s.totalMemSizeGB}g -jar {snpsiftBin} extractFields -v -e "." - {fields} > {output}'
        
        status = os.system(cmd)
        if status != 0: 
            raise RuntimeError("Convert VCF to TSV using SnpSift was terminated unexpectedly.")

        if is_null_or_small(output, 10240):
            raise Exception(f"{output} file(s) may be corrupted.")

    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def rename_tags(input, before=None, after=None):
    s=settings.get_value()
    output=append_id(input, "retag")
    workDir=os.getcwd()

    try:
        lg.info("Renaming tags...")
        if not is_null_or_small(output):
            return output
        os.chdir(s.dataDir)   

        if before != None and after != None:
            status = shell_do_redir_stdout(f"bcftools annotate -c {after}:={before} {input}", output)
        elif before != None and after == None:
            status = shell_do_redir_stdout(f'bcftools annotate -c {before} {input}', output)
        
        if status != 0: 
            raise RuntimeError("VCF tag rename using bcftools was terminated unexpectedly.")

        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")

    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def gnomad_filter(input, af, dummy=None):
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    snpsiftBin=os.path.join(snpeffDir, "snpEff", "SnpSift.jar")
    output=append_id(input, "fltr")
    workDir=os.getcwd()

    try:
        lg.info(f"Filtering variants using mutant allele frequency {af}...")
        if is_null_or_small(output) == False:
            return output
        os.chdir(s.dataDir)

        status = os.system(
            f'cat {input} | java -Xmx{s.totalMemSizeGB}g -jar {snpsiftBin} filter "( AF < {af} ) & ( AF_eas < {af} ) & ( AF_mid < {af} ) & ( AF_sas < {af} ) & ( AF_afr < {af} ) & ( AF_amr < {af} ) & ( AF_nfe < {af} )" > {output}'
        )
        if status != 0: 
            raise RuntimeError("Filtration was terminated unexpectedly.")   
        
        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")
        
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def gnomad_genome(input):
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    snpsiftBin=os.path.join(snpeffDir, "snpEff", "SnpSift.jar")
    #gnomadExome=os.path.join(snpeffDir, "snpEff", "data", "gnomadv4", "gnomad.exomes.v4.1.sites.chr14.vcf.bgz")
    gnomadGenome=os.path.join(snpeffDir, "snpEff", "data", "gnomadv4", "gnomad.genomes.v4.1.sites.chr14.vcf.bgz")
    output=append_id(input, "genome")
    workDir=os.getcwd()

    try:
        lg.info(f"Annotating VCF(s) using SnpSift with gnomAD Genome v4.1...")
        if not is_null_or_small(output):
            return output
        os.chdir(s.dataDir)   

        status=shell_do_redir_stdout(f"java -Xmx{s.totalMemSizeGB}g -jar {snpsiftBin} annotate -info AF,AF_eas,AF_eas,AF_sas,AF_mid,AF_afr,AF_amr,AF_nfe,AF_fin,AF_asj,AC,AC_eas,AC_sas,AC_mid,AC_afr,AC_amr,AC_nfe,AC_fin,AC_asj,AN,nhomalt -v {gnomadGenome} {input}", output)
        if status != 0: 
            raise RuntimeError("Annotation was terminated unexpectedly.")
        
        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")
                    
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output
    
def gnomad_exome(input):
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    snpsiftBin=os.path.join(snpeffDir, "snpEff", "SnpSift.jar")
    gnomadExome=os.path.join(snpeffDir, "snpEff", "data", "gnomadv4", "gnomad.exomes.v4.1.sites.chr14.vcf.bgz")
    #gnomadGenome=os.path.join(snpeffDir, "snpEff", "data", "gnomadv4", "gnomad.genomes.v4.1.sites.chr14.vcf.bgz")
    output=append_id(input, "exome")
    workDir=os.getcwd()

    try:
        lg.info(f"Annotating VCF(s) using SnpSift with gnomAD Exome v4.1...")
        if not is_null_or_small(output):
            return output
        os.chdir(s.dataDir)   

        status=shell_do_redir_stdout(f"java -Xmx{s.totalMemSizeGB}g -jar {snpsiftBin} annotate -info AF,AF_eas,AF_eas,AF_sas,AF_mid,AF_afr,AF_amr,AF_nfe,AF_fin,AF_asj,AC,AC_eas,AC_sas,AC_mid,AC_afr,AC_amr,AC_nfe,AC_fin,AC_asj,AN,nhomalt -v {gnomadExome} {input}", output)
        if status != 0: 
            raise RuntimeError("Annotation was terminated unexpectedly.")
        
        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")
                    
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def annotate():
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    snpeffBin=os.path.join(snpeffDir, "snpEff", "snpEff.jar")
    snpeffDB=os.path.join(snpeffDir, "snpEff", "data")
    input=os.path.join(s.dataDir, f"pheno_GCH1_{s.ancestry}_final.vcf")
    output=append_id(input, "ann")
    workDir=os.getcwd()

    try:
        lg.info(f"Annotating VCF(s) using snpEff with GRCh38 MANE 1.2 RefSeq...")
        if not is_null_or_small(output):
            return output
        os.chdir(s.dataDir)

        status=shell_do_redir_stdout(f"java -Xmx{s.totalMemSizeGB}g -jar {snpeffBin} -v -stats pheno_GCH1_{s.ancestry}_final.ann.html GRCh38.mane.1.2.refseq {input}", output)
        if status != 0: 
            raise RuntimeError("Annotation was terminated unexpectedly.")
        
        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")
                    
    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def merge_vcf():
    s=settings.get_value()
    workDir=os.getcwd()

    try:
        targets=[]
        output=os.path.join(s.dataDir, "ancestry.merged.vcf.gz")

        for path in Path(s.dataDir).rglob("*genome.fltr.retag.vcf"):
            shell_do(f'bcftools sort -Oz -o {str(path)}.gz {str(path)}')
            #shell_do(f'bgzip --threads {multiprocessing.cpu_count()} --force {str(path)}')

        for path in Path(s.dataDir).rglob("*genome.fltr.retag.vcf.gz"):
            lg.info(path.resolve())
            shell_do(f'bcftools index --threads {multiprocessing.cpu_count()} --force {str(path)}')
            targets.append(path.resolve())
            
        cmd = f'bcftools merge --threads {multiprocessing.cpu_count()} -m none -Oz {" ".join([str(file) for file in targets])}'
        lg.info(cmd)
        status = shell_do_redir_stdout(cmd, output)
        if status != 0: 
            raise RuntimeError("Merging VCFs using bcftools was terminated unexpectedly.")
        
        status = shell_do(f'bcftools index --threads {multiprocessing.cpu_count()} --force {output}')
        if status != 0: 
            raise RuntimeError("Indexing the merged VCF using bcftools was terminated unexpectedly.")
        
        if is_null_or_small(output):
            raise Exception(f"{output} file(s) may be corrupted.")

    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
        return output

def extract_gch1():
    s=settings.get_value()
    workDir=os.getcwd()

    try:
        os.chdir(s.dataDir)
        plink2Bin=os.path.join(s.toolDir, "plink2", "plink2")
        shell_do(f'{plink2Bin} --pfile chr14_{s.ancestry}_release7 --chr 14 --from-bp 54842007 --to-bp 54902836 --make-bed --out pheno_GCH1_{s.ancestry}')
        shell_do(f'{plink2Bin} --bfile pheno_GCH1_{s.ancestry} --export vcf --out pheno_GCH1_{s.ancestry}_final')

    except Exception as ex:
        os.removedirs(s.dataDir)
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)
    
def copy_data():
    s=settings.get_value()    
    workDir=os.getcwd()

    try:
        if not os.path.exists(s.dataDir):
            os.makedirs(s.dataDir)
        
        #shell_do(f'gsutil -u {s.billingProjectID} ls {s.gp2ImputedGenoPath}')
        #shell_do(f'gsutil -u {s.billingProjectID} ls -l {s.gp2ImputedGenoPath}/{s.ancestry}/')
        if (not os.path.exists(os.path.join(s.dataDir, f"chr14_{s.ancestry}_release7.log")) or
            not os.path.exists(os.path.join(s.dataDir, f"chr14_{s.ancestry}_release7.pgen")) or
            not os.path.exists(os.path.join(s.dataDir, f"chr14_{s.ancestry}_release7.psam")) or
            not os.path.exists(os.path.join(s.dataDir, f"chr14_{s.ancestry}_release7.pvar"))):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ImputedGenoPath}/{s.ancestry}/chr14_{s.ancestry}_release7* {s.dataDir}')
        #shell_do(f'gsutil -u {s.billingProjectID} ls {s.gp2ImputedGenoPath}')
        if not os.path.exists(os.path.join(s.dataDir, "GP2_clinical_data_dictionary.xlsx")):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ClinicalReleasePath}/GP2_clinical_data_dictionary.xlsx {s.dataDir}')
        if not os.path.exists(os.path.join(s.dataDir, "clinical_data_master_release7.csv")):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ClinicalReleasePath}/clinical_data_master_release7.csv {s.dataDir}')
        if not os.path.exists(os.path.join(s.dataDir, "clinical_data_summary_release7.xlsx")):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ClinicalReleasePath}/clinical_data_summary_release7.xlsx {s.dataDir}')
        if not os.path.exists(os.path.join(s.dataDir, "master_key_release7_data_dictionary.csv")):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ClinicalReleasePath}/master_key_release7_data_dictionary.csv {s.dataDir}')
        if not os.path.exists(os.path.join(s.dataDir, "master_key_release7_final.csv")):
            shell_do(f'gsutil -u {s.billingProjectID} -m cp -r {s.gp2ClinicalReleasePath}/master_key_release7_final.csv {s.dataDir}')

        lg.info(f"Data downloaded in {s.dataDir}")
        
    except Exception as ex:
        os.removedirs(s.dataDir)
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)

def install_snpeff():
    s=settings.get_value()
    snpeffDir=os.path.join(s.toolDir, "snpeff")
    workDir=os.getcwd()

    try:
        lg.info(f"Check whether snpEff is installed and install it if necessary.")
        if not os.path.exists(s.toolDir):
            os.makedirs(s.toolDir)
        if os.path.exists(snpeffDir):
            lg.info(f"snpEff is already installed in {snpeffDir}")
        else:
            os.makedirs(snpeffDir)
            os.chdir(snpeffDir)

            status=shell_do("curl -O https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip")
            if status != 0: 
                raise RuntimeError("Cannot download snpEff.")
            
            status=shell_do("unzip -o snpEff_latest_core.zip")
            if status != 0: 
                raise RuntimeError("Cannot install snpEff.")

            lg.info(f"snpEff downloaded and unzipped in {snpeffDir}")

        snpeffBin=os.path.join(snpeffDir, "snpEff", "snpEff.jar")
        if not os.path.exists(snpeffBin):
            raise RuntimeError("Cannot find a snpEff.jar runtime.")       
        
        snpeffDB=os.path.join(snpeffDir, "snpEff", "data")
        if not os.path.exists(snpeffDB):
            os.makedirs(snpeffDB)

        if not os.path.exists(os.path.join(snpeffDB, "GRCh38.mane.1.2.refseq")):
            lg.info(f"Installing GRCh38 databases...")
            status=shell_do(f"java -jar {snpeffBin} download -v GRCh38.mane.1.2.refseq")
            if status != 0: 
                raise RuntimeError("Cannot install GRCh38.mane.1.2.refseq database for snpEff.")
        
        snpeffGnomAD=os.path.join(snpeffDB, "gnomadv4")
        if not os.path.exists(snpeffGnomAD):
            lg.info(f"Installing gnomAD databases...")
            os.makedirs(snpeffGnomAD)
            os.chdir(snpeffGnomAD)

            if not os.path.exists(os.path.join(snpeffGnomAD, "gnomad.genomes.v4.1.sites.chr14.vcf.bgz")):
                status=shell_do(f"curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr14.vcf.bgz")
                if status != 0: 
                    raise RuntimeError("Cannot download gnomad.genomes.v4.1.sites.chr14.vcf.bgz")
                status=shell_do(f"curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr14.vcf.bgz.tbi")
                if status != 0: 
                    raise RuntimeError("Cannot download gnomad.genomes.v4.1.sites.chr14.vcf.bgz.tbi")
                
            if not os.path.exists(os.path.join(snpeffGnomAD, "gnomad.exomes.v4.1.sites.chr14.vcf.bgz")):
                status=shell_do(f"curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr14.vcf.bgz")
                if status != 0: 
                    raise RuntimeError("Cannot download gnomad.exomes.v4.1.sites.chr14.vcf.bgz")
                status=shell_do(f"curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr14.vcf.bgz.tbi")
                if status != 0: 
                    raise RuntimeError("Cannot download gnomad.exomes.v4.1.sites.chr14.vcf.bgz.tbi")

    except Exception as ex:
        os.removedirs(snpeffDir)
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)

def install_plink2():
    s=settings.get_value()  
    plinkDir=os.path.join(s.toolDir, "plink2")
    workDir=os.getcwd()

    try:
        lg.info(f"Check whether Plink2 is installed and install it if necessary.")
        if not os.path.exists(s.toolDir):
            os.makedirs(s.toolDir)
        if os.path.exists(plinkDir):
            lg.info(f"Plink2 is already installed in {plinkDir}")
        else:
            os.makedirs(plinkDir)
            os.chdir(plinkDir)

            status=shell_do("curl -O https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip")
            if status != 0: 
                raise RuntimeError("Cannot download plink2.")
            
            status=shell_do("unzip -o plink2_linux_avx2_20240105.zip")
            if status != 0: 
                raise RuntimeError("Cannot install plink2.")

            lg.info(f"plink2 downloaded and unzipped in {plinkDir}")

    except Exception as ex:
        os.removedirs(plinkDir)
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)

def main(args):
    ANCESTRY=args.ancestry
    if args.ancestry!=None:
        BILLING_PROJECT_ID=args.billing_project_id if args.billing_project_id is not None else os.environ["GOOGLE_PROJECT"]
        WORKSPACE_NAMESPACE=args.workspace_namespace if args.workspace_namespace is not None else os.environ["WORKSPACE_NAMESPACE"]
        WORKSPACE_NAME=args.workspace_name if args.workspace_name is not None else os.environ["WORKSPACE_NAME"]
        WORKSPACE_BUCKET=args.workspace_bucket if args.workspace_bucket is not None else os.environ["WORKSPACE_BUCKET"]

    GP2_RELEASE_PATH="gs://gp2tier2/release7_30042024"
    GP2_CLINICAL_RELEASE_PATH=f"{GP2_RELEASE_PATH}/clinical_data"
    GP2_RAW_GENO_PATH=f"{GP2_RELEASE_PATH}/raw_genotypes"
    GP2_IMPUTED_GENO_PATH=f"{GP2_RELEASE_PATH}/imputed_genotypes"

    AMP_RELEASE_PATH="gs://amp-pd-data/releases/2022_v3release_1115"
    AMP_CLINICAL_RELEASE_PATH=f"{AMP_RELEASE_PATH}/clinical"

    AMP_WGS_RELEASE_PATH="gs://amp-pd-genomics/releases/2022_v3release_1115/wgs-WB-DWGS"
    AMP_WGS_RELEASE_PLINK_PATH=os.path.join(AMP_WGS_RELEASE_PATH, "plink")
    AMP_WGS_RELEASE_PLINK_PFILES=os.path.join(AMP_WGS_RELEASE_PLINK_PATH, "pfiles")

    s=settings.get_value()
    s.ancestry=args.ancestry
    s.vcfToTsv=args.vcf_to_tsv
    s.homeDir=os.path.expanduser("~")   
    s.toolDir=os.path.join(s.homeDir, "data", "gch1", "tools")
    s.dataDir=os.path.join(s.homeDir, "data", "gch1", "data")    
    if args.ancestry!=None:
        s.billingProjectID=BILLING_PROJECT_ID
        s.workspaceNamespace=WORKSPACE_NAMESPACE
        s.workspaceName=WORKSPACE_NAME
        s.workspaceBucket=WORKSPACE_BUCKET
        s.dataDir=os.path.join(s.homeDir, "data", "gch1", "data", s.ancestry)    
    s.gp2ReleasePath=GP2_RELEASE_PATH
    s.gp2ClinicalReleasePath=GP2_CLINICAL_RELEASE_PATH
    s.gp2RawGenoPath=GP2_RAW_GENO_PATH
    s.gp2ImputedGenoPath=GP2_IMPUTED_GENO_PATH
    s.ampReleasePath=AMP_RELEASE_PATH
    s.ampClinicalReleasePath=AMP_CLINICAL_RELEASE_PATH
    s.ampWGSReleasePath=AMP_WGS_RELEASE_PATH
    s.ampWGSReleasePlinkPath=AMP_WGS_RELEASE_PLINK_PATH
    s.ampWGSReleasePlinkFiles=AMP_WGS_RELEASE_PLINK_PFILES 
    s.totalMemSizeGB=2 ** (int((psutil.virtual_memory().total / (1024 ** 3))).bit_length() - 1)
    settings.set_value(s)

    retag_exome = "AF_exome:=INFO/AF,INFO/AF_eas_exome:=INFO/AF_eas,INFO/AF_eas_exome:=INFO/AF_eas,INFO/AF_sas_exome:=INFO/AF_sas,INFO/AF_mid_exome:=INFO/AF_mid,INFO/AF_afr_exome:=INFO/AF_afr,INFO/AF_amr_exome:=INFO/AF_amr,INFO/AF_nfe_exome:=INFO/AF_nfe,INFO/AF_fin_exome:=INFO/AF_fin,INFO/AF_asj_exome:=INFO/AF_asj,INFO/AC_exome:=INFO/AC,INFO/AC_eas_exome:=INFO/AC_eas,INFO/AC_sas_exome:=INFO/AC_sas,INFO/AC_mid_exome:=INFO/AC_mid,INFO/AC_afr_exome:=INFO/AC_afr,INFO/AC_amr_exome:=INFO/AC_amr,INFO/AC_nfe_exome:=INFO/AC_nfe,INFO/AC_fin_exome:=INFO/AC_fin,INFO/AC_asj_exome:=INFO/AC_asj,INFO/AN_exome:=INFO/AN,INFO/nhomalt_exome:=INFO/nhomalt"
    retag_genome = "AF_genome:=INFO/AF,INFO/AF_eas_genome:=INFO/AF_eas,INFO/AF_eas_genome:=INFO/AF_eas,INFO/AF_sas_genome:=INFO/AF_sas,INFO/AF_mid_genome:=INFO/AF_mid,INFO/AF_afr_genome:=INFO/AF_afr,INFO/AF_amr_genome:=INFO/AF_amr,INFO/AF_nfe_genome:=INFO/AF_nfe,INFO/AF_fin_genome:=INFO/AF_fin,INFO/AF_asj_genome:=INFO/AF_asj,INFO/AC_genome:=INFO/AC,INFO/AC_eas_genome:=INFO/AC_eas,INFO/AC_sas_genome:=INFO/AC_sas,INFO/AC_mid_genome:=INFO/AC_mid,INFO/AC_afr_genome:=INFO/AC_afr,INFO/AC_amr_genome:=INFO/AC_amr,INFO/AC_nfe_genome:=INFO/AC_nfe,INFO/AC_fin_genome:=INFO/AC_fin,INFO/AC_asj_genome:=INFO/AC_asj,INFO/AN_genome:=INFO/AN,INFO/nhomalt_genome:=INFO/nhomalt"
    tsv_fields = 'CHROM POS REF ALT QUAL FILTER "ANN[*].GENE" "LOF[*].GENE" "NMD[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].BIOTYPE" "ANN[*].RANK" "AF_exome" "AF_eas_exome" "AF_eas_exome" "AF_sas_exome" "AF_mid_exome" "AF_afr_exome" "AF_amr_exome" "AF_nfe_exome" "AF_fin_exome" "AF_asj_exome" "AC_exome" "AC_eas_exome" "AC_sas_exome" "AC_mid_exome" "AC_afr_exome" "AC_amr_exome" "AC_nfe_exome" "AC_fin_exome" "AC_asj_exome" "AN_exome" "nhomalt_exome" "AF_genome" "AF_eas_genome" "AF_eas_genome" "AF_sas_genome" "AF_mid_genome" "AF_afr_genome" "AF_amr_genome" "AF_nfe_genome" "AF_fin_genome" "AF_asj_genome" "AC_genome" "AC_eas_genome" "AC_sas_genome" "AC_mid_genome" "AC_afr_genome" "AC_amr_genome" "AC_nfe_genome" "AC_fin_genome" "AC_asj_genome" "AN_genome" "nhomalt_genome" "GEN[*].GT"'
    #reheader = [("CHROM","#CHROM"), ("ANN\[\*\]\.GENE","Gene"), ("ANN\[\*\]\.FEATUREID","Transcript"), ("ANN\[\*\]\.HGVS_C","CDS"), ("ANN\[\*\]\.HGVS_P","AA"), ("ANN\[\*\]\.RANK","Exon"), ("ANN\[\*\].EFFECT", "Effect"), ("ANN\[\*\]\.BIOTYPE","Biotype"), ("LOF\[\*\]\.GENE","LOF"), ("NMD\[\*\]\.GENE","NMD"), ("dbNSFP_", "")]
    reheader_py = [("CHROM","#CHROM"), ("ANN[*].GENE","Gene"), ("ANN[*].FEATUREID","Transcript"), ("ANN[*].HGVS_C","CDS"), ("ANN[*].HGVS_P","AA"), ("ANN[*].RANK","Exon"), ("ANN[*].EFFECT", "Effect"), ("ANN[*].BIOTYPE","Biotype"), ("LOF[*].GENE","LOF"), ("NMD[*].GENE","NMD"), ("dbNSFP_", "")]
    headers = ["#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "Gene", "LOF", "NMD", "Transcript", "Effect", "CDS", "AA", "Biotype", "Exon", "AF_exome", "AF_eas_exome", "AF_eas_exome.1", "AF_sas_exome", "AF_mid_exome", "AF_afr_exome", "AF_amr_exome", "AF_nfe_exome", "AF_fin_exome", "AF_asj_exome", "AC_exome", "AC_eas_exome", "AC_sas_exome", "AC_mid_exome", "AC_afr_exome", "AC_amr_exome", "AC_nfe_exome", "AC_fin_exome", "AC_asj_exome", "AN_exome", "nhomalt_exome", "AF_genome", "AF_eas_genome", "AF_eas_genome.1", "AF_sas_genome", "AF_mid_genome", "AF_afr_genome", "AF_amr_genome", "AF_nfe_genome", "AF_fin_genome", "AF_asj_genome", "AC_genome", "AC_eas_genome", "AC_sas_genome", "AC_mid_genome", "AC_afr_genome", "AC_amr_genome", "AC_nfe_genome", "AC_fin_genome", "AC_asj_genome", "AN_genome", "nhomalt_genome"]

    if s.ancestry != None:
        lg.info(f"Ancestry: {ANCESTRY}")
        print("")
        lg.info("Billing and Workspace")
        lg.info(f"Workspace Namespace: {WORKSPACE_NAMESPACE}")
        lg.info(f"Workspace Name: {WORKSPACE_NAME}")
        lg.info(f"Billing Project: {BILLING_PROJECT_ID}")
        lg.info(f"Workspace Bucket, where you can upload and download data: {WORKSPACE_BUCKET}")
        print("")
        lg.info("GP2 v6.0")
        lg.info(f"Path to GP2 v6.0 Clinical Data: {GP2_CLINICAL_RELEASE_PATH}")
        lg.info(f"Path to GP2 v6.0 Raw Genotype Data: {GP2_RAW_GENO_PATH}")
        lg.info(f"Path to GP2 v6.0 Imputed Genotype Data: {GP2_IMPUTED_GENO_PATH}")
        print("")
        lg.info("AMP-PD v3.0")
        lg.info(f"Path to AMP-PD v3.0 Clinical Data: {AMP_CLINICAL_RELEASE_PATH}")
        lg.info(f"Path to AMP-PD v3.0 WGS Data: {AMP_WGS_RELEASE_PLINK_PATH}")
        lg.info(f"Path to AMP-PD v3.0 WGS Data: {AMP_WGS_RELEASE_PLINK_PFILES}")
        print("")

        install_plink2()
        install_snpeff()
        copy_data()
        extract_gch1()  
    
        result = annotate() | Pipe(gnomad_exome) | Pipe2(gnomad_filter, 0.05) | Pipe2(rename_tags, retag_exome) | \
        Pipe(gnomad_genome) | Pipe2(gnomad_filter, 0.05) | Pipe2(rename_tags, retag_genome)

        lg.info(f"Run for {s.ancestry} was completed.")

    elif s.vcfToTsv == True:        
        lg.info(f"Merging VCFs...")
        result = merge_vcf()
        lg.info(f"Extracts sample names from merged VCF...")
        samples = subprocess.check_output(f'zcat "{result}" | grep "#C" | cut -f 10-', shell=True).decode()
        samples = samples.replace("\n","")        
        result = result | Pipe2(to_tsv, tsv_fields) | Pipe2(reheader_tsv_py, reheader_py+[("GEN[*].GT", samples)])

        lg.info(f"Selecting samples contains alt. het. or alt. homo...")
        output = os.path.join(s.dataDir, "selected.tsv")
        df = pd.read_csv(result, delimiter="\t")
        gt_to_keep = ["0/1", "1/1", "1/0", "0|1", "1|1", "1|0"]
        rows = (df.map(lambda x: str(x) in gt_to_keep)).any(axis=1)
        cols = (df.map(lambda x: str(x) in gt_to_keep)).any(axis=0)
        for colHeader in headers:
            cols[colHeader] = True
        filtered = df.loc[rows]
        filtered = filtered.loc[:, cols]
        filtered.to_csv(output, index=False, sep="\t")
        lg.info(f"Results were saved to: {output}")

        #result = annotate() | Pipe(gnomad_exome) | Pipe2(gnomad_filter, 0.05) | Pipe2(rename_tags, retag_exome) | \
        #Pipe(gnomad_genome) | Pipe2(gnomad_filter, 0.05) | Pipe2(rename_tags, retag_genome) | \
        #Pipe2(to_tsv, tsv_fields) | Pipe2(reheader_tsv, reheader)
        #shell_do(f'gsutil -mu {s.billingProjectID} cp -r {s.dataDir} {s.workspaceBucket}')

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="GCH1 Project on GP2")
    
    group = parser.add_mutually_exclusive_group(required=True)    
    group.add_argument("-tsv", "--vcf-to-tsv", help="Converts ancestry VCFs to tsv", action="store_true")
    group.add_argument("-e", "--ancestry", type=str, help="Enter one of the following: AAC, AFR, AJ, AMR, CAH, CAS, EAS, EUR, FIN, MDE, or SAS.", choices=["AAC", "AFR", "AJ", "AMR", "CAH", "CAS", "EAS", "EUR", "FIN", "MDE", "SAS"])

    parser.add_argument("-id", "--billing-project-id", type=str, help="GP2 Terra Billing Project ID")
    parser.add_argument("-ns", "--workspace-namespace", type=str, help="GP2 Terra Namespace of Workspace")
    parser.add_argument("-ws", "--workspace-name", type=str, help="GP2 Terra Workspace Name")
    parser.add_argument("-gcs", "--workspace-bucket", type=str, help="GP2 Terra Workspace bucket (Google Cloud Storage)")
    args=parser.parse_args()
    
    main(args)