import os
import argparse
from pck.gchutil import *
from pck.log import logger as lg

def install_packages():
    homeDir = os.path.expanduser("~")
    toolDir = os.path.join(homeDir, "tools")
    plinkDir = os.path.join(toolDir, "plink2")
    workDir = os.getcwd()

    try:
        if not os.path.exists(toolDir):
            os.makedirs(toolDir)
        if os.path.exists(plinkDir):
            lg.info(f"Plink2 is already installed in {plinkDir}")
        else:
            os.makedirs(plinkDir)
            os.chdir(plinkDir)

            status = shell_do("wget -N https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240105.zip")
            if status != 0: 
                raise RuntimeError("Cannot download plink2.")
            
            status = shell_do("unzip -o plink2_linux_avx2_20240105.zip")
            if status != 0: 
                raise RuntimeError("Cannot install plink2.")

            lg.info(f"plink2 downloaded and unzipped in {plinkDir}")

    except Exception as ex:
        lg.critical(ex)
        exit(1)

    finally:
        os.chdir(workDir)

def main(args):
    ANCESTRY = args.ancestry
    BILLING_PROJECT_ID = args.billing_project_id if args.billing_project_id is not None else os.environ["GOOGLE_PROJECT"]
    WORKSPACE_NAMESPACE = args.workspace_namespace if args.workspace_namespace is not None else os.environ["WORKSPACE_NAMESPACE"]
    WORKSPACE_NAME = args.workspace_name if args.workspace_name is not None else os.environ["WORKSPACE_NAME"]
    WORKSPACE_BUCKET = args.workspace_bucket if args.workspace_bucket is not None else os.environ["WORKSPACE_BUCKET"]

    GP2_RELEASE_PATH = "gs://gp2tier2/release6_21122023"
    GP2_CLINICAL_RELEASE_PATH = f"{GP2_RELEASE_PATH}/clinical_data"
    GP2_RAW_GENO_PATH = f"{GP2_RELEASE_PATH}/raw_genotypes"
    GP2_IMPUTED_GENO_PATH = f"{GP2_RELEASE_PATH}/imputed_genotypes"

    AMP_RELEASE_PATH = "gs://amp-pd-data/releases/2022_v3release_1115"
    AMP_CLINICAL_RELEASE_PATH = f"{AMP_RELEASE_PATH}/clinical"

    AMP_WGS_RELEASE_PATH = "gs://amp-pd-genomics/releases/2022_v3release_1115/wgs-WB-DWGS"
    AMP_WGS_RELEASE_PLINK_PATH = os.path.join(AMP_WGS_RELEASE_PATH, "plink")
    AMP_WGS_RELEASE_PLINK_PFILES = os.path.join(AMP_WGS_RELEASE_PLINK_PATH, "pfiles")

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

    install_packages()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GCH1 Projcet on GP2")
    parser.add_argument("-e", "--ancestry", type=str, help="Enter one of the following: AAC, AFR, AJ, AMR, CAH, CAS, EAS, EUR, FIN, MDE, or SAS.", required=True)
    parser.add_argument("-id", "--billing-project-id", type=str, help="GP2 Terra Billing Project ID")
    parser.add_argument("-ns", "--workspace-namespace", type=str, help="GP2 Terra Namespace of Workspace")
    parser.add_argument("-ws", "--workspace-name", type=str, help="GP2 Terra Workspace Name")
    parser.add_argument("-gcs", "--workspace-bucket", type=str, help="GP2 Terra Workspace bucket (Google Cloud Storage)")
    args = parser.parse_args()
    
    main(args)