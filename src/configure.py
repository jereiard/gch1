class Settings:
    def __init__(self):
        self.value = None
        self.HOME_DIR = None
        self.PROJECT_DIR = None
        self.PROJECT_DATA_DIR = None
        self.USER_DATA_DIR = None
        self.TOOL_DIR = None
        self.DATA_ROOT = None
        self.CLINICAL_DATA_DIR = None
        self.BILLING_PROJECT_ID = None
        self.WORKSPACE_NAMESPACE = None
        self.WORKSPACE_NAME = None
        self.WORKSPACE_BUCKET = None
        self.GP2_WGS_VCF_PATH = None
        self.GP2_RELEASE_PATH = None
        self.GP2_CLINICAL_RELEASE_PATH = None
        self.GP2_RAW_GENO_PATH = None
        self.GP2_IMPUTED_GENO_PATH = None
        self.AMP_RELEASE_PATH = None
        self.AMP_CLINICAL_RELEASE_PATH = None
        self.AMP_WGS_RELEASE_PATH = None
        self.AMP_WGS_RELEASE_PLINK_PATH = None
        self.AMP_WGS_RELEASE_PLINK_PFILES = None
        self.TOTAL_MEM_SIZE = None
        self.GP2_BUCKET = None
        self.AMP_BUCKET = None
    
__all__ = ['Settings']