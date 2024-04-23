import os, sys
from pathlib import Path
import subprocess
import pandas as pd
from io import StringIO
import urllib.parse
import shlex

# Utility routine for printing a shell command before executing it
def shell_do(command):
    print(f'Executing: {command}', file=sys.stderr)
    cmd = shlex.split(command, posix=False)
    #result = subprocess.run(cmd, capture_output=False, text=True, shell=True)
    returncode = subprocess.call(cmd)
    return returncode 

def shell_do_redir_stdout(command, outputFile):
    print(f'Executing: {command}', file=sys.stderr)
    print(f'Output: {outputFile}', file=sys.stderr)
    cmd = shlex.split(command, posix=False)
    #result = subprocess.run(cmd, capture_output=False, text=True, shell=True)
    with open(outputFile, "w") as f:
        returncode = subprocess.call(cmd, stdout=f)
    return returncode 
    
def shell_return(command):
    print(f'Executing: {command}', file=sys.stderr)    
    cmd = shlex.split(command, posix=False)
    print(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    return result.stdout

# Utility routine for printing a query before executing it
def bq_query(query):
    print(f'Executing: {query}', file=sys.stderr)
    return pd.read_gbq(query, project_id=os.environ['GOOGLE_PROJECT'], dialect='standard')

# Utility routine for display a message and a link
def display_html_link(description, link_text, url):
    html = f'''
    <p>
    </p>
    <p>
    {description}
    <a target=_blank href="{url}">{link_text}</a>.
    </p>
    '''
    print(html)
    #display(HTML(html))

# Utility routines for reading files from Google Cloud Storage
def gcs_read_file(path):
    """Return the contents of a file in GCS"""
    result = subprocess.run(["gsutil" f"-u {os.environ['GOOGLE_PROJECT']}", f"cat {path}"], capture_output=False, text=True, shell=True)
    return result.stdout
    #contents = !gsutil -u {BILLING_PROJECT_ID} cat {path}
    #return '\n'.join(contents)
    
def gcs_read_csv(path, sep=None):
    """Return a DataFrame from the contents of a delimited file in GCS"""
    return pd.read_csv(StringIO(gcs_read_file(path)), sep=sep, engine='python')

# Utility routine for displaying a message and link to Cloud Console
def link_to_cloud_console_gcs(description, link_text, gcs_path):
    url = '{}?{}'.format(
        os.path.join('https://console.cloud.google.com/storage/browser',
                     gcs_path.replace("gs://","")),
        urllib.parse.urlencode({'userProject': os.environ['GOOGLE_PROJECT']}))
    display_html_link(description, link_text, url)

def append_id(filename, id, destdir=None):
    p = Path(filename)
    s = Path(p.stem)
    if destdir == None:
        if p.suffix == ".gz" and s.suffix ==".vcf":
            return os.path.join(
                p.parent,
                f"{s.stem}{'.' if not id == '' else ''}{id}{s.suffix}{p.suffix}"
            )
        else:
            return os.path.join(
                p.parent,
                f"{p.stem}{'.' if not id == '' else ''}{id}{p.suffix}"
            )
    else:
        if not os.path.exists(os.path.join(destdir)):
            os.mkdir(os.path.join(destdir))
        if p.suffix == ".gz" and s.suffix ==".vcf":
            return os.path.join(                
                destdir,
                f"{s.stem}{'.' if not id == '' else ''}{id}{s.suffix}{p.suffix}",
        )
        else:
            return os.path.join(
                destdir,
                f"{p.stem}{'.' if not id == '' else ''}{id}{p.suffix}",
        )

def is_null_or_small(path, size=78848):
    if not os.path.exists(path):
        return True
    if os.path.isfile(path) and os.stat(path).st_size > size: # 적당한 컷오프 77K?
        return False
    else:
        return True
    
def change_ext(filename, ext):
        p = Path(filename)
        s = Path(p.stem)
        if p.suffix == ".gz" and s.suffix == ".vcf":
            return "{0}{1}".format(Path.joinpath(p.parent, s.stem), ext)
        else:
            return "{0}{1}".format(Path.joinpath(p.parent, p.stem), ext)
        