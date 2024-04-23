import os, sys
import subprocess
import pandas as pd
from io import StringIO
import urllib.parse
import shlex

# Utility routine for printing a shell command before executing it
def shell_do(command):
    print(f'Executing: {command}', file=sys.stderr)
    result = subprocess.run(shlex.split(command), capture_output=False, text=True, shell=True)
    return result.returncode 
    
def shell_return(command):
    print(f'Executing: {command}', file=sys.stderr)    
    #result = subprocess.run([command + args], capture_output=False, text=True, shell=True)
    result = subprocess.run(shlex.split(command), capture_output=True, text=True, shell=True)
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