"""
If HEAD is an exact tag match, this is an official release, so
return no suffix. Otherwise it's a development build so return
suffix '+shortsha'.
"""

import subprocess, sys                                                                                                           
  
try:                                                                                                                             
    tag = subprocess.check_output(
        ['git', 'describe', '--exact-match', '--tags', 'HEAD'],                                                                  
        stderr=subprocess.DEVNULL
    ).decode().strip()                                                                                                           
    print('', end='')                                                                                                            
except subprocess.CalledProcessError:                                                                                            
    sha = subprocess.check_output(                                                                                               
        ['git', 'rev-parse', '--short', 'HEAD']
    ).decode().strip()                                                                                                           
    print(f'+{sha}', end='')