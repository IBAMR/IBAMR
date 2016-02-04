#!/usr/bin/env python

# clang-format is in the user's bin directory; however apps don't inherit the same paths that
# shells do. so we need to manually add it. 
import os

default_path = os.environ['PATH']
os.environ['PATH'] =  os.path.expanduser("~/bin") + ":/usr/local/bin:" + default_path 
child_env = os.environ.copy()

# when run from within an app (e.g. SourceTree), subprocess will throw an exception because
# of the error from git clang-format. we'll ignore that to carry on
import subprocess
p = subprocess.Popen(["git", "clang-format", "--diff"], stdout=subprocess.PIPE, env=child_env)
output, err = p.communicate()

if output not in ['no modified files to format\n', 'clang-format did not modify any files\n']:
    print "Run git clang-format, then commit.\n"        
    exit(1)
else:
    exit(0)
