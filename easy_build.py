#!/usr/bin/env python
import os, sys
import subprocess
from timeit import time

ibamr_config_options = ''
petsc_config_options = ''
SKIP_PETSC           = False
QUIET                = False

if 'LC_LOCAL' in os.environ and os.environ['LC_LOCAL'] != '' and os.environ['LC_LOCAL'] != 'en_US' and os.environ['LC_LOCAL']!= 'en_US.UTF-8': os.environ['LC_LOCAL'] = 'en_US.UTF-8'
if 'LANG' in os.environ and os.environ['LANG'] != '' and os.environ['LANG'] != 'en_US' and os.environ['LANG'] != 'en_US.UTF-8': os.environ['LANG'] = 'en_US.UTF-8'

if not hasattr(sys, 'version_info') or not sys.version_info[0] == 2 or not sys.version_info[1] >= 4:
  print '*** You must have Python2 version 2.4 or higher to run ./configure        *****'
  print '*          Python is easy to install for end users or sys-admin.              *'
  print '*                  http://www.python.org/download/                            *'
  print '*                                                                             *'
  print '*           You CANNOT configure IBAMR without Python                         *'
  print '*******************************************************************************'
  sys.exit(4)

if 'IBAMR_DIR' not in os.environ:
    print '****************************************************************************'
    print '*               $IBAMR_DIR must be defined.                                *'
    print '*    Navigate to root directory of IBAMR and set IBAMR_DIR = $PWD          *'
    print '****************************************************************************'
    sys.exit(4)

if os.environ['PWD'] != os.environ['IBAMR_DIR']:
    print '****************************************************************************'
    print '*              `easy_build.py` must be invoked from \$IBAMR_DIR                     *'
    print '****************************************************************************'
    sys.exit(4)

if 'IBAMR_ARCH' not in os.environ:
    print '****************************************************************************'
    print '*                    No \$IBAMR_ARCH specified.                            *'
    while ('IBAMR_ARCH' not in os.environ):
        debug_response = raw_input("Enter 'yes' for a DEBUG BUILD  >  ")
        if debug_response.lower().startswith("yes"):
            os.environ['IBAMR_ARCH'] = 'debug'
            os.environ['DEBUG_OPTION'] = '--with-debugging'
        else:
            opt_response = raw_input("Enter 'yes' for and OPTIMIZED BUILD (for production runs only)  >  ")
            if opt_response.lower().startswith("yes"):
                os.environ['IBAMR_ARCH'] = 'opt'
                os.environ['DEBUG_OPTION'] = ''
            else:
                print "Try again. Please enter 'yes' for either DEBUG BUILD or OPTIMIZED BUILD"
    print '****************************************************************************'

IBAMR_DIR = os.environ['IBAMR_DIR']
IBAMR_ARCH = os.environ['IBAMR_ARCH']

print '****************************************************************************'
print '  BUILDING IBAMR and its dependencies in %s/%s' %(IBAMR_DIR, IBAMR_ARCH)
print '****************************************************************************'

for item in sys.argv:
    if 'debug' in item:
        print "Found option %s, building libraries with DEBUGGING ENABLED" %item
        petsc_config_options += ' --with-debugging=1 '
        ibamr_config_options += ' --with-debugging=1 '
        sys.argv.remove(item)
    if 'opt' in item:
        print "Found option %s, building libraries WITHOUT debugging flags (OPTIMIZED BUILD)" %item
        petsc_config_options += ' --with-debugging=0 '
        ibamr_config_options += ' --with-debugging=0 '
        sys.argv.remove(item)
    if 'CC' in item:
        petsc_config_options += ' ' + item + ' '
        sys.argv.remove(item)
    if 'CXX' in item:
        petsc_config_options += ' ' + item + ' '
        sys.argv.remove(item)
    if 'FC' in item:
        petsc_config_options += ' ' + item + ' '
        sys.argv.remove(item)
    if 'skip-petsc' in item:
        SKIP_PETSC = True
    if '--quiet' in item:
        QUIET = True

additional_options =' '.join(sys.argv[1:])
petsc_config_options += additional_options

if 'debug' not in petsc_config_options:
    print '****************************************************************************'
    print '*                     Is this a DEBUG build?                               *'
    print '*        (only answer "no" for production or "optimized" builds)           *'
    print '* To avoid this prompt in the future, invoke this script with the options: *'
    print '*          ./easy_build.py --with-debugging for DEBUG BUILD                *'
    print '*          ./easy_build.py --optimized      for OPTIMIZED BUILD            *'
    print '*                                                                          *'
    debug_response = raw_input("Is this a DEBUG build? (yes/no) >  ")
    if debug_response.lower().startswith("yes"):
        petsc_config_options += ' --with-debugging=1 '
        ibamr_config_options += ' --with-debugging=1 '
    print '****************************************************************************'

print '****************************************************************************'
if '--CC=' not in petsc_config_options:
    petsc_config_options +=' --CC=' +raw_input("Please enter name of C compiler        >  ") + ' '
    print ' In the future, avoid this prompt by passing --CC=$NAME_OF_C_COMPILER to easy_build.py'
if '--CXX=' not in petsc_config_options:
    petsc_config_options +=' --CXX=' +raw_input("Please enter name of C++ compiler        >  ") + ' '
    print ' In the future, avoid this prompt by passing --CXX=$NAME_OF_CXX_COMPILER to easy_build.py'
if '--FC=' not in petsc_config_options:
    petsc_config_options +=' --FC='+raw_input("Please enter name of Fortran compiler        >  ") + ' '
    print ' In the future, avoid this prompt by passing --FC=$NAME_OF_FORTRAN_COMPILER to easy_build.py'
print '****************************************************************************'

def execute(command, quiet=QUIET):
    '''https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running'''
    if quiet:
        return subprocess.check_output(command, shell=True)
    else:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output  = ''
        for line in iter(process.stdout.readline, ''):
            print line
            output += line

        process.wait()
        exitCode = process.returncode
        return output

get_petsc       = 'mkdir -p $PETSC_DIR && git clone -b maint https://bitbucket.org/petsc/petsc $PETSC_DIR'
fetch_petsc     = 'mkdir -p $PETSC_DIR && cd $PETSC_DIR && git fetch origin && git pull origin maint'
configure_petsc = 'cd $PETSC_DIR && ./configure --download-openmpi --download-fblaslapack ' + petsc_config_options
make_petsc      = 'cd $PETSC_DIR && make'
# NOTE: libmesh breaks when HDF5 is built by PETSc, but works when built by IBAMR.
configure_ibamr = './configure.new --IBAMR_ARCH=$IBAMR_ARCH --download-muparser --download-silo '
configure_ibamr += ' --download-hdf5 --download-libmesh --download-boost --boost-headers-only '
configure_ibamr += ' --download-samrai --download-fblaslapack --with-mpi-dir=$PETSC_DIR/$PETSC_ARCH '
configure_ibamr += ' --download-eigen ' + ibamr_config_options
make_ibamr      = 'make all'

if not SKIP_PETSC:
    os.environ['PETSC_DIR'] = os.path.join(os.path.join(IBAMR_DIR, IBAMR_ARCH), 'PETSC')
    os.environ['PETSC_ARCH'] = IBAMR_ARCH
    if os.path.exists(os.environ['PETSC_DIR']) and os.path.exists(os.path.join(os.environ['PETSC_DIR'], 'include/petsc.h')):
        execute(fetch_petsc)
    else:
        execute(get_petsc)

    execute(configure_petsc)
    execute(make_petsc)

if ('PETSC_DIR' not in os.environ):
    os.environ['PETSC_DIR'] = raw_input("PETSC_DIR not defined, please provide path to petsc: ")
if ('PETSC_ARCH' not in os.environ):
    os.environ['PETSC_ARCH'] = raw_input("PETSC_ARCH not defined, please provide name of PETSC_ARCH: ")

execute(configure_ibamr)
execute(make_ibamr)
