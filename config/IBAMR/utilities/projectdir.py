import config.base
import os
import re

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    self.Project      = self.framework.Project
    self.project      = self.Project.lower()
    self.PROJECT      = self.Project.upper()
    self.headerPrefix = self.PROJECT
    self.substPrefix  = self.PROJECT
    return

  def __str1__(self):
    if hasattr(self, 'dir'):
      return '  '+self.PROJECT+'_DIR: '+str(self.dir)+'\n'
    return ''

  def setupHelp(self, help):
    import nargs
    help.addArgument(self.Project, '-'+self.PROJECT+'_DIR=<root-dir>', nargs.Arg(None, None, 'The root directory of the '+self.Project+' installation'))
    return

  def configureDirectories(self):
    '''Checks PROJECT_DIR and sets if not set'''
    pdir = self.PROJECT+'_DIR'
    if pdir in self.framework.argDB:
      self.dir = self.framework.argDB[pdir]
      if self.dir == 'pwd':
        raise RuntimeError('You have set -'+pdir+'=pwd, you need to use back quotes around the pwd\n  like -'+pdir+'=`pwd`')
      if not os.path.isdir(self.dir):
        raise RuntimeError('The value you set with -'+pdir+'='+self.dir+' is not a directory')
    elif pdir in os.environ:
      self.dir = os.environ[pdir]
      if self.dir == 'pwd':
        raise RuntimeError ('''
The environmental variable %s is set incorrectly. Please use the following: [notice backquotes]
  For sh/bash  : %s=`pwd`; export %s
  for csh/tcsh : setenv %s `pwd`''' % (pdir, pdir, pdir, pdir))
      elif not os.path.isdir(self.dir):
        raise RuntimeError('The environmental variable '+pdir+' '+self.dir+' is not a directory')
    else:
      self.dir = os.getcwd()
    if not os.path.realpath(self.dir) == os.path.realpath(os.getcwd()):
      raise RuntimeError('The environmental variable '+pdir+' '+self.dir+' MUST be the current directory '+os.getcwd())
    self.version  = 'Unknown'
    versionHeader = os.path.join(self.dir, 'include', self.project+'version.h')
    versionInfo = []
    if os.path.exists(versionHeader):
      f = file(versionHeader)
      for line in f:
        if line.find('define '+self.PROJECT+'_VERSION') >= 0:
          versionInfo.append(line[:-1])
      f.close()
    else:
      raise RuntimeError('Invalid '+self.Project+' directory '+str(self.dir)+'. Could not locate '+versionHeader)
    self.version = '.'.join([line.split(' ')[-1] for line in versionInfo[1:4]])
    self.logPrint('Version Information:')
    for line in versionInfo:
      self.logPrint(line)
    self.addMakeMacro('DIR', self.dir)
    self.framework.argDB['search-dirs'].append(os.path.join(self.dir, 'bin', 'win32fe'))

    return

  def configure(self):
    self.executeTest(self.configureDirectories)
    return
