import config.base
import os

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    return

  def getArch(self):
    '''The architecture identifier'''
    if hasattr(self, 'archProvider'):
      if hasattr(self.archProvider, 'arch'):
        return self.archProvider.arch
    return self._arch
  def setArch(self, arch):
    '''The architecture identifier'''
    self._arch = arch
    return
  arch = property(getArch, setArch, doc = 'The architecture identifier')

  def printSudoPasswordMessage(self, needsudo = 1):
    '''Prints a message that sudo password will be needed for installs of packages'''
    '''Packages like sowing and make that are never installed in an sudo location would pass 0 for needsudo'''
    if needsudo and self.installSudoMessage:
      self.logPrintBox(self.installSudoMessage)
      self.installSudoMessage = ''
    return

  def setInstallDir(self):
    ''' setup installDir to either prefix or if that is not set to the chosen ARCH'''
    self.installSudo        = ''
    self.installSudoMessage = ''
    if self.framework.argDB['prefix']:
      self.dir = self.framework.argDB['prefix']
      try:
        os.makedirs(os.path.join(self.dir, 'BuildSystemTestDirectory'))
        os.rmdir(os.path.join(self.dir, 'BuildSystemTestDirectory'))
      except:
        self.installSudoMessage = 'You do not have write permissions to the --prefix directory '+self.dir+'\nYou will be prompted for the sudo password for any external package installs'
        self.installSudo = 'sudo '
    else:
      self.dir = os.path.abspath(os.path.join(self.arch))
    self.confDir = os.path.abspath(os.path.join(self.arch))
    return

  def configureInstallDir(self):
    '''Makes  installDir subdirectories if it does not exist for both prefix install location and work install location'''
    dir = os.path.abspath(os.path.join(self.arch))
    if not os.path.exists(dir):
      os.makedirs(dir)
    for i in ['include','lib','bin','conf']:
      newdir = os.path.join(dir,i)
      if not os.path.exists(newdir):
        os.makedirs(newdir)
    if os.path.isfile(self.framework.argDB.saveFilename):
      os.remove(self.framework.argDB.saveFilename)
    confdir = os.path.join(dir,'conf')
    self.framework.argDB.saveFilename = os.path.abspath(os.path.join(confdir, 'RDict.db'))
    self.framework.logPrint('Changed persistence directory to '+confdir)
    return

  def cleanInstallDir(self):
    if 'with-clean' in self.framework.argDB and self.framework.argDB['with-clean'] and os.path.isdir(self.dir):
      import shutil
      self.logPrintBox('Warning: "with-clean" is specified. Removing all build files from '+ self.dir)
      shutil.rmtree(self.dir)
    return

  def saveReconfigure(self):
    self.reconfigure_file = os.path.join(self.dir,'conf','reconfigure-'+self.arch+'.py')
    self.save_reconfigure_file = None
    if 'with-clean' in self.framework.argDB and self.framework.argDB['with-clean'] and os.path.exists(self.reconfigure_file):
      self.save_reconfigure_file = '.save.reconfigure-'+self.arch+'.py'
      try:
        if os.path.exists(self.save_reconfigure_file): os.unlink(self.save_reconfigure_file)
        os.rename(self.reconfigure_file,self.save_reconfigure_file)
      except Exception, e:
        self.save_reconfigure_file = None
        self.framework.logPrint('error in saveReconfigure(): '+ str(e))
    return

  def restoreReconfigure(self):
    if 'with-clean' in self.framework.argDB and self.framework.argDB['with-clean'] and self.save_reconfigure_file:
      try:
        os.rename(self.save_reconfigure_file,self.reconfigure_file)
      except Exception, e:
        self.framework.logPrint('error in restoreReconfigure(): '+ str(e))
    return

  def configure(self):
    self.executeTest(self.setInstallDir)
    self.executeTest(self.saveReconfigure)
    self.executeTest(self.cleanInstallDir)
    self.executeTest(self.configureInstallDir)
    self.executeTest(self.restoreReconfigure)
    return
