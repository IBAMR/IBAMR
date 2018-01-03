import config.base
import os

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    return

  def setupHelp(self, help):
    import nargs
    help.addArgument('External Packages', '-with-clean', nargs.ArgBool(None, 0, 'Remove all external packages and build files before reinstalling'))
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

  def setExternalPackagesDir(self):
    if self.framework.externalPackagesDir is None:
      self.dir = os.path.join(os.path.abspath(os.path.join(self.arch)), 'externalpackages')
    else:
      self.dir = os.path.join(self.framework.externalPackagesDir, self.arch)
    return

  def cleanExternalpackagesDir(self):
    if self.framework.argDB['with-clean'] and os.path.isdir(self.dir):
      import shutil
      self.logPrintBox('Warning: "with-clean" is specified. Removing all externalpackage files from '+ self.dir)
      shutil.rmtree(self.dir)
    return

  def configure(self):
    self.executeTest(self.setExternalPackagesDir)
    self.executeTest(self.cleanExternalpackagesDir)
    return
