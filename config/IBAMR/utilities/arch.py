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
    if not hasattr(self, 'arch'):
      return ''
    desc = [self.Project+':']
    desc.append('  '+self.PROJECT+'_ARCH: '+str(self.arch))
    return '\n'.join(desc)+'\n'

  def setupHelp(self, help):
    import nargs
    help.addArgument(self.Project, '-'+self.PROJECT+'_ARCH=<string>',     nargs.Arg(None, None, 'The configuration name'))
    help.addArgument(self.Project, '-with-'+self.project+'-arch=<string>',nargs.Arg(None, None, 'The configuration name'))
    return

  def configureArchitecture(self):
    '''Checks PROJECT_ARCH and sets if not set'''


    # Warn if PROJECT_ARCH doesnt match env variable
    if self.PROJECT+'_ARCH' in self.framework.argDB and self.PROJECT+'_ARCH' in os.environ and self.framework.argDB[self.PROJECT+'_ARCH'] != os.environ[self.PROJECT+'_ARCH']:
      self.logPrintBox('''\
Warning: %s_ARCH from environment does not match command-line or name of script.
Warning: Using from command-line or name of script: %s, ignoring environment: %s''' % (self.PROJECT, str(self.framework.argDB[self.PROJECT+'_ARCH']), str(os.environ[self.PROJECT+'_ARCH'])))
    if 'with-'+self.project+'-arch' in self.framework.argDB:
      self.arch = self.framework.argDB['with-'+self.project+'-arch']
    elif self.PROJECT+'_ARCH' in self.framework.argDB:
      self.arch = self.framework.argDB[self.PROJECT+'_ARCH']
    else:
      if self.PROJECT+'_ARCH' in os.environ:
        if not len(os.environ[self.PROJECT+'_ARCH']):
          raise RuntimeError(self.PROJECT+'_ARCH is the empty string in your environment. It must either be a valid string, or not be defined in the environment at all.')
        self.arch = os.environ[self.PROJECT+'_ARCH']
      else:
        import sys
        self.arch = 'arch-' + sys.platform.replace('cygwin','mswin')
        # use opt/debug, c/c++ tags.s
        self.arch+= '-'+self.framework.argDB['with-clanguage'].lower().replace('+','x')
        if self.framework.argDB['with-debugging']:
          self.arch += '-debug'
        else:
          self.arch += '-opt'
    if self.arch.find('/') >= 0 or self.arch.find('\\') >= 0:
      raise RuntimeError(self.PROJECT+'_ARCH should not contain path characters, but you have specified: '+str(self.arch))
    self.archBase = re.sub(r'^(\w+)[-_]?.*$', r'\1', self.arch)
    self.addDefine('ARCH', '"'+self.arch+'"')
    return

  def configure(self):
    self.executeTest(self.configureArchitecture)
    # required by top-level configure.py
    self.framework.arch = self.arch
    return
