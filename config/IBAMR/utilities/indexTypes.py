#!/usr/bin/env python
from __future__ import generators
import user
import config.base

class Configure(config.base.Configure):
  def __init__(self, framework):
    config.base.Configure.__init__(self, framework)
    self.headerPrefix = ''
    self.substPrefix  = ''
    return

  def __str1__(self):
    desc = []
    if hasattr(self, 'integerSize'):
      desc.append('  Integer size: ' + str(self.integerSize))
    return '\n'.join(desc)+'\n'

  def setupHelp(self, help):
    import nargs
    help.addArgument('IBAMR', '-with-64-bit-indices=<bool>',   nargs.ArgBool(None, 0, 'Use 64 bit integers (long long) for indexing in vectors and matrices'))
    return

  def setupDependencies(self, framework):
    config.base.Configure.setupDependencies(self, framework)
    self.setCompilers = framework.require('config.setCompilers', self)
    self.libraries    = framework.require('config.libraries', None)
    self.compilers    = framework.require('config.compilers', None)
    return

  def configureIndexSize(self):
    if self.framework.argDB['with-64-bit-indices']:
      self.integerSize = 64
      self.addDefine('USE_64BIT_INDICES', 1)
      if self.libraries.check('-lgcc_s.1', '__floatdidf'):
        self.compilers.LIBS += ' '+self.libraries.getLibArgument('-lgcc_s.1')
      self.addMakeMacro('IBAMR_INDEX_SIZE', '64')
    else:
      self.integerSize = 32
      self.addMakeMacro('IBAMR_INDEX_SIZE', '32')
    return

  def configure(self):
    self.executeTest(self.configureIndexSize)
    return
