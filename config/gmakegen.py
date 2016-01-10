#!/usr/bin/env python
import os, sys

sys.path.insert(0, os.path.join(os.environ['IBAMR_DIR'], 'config'))
sys.path.insert(0, os.path.join(os.environ['PETSC_DIR'], 'config', 'BuildSystem'))

import script

class ConfigReader(script.Script):
  def __init__(self):
    import RDict
    import os

    self.archDir = os.path.join(os.environ['IBAMR_DIR'], os.environ['IBAMR_ARCH'])
    argDB = RDict.RDict(None, None, 0, 0)
    argDB.saveFilename = os.path.join(self.archDir, 'conf', 'RDict.db')
    argDB.load()
    script.Script.__init__(self, argDB = argDB)
    return

  def run(self):
    self.setup()
    self.framework = self.loadConfigure()
    return

class GMakefileGenerator(ConfigReader):
  def __init__(self):
    ConfigReader.__init__(self)
    return

  def run(self):
    ConfigReader.run(self)
    with file(os.path.join(self.archDir, 'gmakefile'), 'w') as f:
      self.framework.outputMakeMacro(f, 'PYTHON', sys.executable)

      self.framework.outputMakeMacro(f, 'PETSC_DIR',  os.environ['PETSC_DIR'])
      self.framework.outputMakeMacro(f, 'PETSC_ARCH', os.environ['PETSC_ARCH'])
      f.write('include ${PETSC_DIR}/lib/petsc/conf/variables\n\n')

      conf = self.framework.require('IBAMR.Configure', None)
      self.framework.outputMakeMacro(f, 'BOOST_INCLUDE', '-I'+os.path.join(conf.projectdir.dir, 'ibtk', 'contrib', 'boost'))
      self.framework.outputMakeMacro(f, 'IBAMR_DIR',  conf.projectdir.dir)
      self.framework.outputMakeMacro(f, 'IBAMR_ARCH', conf.arch.arch)
      inc = []
      lib = []
      for packageName in ['samrai', 'hdf5', 'eigen', 'silo', 'muparser']:
        package = getattr(conf, packageName)
        NAME    = package.PACKAGE.replace('-','_')
        if hasattr(package, 'include'):
          self.framework.outputMakeMacro(f, NAME+'_INCLUDE', conf.headers.toStringNoDupes(package.include))
          inc.append('${'+package.PACKAGE+'_INCLUDE}')
        if hasattr(package, 'lib'):
          self.framework.outputMakeMacro(f, NAME+'_LIB',     conf.libraries.toStringNoDupes(package.lib))
          lib.append('${'+package.PACKAGE+'_LIB}')
      self.framework.outputMakeMacro(f, 'IBAMR_INCLUDE', '-I${IBAMR_DIR}/include -I${IBAMR_DIR}/ibtk/include '+' '.join(inc))
      self.framework.outputMakeMacro(f, 'IBAMR_LIB',  '-L${IBAMR_DIR}/${IBAMR_ARCH}/lib -libamr '+' '.join(lib))
      f.write('CFLAGS += ${IBAMR_INCLUDE} ${BOOST_INCLUDE} -DNDIM=2\n')
      f.write('include ${IBAMR_DIR}/${IBAMR_ARCH}/lib/petsc/conf/ibamrvariables\n\ninclude ${IBAMR_DIR}/base.mk\n')
    return

if __name__ == '__main__':
  GMakefileGenerator().run()
