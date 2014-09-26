import config.package

class Configure(config.package.GNUPackage):
  def __init__(self, framework):
    config.package.GNUPackage.__init__(self, framework)
    self.download  = ['https://computation-rnd.llnl.gov/SAMRAI/download/SAMRAI-v2.4.4.tar.gz']
    self.functions = []
    self.includes  = ['muParser.h']
    self.liblist   = [['libmuparser.a']]
    self.pkgname   = 'muparser-2.2.3'
    self.cxx       = 1
    return

  def setupDependencies(self, framework):
    config.package.GNUPackage.setupDependencies(self, framework)
    self.setCompilers = framework.require('config.setCompilers', self)
    self.mpi  = framework.require('config.packages.MPI', self)
    self.hdf5 = framework.require('config.packages.hdf5', self)
    self.silo = framework.require('IBAMR.packages.silo', self)
    self.deps = [self.mpi, self.hdf5, self.silo]
    return

  def formGNUConfigureArgs(self):
    args = config.package.GNUPackage.formGNUConfigureArgs(self)
    args.append('--with-MPICC='+self.setCompilers.CC)
    args.append('--with-hdf5='+self.hdf5.installDir)
    args.append('--with-silo='+self.silo.installDir)
    args.append('--enable-shared')
    args.append('--enable-debug')
    args.append('--disable-opt')
    args.append('--enable-implicit-template-instantiation')
    args.append('--disable-deprecated')
    args.append('--without-hypre')
    args.append('--without-blaslapack')
    args.append('--without-cubes')
    args.append('--without-eleven')
    args.append('--without-kinsol')
    args.append('--without-petsc')
    args.append('--without-sundials')
    args.append('--without-x')
    args.append('--with-doxygen')
    args.append('--with-dot')
    return args
