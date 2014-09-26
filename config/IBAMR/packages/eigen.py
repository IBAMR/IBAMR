import config.package

class Configure(config.package.CMakePackage):
  def __init__(self, framework):
    config.package.CMakePackage.__init__(self, framework)
    self.download  = ['http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2']
    self.functions = []
    self.includes  = ['eigen3/Eigen/Core']
    self.liblist   = [[]]
    self.pkgname   = 'eigen-3.2.2'
    self.cxx       = 1
    return
