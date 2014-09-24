import config.package

class Configure(config.package.GNUPackage):
  def __init__(self, framework):
    config.package.GNUPackage.__init__(self, framework)
    self.download  = ['https://docs.google.com/uc?export=download&confirm=no_antivirus&id=0BzuB-ydOOoduZjlFOEFRREZrT2s']
    self.functions = []
    self.includes  = ['muParser.h']
    self.liblist   = [['libmuparser.a']]
    self.pkgname   = 'muparser-2.2.3'
    self.cxx       = 1
    return

  def formGNUConfigureArgs(self):
    args = config.package.GNUPackage.formGNUConfigureArgs(self)
    args.append('--enable-shared=no')
    args.append('--enable-samples=no')
    args.append('--enable-debug=no')
    return args
