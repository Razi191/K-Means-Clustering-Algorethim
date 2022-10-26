from distutils.core import setup, Extension

setup(name='myspkmeans',
      version='1.0',
      description='myspkmeans for sp class',
      ext_modules=[Extension('myspkmeans', sources=['spkmeansmodule.c', 'spkmeans.c'])])
