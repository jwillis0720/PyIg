#!/usr/bin env python
from distutils.core import setup
import os
import pickle
from shutil import copytree
import sys

library_dir = "/usr/local/lib/"
pyig_lib = os.path.join(library_dir, 'pyig')

if not os.path.exists(pyig_lib):
    try:
      print "Trying...copying data to {0}".format(pyig_lib)
      copytree('data_files', pyig_lib)
      os.putenv('PYIGLIB',library_dir)
    except OSError:
        print "Need root to install in {0} directory location," \
              "use --library-path flag to specify where you want" \
              "to install the library".format(pyig_lib)
        sys.exit(1)


pickle.dump(pyig_lib,open('src/pyig/library_dir.txt','w'))
setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.gui', 'pyig.commandline', 'pyig'],
      package_dir={'': 'src'},
      package_data={'pyig.gui':['library_dir.txt']},
      scripts=['scripts/PyIg_GUI', 'scripts/PyIg'])
