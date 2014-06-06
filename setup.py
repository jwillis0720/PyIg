#!/usr/bin env python
from distutils.core import setup

import fnmatch
import os

matches = []
for root, dirnames, filenames in os.walk(''):
    for filename in fnmatch.filter(filenames, '*'):
        matches.append(os.path.join(root, filename))


setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.gui', 'pyig.commandline', 'pyig'],
      package_dir={'': 'src'},
      #package_data={'': ['database/*']},
      scripts=['scripts/PyIg_GUI', 'scripts/PyIg'])
