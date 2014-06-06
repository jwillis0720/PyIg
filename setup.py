#!/usr/bin env python
from distutils.core import setup

setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.gui', 'pyig.commandline', 'pyig'],
      package_dir={'': 'src'},
      scripts=['scripts/PyIg_GUI', 'scripts/PyIg'])
