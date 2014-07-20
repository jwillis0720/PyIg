#!/usr/bin env python
from distutils.core import setup
import os
import pickle
from shutil import copytree, rmtree
import sys
import platform
# library_dir = "/usr/local/lib/"
# pyig_lib = os.path.join(library_dir, 'pyig')

# if os.path.exists(pyig_lib):
#   print "Deleting old copy of".format(pyig_lib)
#   rmtree(pyig_lib)
# try:
#   print "Trying...copying data to {0}".format(pyig_lib)
#   copytree('data_files', pyig_lib)
#   os.putenv('PYIGLIB',library_dir)
# except OSError:
#     print "Need root to install in {0} directory location," \
#           "use --library-path flag to specify where you want" \
#           "to install the library".format(pyig_lib)
#     sys.exit(1)

# #my hack way of knowking the directory location. This is supernice for dynamic code. THERE HAS TO BE A BETTER WAY OF DOING THIS
# pickle.dump(pyig_lib,open('src/pyig/gui/library_dir.txt','w'))


# def get_os():
#   _os, _op = sys.platform, platform.system()
#   if _os == 'darwin':
#     return 'darwin'
  
#   elif _os == 'Linux2':
#     if _op[0] == 'centos':
#       return 'centos'
#     elif _op[0] == 'redhat':
#       return 'redhat'
#     else:
#       raise Exception("Can't find your operating system...I only have redhat, mac, and centos in my files for now. Contact the higherups.")
#   else:
#     raise Exception("Can't find your architecture...I only have mac, and linux supported in my files for now. Contact the higherups.")


setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.commandline', 'pyig'],
      package_dir={'pyig': 'src/pyig'},
      package_data = {'pyig':['data_files/*.txt',
                              'data_files/*.py',
                              'data_files/database/Ig/human/*',
                              'data_files/database/Ig/mouse/*',
                              'data_files/database/TCR/human/*',
                              'data_files/database/TCR/mouse/*',
                              'data_files/internal_data/human/*',
                              'data_files/internal_data/mouse/*',
                              'data_files/internal_data/rabbit/*',
                              'data_files/internal_data/rat/*',
                              'data_files/optional_file/*',
                              ]},
      #package_data={'pyig.gui':['library_dir.txt']},
      scripts=['scripts/PyIg'])
      #'igblast_source/'+get_os()+"/bin/igblastn"])