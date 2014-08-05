#!/usr/bin env python
from distutils.core import setup
import os
from shutil import copytree, rmtree, copyfile
import sys
import platform

library_dir = "/usr/local/pyig/data_dir"
bin = "/usr/local/bin/"

if os.path.exists(library_dir):
      print "Deleting old copy of".format(library_dir)
      rmtree(library_dir)
try:
      print "Trying...copying data to {0}".format(library_dir)
      copytree('data_dir', library_dir)
except OSError:
      print "Need root to install in {0} directory location," \
             "to install the library".format(library_dir)
      sys.exit(1)


def get_os():
      _os, _op = sys.platform, platform.system()
      if _os == 'darwin':
            return 'darwin'
      elif _os == 'Linux2':
            if _op[0] == 'centos':
                  return 'centos'
            elif _op[0] == 'redhat':
                  return 'redhat'
            else:
                  raise Exception("Can't find your operating system...Compile yourself")
      else:
            raise Exception("Can't find your architecture...I only have mac, and linux supported in my files for now. Contact the higherups.")

#copy igblastn
src = "igblast/"+get_os()+"/bin/igblastn"
src = os.path.abspath(src)
try:
      copyfile(src,os.path.abspath(bin+"igblastn"))
      os.chmod(os.path.abspath(bin+"igblastn"), 0755)
except IOError:
      raise IOError("You don't have correct permissions, please run as admin")

setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.commandline', 'pyig'],
      package_dir={'pyig': 'src/pyig'},
      scripts=['scripts/PyIg'])
      




      #'igblast_source/'+get_os()+"/bin/igblastn"])



      # #package_data = {'pyig':['data_files/*.txt',
      #  #                       'data_files/*.py',
      #                         'data_files/database/Ig/human/*',
      #                         'data_files/database/Ig/mouse/*',
      #                         'data_files/database/TCR/human/*',
      #                         'data_files/database/TCR/mouse/*',
      #                         'data_files/internal_data/human/*',
      #                         'data_files/internal_data/mouse/*',
      #                         'data_files/internal_data/rabbit/*',
      #                         'data_files/internal_data/rat/*',
      #                         'data_files/optional_file/*',
                            #  ]},
      #package_data={'pyig.gui':['library_dir.txt']},