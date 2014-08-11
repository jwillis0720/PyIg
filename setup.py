#!/usr/bin env python
from shutil import copytree, rmtree, copyfile
from distutils.core import setup
import platform
import os
import sys
import glob
import subprocess

library_dir = "/usr/local/pyig/data_dir"
bin = "/usr/local/bin/"


if sys.version_info < (2, 7):
    raise OSError("You need python 2.7.x")

try:
    import Bio
except ImportError:
    raise ImportError("You need the Biopython Package...see documentation")


if os.path.exists(library_dir):
    print "Deleting old copy of".format(library_dir)
    if os.access(library_dir, os.W_OK):
        rmtree(library_dir)
        copytree('data_dir', library_dir)
    else:
        print "Need to root to install {0}".format(library_dir)
else:
    try:
        copytree('data_dir', library_dir)
    except OSError:
        print "Need root to install in {0} directory location to install the library".format(library_dir)
        sys.exit(1)

def get_igblast():
    print "Determining OS"
    igblasts = glob.glob('igblast/igblastn_*')
    for binary in igblasts:
      try:
        if subprocess.check_call([binary,'-h'],stdout=subprocess.PIPE) == 0:
          return os.path.abspath(binary)
      except OSError:
          continue
      return ""

#copy igblastn
igblast = get_igblast()
if igblast:
  try:
    new_igblast = os.path.abspath(bin + "igblastn")
    print "Copying {0} to {1}".format(igblast,new_igblast)
    copyfile(igblast, new_igblast)
    print "Changing directory permissions of {0}".format(new_igblast)
    os.chmod(new_igblast, 0755)
  except IOError:
    raise IOError("You don't have correct permissions, please run as admin")

else:
  print "We don't have a IgBlast that will run, please see documentation to compile yourself, press anykey to continue"
  raw_input()

setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.commandline', 'pyig'],
      package_dir={'pyig': 'src/pyig'},
      scripts=['src/pyig/commandline/PyIg'])

