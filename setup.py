#!/usr/bin env python
from shutil import copytree, rmtree, copyfile
from distutils.core import setup
import platform
import os
import sys

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
            print "Can't find your operating system...Compile yourself, see documentation press any key to continue without installation of IgBlastN"
            raw_input()
    else:
        print "Can't find your architecture...I only have mac, \
              and linux supported in my files for now. Compile yourself. Press any key to continue without installation of Igblastn")
        raw_input()


#copy igblastn
os = get_os()
if os:
  src = "igblast/" + get_os() + "/bin/igblastn"
  src = os.path.abspath(src)
  try:
    copyfile(src, os.path.abspath(bin + "igblastn"))
    os.chmod(os.path.abspath(bin + "igblastn"), 0755)
  except IOError:
    raise IOError("You don't have correct permissions, please run as admin")

setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.commandline', 'pyig'],
      package_dir={'pyig': 'src/pyig'},
      scripts=['src/pyig/commandline/PyIg'])

