#!/usr/bin env python
from shutil import copytree, rmtree, copyfile
from distutils.core import setup
from os.path import expanduser
import os
import sys
import glob
import subprocess

#expand my home directory
home = expanduser("~")
library_dir = home + "/lib/pyig/data_dir"
bin_path = home + "/bin/"

print "Default Library Path is {0}".format(library_dir)
print "Default Bin Path is {0}".format(bin_path)

print "Checking Permissions for {0} and {1}".format(library_dir, bin_path)
if not os.access(library_dir, os.W_OK) and not os.access(bin_path, os.W_OK):
    raise OSError("Can't install libraries or binaries at {0} and {1}, permission error, seek administrator".format(
        library_dir, bin_path))


if sys.version_info < (2, 7,5):
    raise OSError("You need python 2.7.5 or greater")

try:
    import Bio
    print "Have {0}".format(Bio.__file__)
except Exception:
    raise ImportError("You need the Biopython Package...see documentation")


if os.path.exists(library_dir):
    print "Deleting old copy of {0}".format(library_dir)
    rmtree(library_dir)
    copytree('data_dir', library_dir)
else:
    copytree('data_dir', library_dir)


def get_igblast():
    print "Determining OS"
    igblasts = glob.glob('igblast/igblastn_*')
    for binary in igblasts:
        try:
            if subprocess.check_call([binary, '-h'],
                                     stdout=subprocess.PIPE) == 0:
                return os.path.abspath(binary)
        except OSError:
            continue
        return ""

# copy igblastn
igblast = get_igblast()
if igblast:
    new_igblast = bin_path + "/igblastn"
    print "Copying {0} to {1}".format(igblast, new_igblast)
    copyfile(igblast, new_igblast)
    print "Changing directory permissions of {0}".format(new_igblast)
    os.chmod(new_igblast, 0755)
else:
    print "Don't have a Igblastn that will run, \
           please see documentation to compile yourself, \
           press anykey to continue"
    raw_input()

#use python site_vax.setup.py --instal-lib ~/python_lib/ --install_script ~/bin/
setup(name='PyIg',
      version='1.1',
      description='Python Immunoglobulin Analysis Tools',
      author='Jordan Willis',
      author_email='jwillis0720@gmail.com',
      packages=['pyig.backend', 'pyig.commandline', 'pyig'],
      package_dir={'pyig': 'src/pyig'},
      scripts=['src/pyig/commandline/PyIg'])
