from setuptools import setup, find_packages
from distutils.core import Extension
import os
import re

DISTNAME = 'WDRT'
PACKAGES = find_packages()
EXTENSIONS = []
DESCRIPTION = 'WEC Design Response Toolbox'
AUTHOR = 'WDRT developers'
MAINTAINER_EMAIL = ''
LICENSE = 'Revised BSD'
URL = 'http://wec-sim.github.io/WDRT/'
CLASSIFIERS=['Development Status :: 0 - Alpha',
             'Programming Language :: Python :: 3',
             'Topic :: Scientific/Engineering',
             'Intended Audience :: Science/Research',
             'Operating System :: OS Independent',
            ],
DEPENDENCIES = ['pandas', 
                'numpy', 
                'scipy',
                'matplotlib', 
                'requests', 
                'bs4', 
				'sklearn',
				'netCDF4', 
				'statsmodels', 
				'lxml', 
				'h5py']

# use README file as the long description
file_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(file_dir, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

# get version from __init__.py
with open(os.path.join(file_dir, 'WDRT', '__init__.py')) as f:
    version_file = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        VERSION = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")
        
setup(name=DISTNAME,
      version=VERSION,
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      maintainer_email=MAINTAINER_EMAIL,
      license=LICENSE,
      url=URL,
      classifiers=CLASSIFIERS,
      zip_safe=False,
      install_requires=DEPENDENCIES,
      scripts=[],
      include_package_data=True
  )
