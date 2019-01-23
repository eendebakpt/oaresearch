#!/usr/bin/env python

"""
setup.py file for OApackage
"""

# %% Load packages
from setuptools import setup, find_packages

import os
import re

def readme():
    with open('README.md') as f:
        return f.read()


def get_version(verbose=1, filename='oaresearch/version.py'):
    """ Extract version information from source code """

    with open(filename, 'r') as f:
        ln = f.readline()
        m = re.search('.* ''(.*)''', ln)
        version = (m.group(1)).strip('\'')
    if verbose:
        print('get_version: %s' % version)
    return version


packages = find_packages()

long_description = readme()

version = get_version()
print('package: version %s' % version)

setup(name='oaresearch',
      version=version,
      author="Pieter Eendebak",
      description="Research part of the OApackage",
      long_description=long_description,
      long_description_content_type='text/markdown',
      author_email='pieter.eendebak@gmail.com',
      license="BSD",
      url='http://www.pietereendebak.nl/oapackage/index.html',
      keywords=["orthogonal arrays, design of experiments, conference designs, isomorphism testing"],
      packages=packages,
      tests_require=['numpy', 'nose>=1.3', 'coverage>=4.0', 'mock'],
      zip_safe=False,
      install_requires=['numpy>=1.13', 'scanf', 'OApackage', 'json_tricks', 'jsmin', 'htmlmin'],
      extras_require={
          'GUI':  ["qtpy", 'matplotlib'],
      },
      requires=['numpy', 'matplotlib'],
      classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'License :: OSI Approved :: BSD License'
                   ]
      )
