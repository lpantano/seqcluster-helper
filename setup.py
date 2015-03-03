#!/usr/bin/env python

import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="seqcluster-helper",
    version="0.99.05",
    packages=["sqhelper"],
    install_requires=["seqcluster",
                      "ipython-cluster-helper"],
    package_data={'sqhelper.templates': ['*.rmd']},
    include_package_data=True,
    scripts=['scripts/seqcluster-helper.py', 'scripts/seqcluster-installer.py']
)
