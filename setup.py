#!/usr/bin/env python

import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="seqcluster-helper",
    version="0.0.1",
    packages=["tailseq"],
    install_requires=["seqcluster",
                      "ipython-cluster-helper"],
    scripts=['scripts/seqcluster-helper.py']
)
