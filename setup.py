#!/usr/bin/env python

# From OLCTools at https://github.com/OLC-Bioinformatics/OLCTools

from setuptools import setup, find_packages

setup(
    name="GenomicTools",
    version="0.5.0",
    packages=find_packages(),
    include_package_data=True,
    author="Cyril Conde",
    author_email="cyril.conde@gmail.com",
    url="https://github.com/ednocc/GenomicTools",
    python_requires='>=3.7',
    install_requires=['biopython',
                      'pyvcf',
                      'pysam',
                      'seaborn',
                      'pandas',
                      'numpy',
                      'pytest',
                      'xlsxwriter',
                      'joblib',
                      'dna_features_viewer']
)
