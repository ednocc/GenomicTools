#!/usr/bin/env python

# From OLCTools at https://github.com/OLC-Bioinformatics/OLCTools

import os
from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="GenomicTools",
    version="0.7.0",
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

util_dir = Path(__file__).resolve().parent / "genomictools/utils/"
usr_bin_dir = Path.home() / "bin"

if usr_bin_dir.exists():
    os.symlink(util_dir / "change_origin.py", usr_bin_dir / "change_origin")
    print("Some scripts have been symlinked to your ~/bin directory")
