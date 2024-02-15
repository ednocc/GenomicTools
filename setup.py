#!/usr/bin/env python

# From OLCTools at https://github.com/OLC-Bioinformatics/OLCTools

import os
from pathlib import Path
from setuptools import setup, find_packages

#'pyvcf'],
setup(
    name="GenomicTools",
    version="0.9.1",
    packages=find_packages(),
    include_package_data=True,
    author="Cyril Conde",
    author_email="cyril.conde@gmail.com",
    url="https://github.com/ednocc/GenomicTools",
    python_requires='>=3.7',
    install_requires=['biopython',
                      'pysam',
                      'seaborn',
                      'reportlab',
                      'pandas',
                      'numpy',
                      'pytest',
                      'xlsxwriter',
                      'joblib',
                      'dna_features_viewer'],
    entry_points={
        'console_scripts': [
            'change_origin = genomictools.utils.change_origin:main',
            'extract_region = genomictools.utils.extract_region:main'
            ]
        }
)

#util_dir = Path(__file__).resolve().parent / "genomictools/utils/"
#
## Edit for your configuration
#usr_bin_dir = Path.home() / "bin"
#
#if usr_bin_dir.exists():
#    try:
#        os.symlink(util_dir / "change_origin.py", usr_bin_dir / "change_origin")
#    except FileExistsError:
#        pass
#    print("Some scripts have been symlinked to your ~/bin directory")
