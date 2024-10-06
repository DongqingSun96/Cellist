# -*- coding: utf-8 -*-
# @Author: dongqing
# @Date:   2023-10-14 12:54:35
# @Last Modified by:   dongqing
# @Last Modified time: 2023-10-14 16:07:36


# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2021-06-10 15:25:08
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-08-03 21:36:20


import sys,os

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.")

exec(open('src/Cellist/version.py').read())

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

def main():
    setup(
        name = "cellist",
        package_dir = {'':'src'},
        version = __version__,
        packages = find_packages(where="src"),
        scripts = ['bin/cellist'],
        package_data={
            "":["*.txt"]
        },
        install_requires = requirements,
        setup_requires = requirements,
        include_package_data = True,
        author = "Dongqing Sun",
        author_email = "Dongqingsun96@gmail.com",
        description = "Cellist (Cell identification in high-resolution Spatial Transcriptomics) is a cell segmentation tool for high-resolution spatial transcriptomics. ",
        license = "GPL-3.0",
        url = "https://github.com/dongqingsun96/Cellist",
        
        classifiers = [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
        python_requires=">=3.8",
    )

if __name__ == "__main__":
    main()
