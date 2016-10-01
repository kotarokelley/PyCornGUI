#!/usr/bin/python

try:
    from setuptools import setup
except ImportError:
    print("setuptools not found, falling back to distutils")
    from distutils.core import setup

setup(
    name='PyCornGUI',
    version='0.1dev',
    author='Kotaro Kelley',
    packages=['pycorngui'],
    
    extras_require = {'plotting':  ["matplotlib"], 'gui': ['wx']},
    platforms=['MacOSX'],
    zip_safe=False,
    classifiers=["License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 2.7",],
      
    package_data={'pycorngui': ['docs/*.*']},
    license='GNU General Public License v2 (GPLv2)',
    description='A GUI editor to extract and edit data from .res files generated \
        by UNICORN Chromatography software supplied with AKTA Systems.',
    url='https://github.com/kotarokelley/PyCornGUI',
    long_description = open('README.md').read()
)