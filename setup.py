# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 18:22:10 2023

@author: idajandl@gmx.de
"""

from setuptools import setup

setup(
    name='prim',
    version='0.1.0',    
    description='A python package to perform physical based retrieval for PRISMA datacubes',
    author='Ida Jandl',
    author_email='idajandl@gmx.de',
    license='MIT License',
    packages=['prim'],
    include_package_data=True,
    package_data={'':['gui/*']},
    install_requires=['matplotlib', 'numpy', 'netCDF4'],
    entry_points={
        'console_scripts': [
            'PRIM = prim.gui.run:run',
        ],
    },
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: MIT License',  
        'Operating System :: POSIX :: Linux/Windows',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)