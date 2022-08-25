from setuptools import setup, find_packages, Extension

setup(
    name="mykmeanssp",
    version = '1.0',
    author='Guy Bilitzki and Sagi Ahrac',
    author_email='sagiahrac@mail.tau.ac.il',
    description='C-API for k-means func',
    packages=find_packages(),
    license='GPL-2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['kmeans.c']
        ),
    ]
)