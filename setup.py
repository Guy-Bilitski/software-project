from setuptools import setup, find_packages, Extension

setup(
    name="mykmeanssp",
    version = '1.0',
    author='Guy Bilitzki and Sagi Ahrac',
    author_email='sagiahrac@mail.tau.ac.il',
    description='C-API for spkmeans methods',
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
            name='mykmeanssp',
            sources=['spkmeansmodule.c'],
            depends=['spkmeans.h', 'spkmeans.c', 'kmeans.c', 'eigenvector.c', 'kmeans_io.c', 'matrix.c', 'point.c', 's_and_c.c', 'jacobi_output.c']
        )
    ]
)