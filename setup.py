"""A setuptools based setup module.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from sphinx.setup_command import BuildDoc

cmdclass = {'build_sphinx': BuildDoc}

name = 'metagen_utils',
version = '1.21'

setup(
name='metagen_utils',
version='1.21',
description='A collection of tools and functions work with [ancient] metagenomics data',
url='https://github.com/mpieva/sediment_pipeline',
author='Frederic Romagne',
author_email='frederic.romagne@eva.mpg.de',
license='MIT',
classifiers=[

    'Development Status :: 3 - Beta',

    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',

    'License :: OSI Approved :: MIT License',

    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3 :: Only',
],
keywords='simulation metagenomic ancient genome',

packages=find_packages(include=['metagen_utils']),

install_requires=['biopython', 'numpy', 'pandas', 'pysam', 'sphinx'],
python_requires='>=3.6',
entry_points={
    'console_scripts': [
        'chunk_genome=chunk_genome:main',
    ],
},
package_data={  # Optional
        'metagen_utils': ['data/length_distribution.tsv', 'data/substitution_matrix.tsv'],
        '': ['docs/*.rst'],
    },
    cmdclass=cmdclass,
    include_package_data=True,
)
