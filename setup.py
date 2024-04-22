from os import path
from setuptools import setup, find_packages

VERSION = '1.0'
DESCRIPTION = 'Sparse Neighbor Joining: rapid phylogenetic inference using a sparse distance matrix'

directory = path.dirname(path.abspath(__file__))

with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

setup(
    name = 'SNJ',
    version = VERSION,
    description = DESCRIPTION,
    url = 'https://github.com/kurtsemih/SNJ',
    author = 'Semih Kurt',
    packages = find_packages(include = ['snj', 'snj.*']),
    entry_points = {
        'console_scripts': [
            'snj = snj.sparseNJ:main',
        ],
    },
    install_requires = required,
)