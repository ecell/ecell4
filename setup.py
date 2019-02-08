import os.path
from setuptools import setup

DESCRIPTION = (
    "A software platform for modeling, simulation and analysis of complex, "
    "heterogeneous and multi-scale systems like the cell. E-Cell has "
    "multi-algorithm, multi-timescale and multi-spatial-representation as "
    "its central feature.")

LONG_DESCRIPTION = open(os.path.join(".", "README.md")).read()

setup(
    name='ecell4',
    version='1.0.dev3',
    packages=['ecell4', 'ecell4.util', 'ecell4.extra', 'ecell4.datasource'],
    url='https://github.com/ecell/ecell4',
    license='the GNU General Public License v2',
    author='Kazunari Kaizu',
    author_email='kaizu@riken.jp',
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    install_requires=['ecell4-base']
)
