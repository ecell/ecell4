from setuptools import setup

setup(
    name='ecell4-extras',
    version='0.1.0',
    packages=['ecell4_extras', 'ecell4_extras.util', 'ecell4_extras.extra', 'ecell4_extras.datasource'],
    url='https://github.com/ecell/ecell4-extras',
    license='the GNU General Public License v2',
    author='Kazunari Kaizu',
    author_email='kaizu@riken.jp',
    description='An extra package other than the core part of the cell simulation platform E-Cell4',
    install_requires=['ecell', 'SPARQLWrapper']
)
