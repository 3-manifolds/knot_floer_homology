from setuptools import setup, Command, Extension
import re, sys, os

# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

hfk = Extension(
    name = 'zs_hfk/_hfk',
    sources = ['cython_src/_hfk.c'],
    include_dirs = ['ComputeHFKv2'],
    extra_link_args = ['-Llib', '-lhfk']
)

from Cython.Build import cythonize
if 'clean' not in sys.argv:
    file = 'cython_src/_hfk.pyx'
    if os.path.exists(file):
        cythonize([file])

setup(
    name='zs_hfk',
    version=version,
    author='Zoltán Szabó and Nathan M. Dunfield',
    author_email='nathan@dunfield.info',
    url='https://github/NathanDunfield/ZS_HFK',
    packages=['zs_hfk'],
    package_dir={'zs_hfk':'python_src'},
    package_data={'zs_hfk':['zs_hfk_binary']},
    ext_modules = [hfk],
    zip_safe = False,
    requires = ['spherogram']
)
                
