from setuptools import setup, Command, Extension
from Cython.Build import cythonize
import re, sys, os, shutil, glob, subprocess

# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

cpp_dir = 'ComputeHFKv2'
cpp_sources = glob.glob(cpp_dir + '/*.cpp')
cpp_sources.remove(cpp_dir + '/Main.cpp')

hfk = Extension(
    name = 'zs_hfk/hfk',
    sources = ['cython_src/hfk.cpp'] + cpp_sources,
    include_dirs = [cpp_sources],
    extra_link_args = [],
    extra_compile_args = ['-O3', '-std=c++11'],
)

class HFKClean(Command):
    """
    Clean *all* the things!
    """
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        for dir in ['build', 'dist', 'lib']:
            shutil.rmtree(dir, ignore_errors=True)
        for file in glob.glob('*.pyc') + glob.glob('cython_src/*.c'):
            if os.path.exists(file):
                os.remove(file)

if 'clean' not in sys.argv:
    file = 'cython_src/hfk.pyx'
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
    ext_modules = [hfk],
    cmdclass = {
        'clean':HFKClean,
        },
    zip_safe = False,
    requires = ['spherogram']
)
                
