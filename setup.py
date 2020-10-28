from setuptools import setup, Command, Extension
from Cython.Build import cythonize
import re, sys, os, shutil, glob, subprocess

# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

cpp_dir = 'ComputeHFKv2'
cpp_sources = glob.glob(cpp_dir + os.sep + '*.cpp')
cpp_sources.remove(cpp_dir + os.sep + 'Main.cpp')

if sys.platform.startswith('win'):
    extra_compile_args = ['/Ox']
else:
    extra_compile_args = ['-O3', '-std=c++11']

hfk = Extension(
    name = 'knot_floer_homology.hfk',
    sources = ['cython_src/hfk.cpp'] + cpp_sources,
    include_dirs = [cpp_dir],
    extra_compile_args = extra_compile_args
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
        for file in glob.glob('*.pyc') + glob.glob('cython_src/*.cpp'):
            if os.path.exists(file):
                os.remove(file)

if 'clean' not in sys.argv:
    file = 'cython_src/hfk.pyx'
    if os.path.exists(file):
        cythonize([file])

setup(
    name='knot_floer_homology',
    version=version,
    author='Zoltán Szabó, Marc Culler, Nathan M. Dunfield, and Matthias Goerner',
    author_email='snappy-help@computop.org',
    url='https://github/3-manifolds/knot_floer_homology',
    packages=['knot_floer_homology'],
    package_dir={'knot_floer_homology':'python_src'},
    ext_modules = [hfk],
    cmdclass = {
        'clean':HFKClean,
        },
    zip_safe = False,
    requires = ['spherogram']
)
                
