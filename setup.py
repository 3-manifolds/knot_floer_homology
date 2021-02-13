from setuptools import setup, Command, Extension
import re, sys, os, shutil, glob, subprocess


# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

cpp_dir = 'ComputeHFKv2'
cpp_sources = glob.glob(cpp_dir + os.sep + '*.cpp')
unused_main_cpp = cpp_dir + os.sep + 'Main.cpp'
if unused_main_cpp in cpp_sources:
    cpp_sources.remove(unused_main_cpp)

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

try:
    from Cython.Build import cythonize
    if 'clean' not in sys.argv and cythonize is not None:
        file = 'cython_src/hfk.pyx'
        if os.path.exists(file):
            cythonize([file])
except ImportError:
    pass

setup(
    name='knot_floer_homology',
    version=version,
    author='Zoltán Szabó, Marc Culler, Nathan M. Dunfield, and Matthias Goerner',
    author_email='snappy-help@computop.org',
    url='https://github.com/3-manifolds/knot_floer_homology',
    description="Python wrapper for Zoltán Szabó's HFK Calculator",
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    packages=['knot_floer_homology'],
    package_dir={'knot_floer_homology':'python_src'},
    package_data={'knot_floer_homology':['HFK_data.json']},
    ext_modules = [hfk],
    cmdclass = {
        'clean':HFKClean,
        },
    zip_safe = False,
    requires = ['spherogram']
)
                
