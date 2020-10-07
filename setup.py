from setuptools import setup, Command, Extension
from Cython.Build import cythonize
import re, sys, os, shutil, glob, subprocess

# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

hfk = Extension(
    name = 'zs_hfk/_hfk',
    sources = ['cython_src/_hfk.c'],
    include_dirs = ['ComputeHFKv2'],
    extra_link_args = ['-Llib', '-lhfk']
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
    file = 'cython_src/_hfk.pyx'
    library = 'lib/libhfk.a'
    if os.path.exists(file):
        cythonize([file])
    if not os.path.exists(library):
        subprocess.check_call(['make', '-C', 'ComputeHFKv2'])

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
    cmdclass = {
        'clean':HFKClean,
        },
    zip_safe = False,
    requires = ['spherogram']
)
                
