from setuptools import setup
import re

# Get version number from module
version = re.search("__version__ = '(.*)'",
                    open('python_src/__init__.py').read()).group(1)

setup(
    name='zs_hfk',
    version=version,
    author='Zoltán Szabó and Nathan M. Dunfield',
    author_email='nathan@dunfield.info',
    url='https://github/NathanDunfield/ZS_HFK',
    packages=['zs_hfk'],
    package_dir={'zs_hfk':'python_src'},
    package_data={'zs_hfk':['zs_hfk_binary']},
    zip_safe = False,
    requires = ['spherogram']
)
                
