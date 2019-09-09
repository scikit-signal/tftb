import builtins
from setuptools import setup, find_packages


builtins.__TFTB__SETUP__ = True

# Setuptools config
NAME = "tftb"
DESCRIPTION = "Python module for time-frequency analysis."
with open('README.md') as f:
    LONG_DESCRIPTION = f.read()
MAINTAINER = 'Jaidev Deshpande'
MAINTAINER_EMAIL = 'deshpande.jaidev@gmail.com'
URL = "https://github.com/scikit-signal/tftb"
DOWNLOAD_URL = 'https://pypi.org/project/tftb/#files'
LICENSE = 'new BSD'
PROJECT_URLS = {
    'Bug Tracker': 'https://github.com/scikit-signal/tftb/issues',
    'Documentation': 'https://tftb.readthedocs.io',
    'Source Code': 'https://github.com/scikit-signal/tftb'
}

# Requirements
install_requires = [
    'numpy',
    'scipy',
    'matplotlib'
]

# Setup
import tftb  # NOQA: E402
setup(
    name=NAME,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    version=tftb.__version__,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=install_requires
)
