from setuptools import setup

NAME = "tftb"

install_requires = [
    'scipy',
    'matplotlib',
    'scikit-image'
]

setup(
    name=NAME,
    version='0.1',
    author='Jaidev Deshpande',
    author_email='deshpande.jaidev@gmail.com',
    packages=['tftb',
              'tftb.processing',
              'tftb.generators',
              'tftb.tests'],
    install_requires=install_requires,
)
