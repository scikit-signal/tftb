from setuptools import setup

NAME = "tftb"

setup(
    name=NAME,
    version='0.0.1',
    author='Jaidev Deshpande',
    author_email='deshpande.jaidev@gmail.com',
    packages=['tftb',
              'tftb.processing',
              'tftb.generators',
              'tftb.tests']
)
