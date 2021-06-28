import sys

__version__ = '0.1.3'

try:
    __TFTB__SETUP__
except NameError:
    __TFTB__SETUP__ = False

if __TFTB__SETUP__:
    sys.stderr.write('Partial import of tftb during the build process.\n')
else:
    from tftb import generators, processing, utils
    __all__ = ["generators", "processing", "utils"]
