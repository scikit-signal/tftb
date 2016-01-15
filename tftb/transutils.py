#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import re
import textwrap
import os.path as op
import subprocess
from scipy.io import loadmat
import numpy as np

OUTPUTARGPATTERN = r"\[.+\]"
INPUTARGPATTERN = r"\(.+\)"
FUNCNAMEPATTERN = r"\=.+\("


def compare_mlabarray(x, scriptname, matfile, varname, path=None, rtol=1e-5,
                      atol=1e-8):
    """Compare a MATLAB generated array with a NumPy ndarray.

    :param x: NumPy ndarray
    :param scriptname: MATLAB code that generates an array.
    :param varname: MATLAB variable name
    :param path: path to directory containing script
    :type x: array-like
    :type scriptname: str
    :type varname: str
    :type path: str
    :return: whether the arrays are identical
    :rtype: bool
    """
    script = op.basename(scriptname).replace(".m", "")
    command = "matlab -nodesktop -nojvm -nosplash -nodisplay -r {}"
    cmd = command.format(script).split()
    if path is not None:
        dirname = op.abspath(op.dirname(path))
    else:
        dirname = op.abspath(op.dirname(scriptname))
    subprocess.check_output(cmd, cwd=dirname)
    matarray = loadmat(matfile).get(varname)
    return np.allclose(matarray, x, rtol=rtol, atol=atol)


def get_function_signature(mpath):
    """Translate a MATLAB function signature into a Python function signature.

    :param mpath: path to MATLAB file
    :type mpath: str
    :return: Python function signature
    :rtype: str
    """
    with open(mpath, 'r') as fid:
        mcode = fid.readlines()
    if mcode[0].startswith('function'):
        outpath = op.basename(mpath).replace('.m', '.py')
        pyfuncstring = translate_funcline(mcode[0])
        with open(outpath, 'w') as fid:
            fid.write(pyfuncstring)
        return pyfuncstring


def translate_funcline(funcline):
    pyfuncstring = textwrap.dedent("""\
            def {funcname}({input_arguments}):


                return {output_arguments}
            """)
    match = re.search(OUTPUTARGPATTERN, funcline).group()
    outargs = match.replace('[', '').replace(']', '').replace(',', ', ')
    match = re.search(FUNCNAMEPATTERN, funcline).group()
    funcname = match.replace('=', '').replace('(', '')
    match = re.search(INPUTARGPATTERN, funcline).group()
    inargs = match.replace('(', '').replace(')', '').replace(',', ', ')
    return pyfuncstring.format(funcname=funcname, output_arguments=outargs,
                               input_arguments=inargs)


if __name__ == '__main__':
    print(compare_mlabarray(np.eye(3), '/tmp/foo.m', '/tmp/foo.mat', 'x'))
