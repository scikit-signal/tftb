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

OUTPUTARGPATTERN = r"\[.+\]"
INPUTARGPATTERN = r"\(.+\)"
FUNCNAMEPATTERN = r"\=.+\("


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
    import sys
    print get_function_signature(sys.argv[1])
