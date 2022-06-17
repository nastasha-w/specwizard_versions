#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 14:50:21 2021

@author: Nastasha
"""

import sys
import pandas as pd

basedir = '../'

def genparfiles_sarah(losnums, zqsos, setnum=1):
    if setnum == 1:
        templatename = 'L0100N1504_longspectrum_sarah_template1.par'
        outdir = basedir + 'parameterfiles/longspectra_sarah_set1/'
        outfile = outdir + 'L0100N1504_longspectrum_sarah_set1_los{num}.par'
        with open(basedir + 'parameterfiles/' + templatename, 'r') as fi:
            template = fi.read()
        for ln, zq in zip(losnums, zqsos):
            _outfile = outfile.format(num=ln)
            out = template.format(losnum=ln, zqso=zq)
            with open(_outfile, 'w') as fo:
                fo.write(out)
        print('done')
    else:
        raise ValueError('No prescription for setnum {}'.format(setnum))
        
def getparset(setnum=1):
    if setnum == 1:
        zfile = basedir + 'parameterfiles/longspectra_sarah_set1/' + \
                'random_z_set1.txt'
        zs = pd.read_csv(zfile, comment='#')
        print('generating files...')
        genparfiles_sarah(zs.index + 1, zs['redshift'], setnum=setnum)
    else:
        raise ValueError('No prescription for setnum {}'.format(setnum))

if __name__ == '__main__':
    setnum = int(sys.argv[1])
    getparset(setnum)
