#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 19:43:55 2020

@author: Nastasha
"""

import numpy as np
import sys

losdir = '../los/'

def generate_random_los(number, filename, losdir=losdir):
    '''
    generate <number> random sightlines, save in specwizard los file format
    in file <filename> (directory: <losdir>)
    '''   
    fmtst = '{x:f}\t{y:f}\t0\n'
    vals = np.random.uniform(low=0.0, high=1.0, size=(number, 2))    
    with open(losdir + filename, 'w') as _f:
        _f.write('{num}\n'.format(num=number))
        for i in range(number):
            _f.write(fmtst.format(x=vals[i, 0], y=vals[i, 1])) 

if __name__ == '__main__':
    number = int(sys.argv[1])
    filename = sys.argv[2]
    if len(sys.argv) > 3:
        losdir = sys.argv[3]
    generate_random_los(number, filename, losdir=losdir)