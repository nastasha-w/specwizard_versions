#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import h5py 
import matplotlib.pyplot as plt

import testsame as ts

# laptop NU, test 1
ddir = '/Users/nastasha/ciera/tests/specwizard_vbugfix/data_test1/'
mdir = '/Users/nastasha/ciera/tests/specwizard_vbugfix/images/'
testn = 'test1'

# should be the same in old and new versions
def test_snapshots():
    '''
    returns:
    --------
    bool: all old/new spectrum pairs are the same, or not
    '''
    z0_old = ddir + 'spec.snap_028_z000p000.0_old.hdf5'
    z3_old = ddir + 'spec.snap_012_z003p017.0_old.hdf5'
    z0_new = ddir + 'spec.{:03d}.snap_028_z000p000.0.hdf5'
    z3_new = ddir + 'spec.{:03d}.snap_012_z003p017.0.hdf5'
    
    numsl = 32
    numf = 16
    allsame_z0 = True
    allsame_z3 = True
    for first in range(numf):
        sls = range(first, numsl, numf)

        fno = z0_old
        fnn = z0_new.format(first)
        flo = 'snap_z0_old'
        fln = 'snap_z0_new'
        same = ts.testsame_shortspectra(fno, fnn, specnums=sls, 
                                        name1=flo, name2=fln,
                                        ignore_projection=False)
        print('Match {} for z0 spectra, file {}'.format(same, first))
        allsame_z0 &= same

        fno = z3_old
        fnn = z3_new.format(first)
        flo = 'snap_z3_old'
        fln = 'snap_z3_new'
        same = ts.testsame_shortspectra(fno, fnn, specnums=sls, 
                                        name1=flo, name2=fln,
                                        ignore_projection=False)
        print('Match {} for z3 spectra, file {}'.format(same, first))
        allsame_z3 &= same
    print('All same for z0: {}'.format(allsame_z0))
    print('All same for z3: {}'.format(allsame_z3))
    return allsame_z0 and allsame_z3



