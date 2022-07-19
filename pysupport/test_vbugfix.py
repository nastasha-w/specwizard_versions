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
                                        name1=flo, name2=fln)
        print('Match {} for z0 spectra, file {}'.format(same, first))
        allsame_z0 &= same

        fno = z3_old
        fnn = z3_new.format(first)
        flo = 'snap_z3_old'
        fln = 'snap_z3_new'
        same = ts.testsame_shortspectra(fno, fnn, specnums=sls, 
                                        name1=flo, name2=fln)
        print('Match {} for z3 spectra, file {}'.format(same, first))
        allsame_z3 &= same
    print('All same for z0: {}'.format(allsame_z0))
    print('All same for z3: {}'.format(allsame_z3))
    return allsame_z0 and allsame_z3


def testrightdiffs_longspectra(filen1, filen2, specnums='all',
                               name1='file1', name2='file2'):
    '''
    test if the specwizard outputs file1 and file2 are the same in the
    right places and different in the right places 
    (for the selected spectra) 
    useful to check if debugging had unwanted/unexpected side effects

    note that in principle different arrays could be similar if a sightline
    is pretty much empty, for example.
    
    input:
    ------
    filen1:   name of the first hdf5 output file
    filen2:   name of the second hdf5 output file
    specnums: list-like of integers or 'all'; which spectrum numbers to compare
              number are assumed to match between the files
    name1:    name for the first file (used in reporting differences)
    name2:    name for the second file (used in reporting differences)
    
    returns:
    --------
    True/False: whether all outputs match as they should(n't)
    
    more detailed differences are printed
    '''
    
    # ...same : should be the same, ...diff: should be different
    checkpaths_attrs_same = ['Constants',
                             'Header',
                             'Header/ModifyMetallicityParameters',
                             'Parameters/ChemicalElements',
                             'Parameters/SpecWizardRuntimeParameters',
                             'Units']
    checkpaths_attrs_diff = []
    checksame_arns = ['Wavelength_Ang']
    checkdiff_arns = []
    
    sgrps_rsmass_same = ['MetalMassFraction',
                         'OverDensity',
                         'Temperature_K']   
    sgrps_rsmass_diff = ['LOSPeculiarVelocity_KMpS'] 
    sgrps_ion_same = ['LogTotalIonColumnDensity']
    sgrps_ion_diff = ['RedshiftSpaceOpticalDepthOfStrongestTransition']
    sgrps_ionw_same = ['OverDensity',
                      'Temperature_K']
    sgrps_ionw_diff = ['LOSPeculiarVelocity_KMpS']
    gen_rsnionw_same = ['NIon_CM3']
    gen_rsnionw_diff = []

    sgrps_shortspecinfo_same =['FileUsed',
                               'Ibfactor',
                               'Icshift',
                               'LosUsed',
                               'NumberOfShortSpectra',
                               'RandomSeeds',
                               'x-axis',
                               'x_simunits',
                               'y-axis',
                               'y_simunits',
                               'z-axis',
                               ]
    sgrps_shortspecinfo_diff = []
    ionw_st_same = 'RealSpaceNionWeighted'
    tauw_st_diff = 'RedshiftSpaceOpticalDepthWeighted'
    
    asame_right = True
    didentical_right = True
    dsimilar_right = True
    adiff_right = True
    dnonidentical_right = True
    dnonsimilar_right = True

    msg_arrdiff_miss_id = '{n1} and {n2} unexpectedly had identical {path}' + \
                          'for spectrum {sn}'
    msg_arrdiff_miss_sm = '{n1} and {n2} unexpectedly had similar {path}' + \
                          'for spectrum {sn}'
    
    kwfmt = {'f1': name1, 'f2': name2}
    kw_allclose = {}
    kwa = {'name1': name1, 'name2': name2, 'kw_allclose': kw_allclose}
    
    with h5py.File(filen1, 'r') as f1, h5py.File(filen2, 'r') as f2:
        # check attributes
        for path in checkpaths_attrs_same:
            gn1 = name1 + ' ' + path
            gn2 = name2 + ' ' + path
            asame_right &= ts.checksame_attrs(f1[path], f2[path], 
                                              name1=gn1, name2=gn2)
        for path in checkpaths_attrs_diff:
            gn1 = name1 + ' ' + path
            gn2 = name2 + ' ' + path
            adiff_right &= ts.checkdiff_attrs(f1[path], f2[path], 
                                              name1=gn1, name2=gn2)
        # check arrays not specific to spectra
        for aname in checksame_arns:
            identical, similar = ts.checksame_arrays(f1, f2, aname, **kwa)
            didentical_right &= identical
            dsimilar_right &= similar
        for aname in checkdiff_arns:
            identical, similar = ts.checksame_arrays(f1, f2, aname, **kwa)
            dnonidentical_right &= not identical
            dnonsimilar_right &= not similar
            if identical:
                print(msg_arrdiff_miss_id.format(n1=name1, n2=name2, 
                                                 path=aname, 
                                                 specnum='general'))
            elif similar:
                print(msg_arrdiff_miss_sm.format(n1=name1, n2=name2, 
                                                 path=aname, 
                                                 specnum='general'))
        
        # find and check specnums
        if specnums == 'all':
            sks1 = set(f1.keys())
            sks2 = set(f2.keys())
            if not sks1 == sks2:
                print('{f1} and {f2} do not contain the same spectra'.format(**kwfmt))
                print('{f1} but not {f2}: {k1}'.format(k1=sorted(list(sks1 - sks2)), **kwfmt))
                print('{f2} but not {f1}: {k2}'.format(k2=sorted(list(sks2 - sks1)), **kwfmt))
                return False
            specgroups = set(sk if 'Spectrum' in sk else None for sk in sks1)
            specgroups -= {None}
            specgroups = sorted(list(specgroups), key=lambda grn: int(grn[8:]))
            specnums = np.array([int(grn[8:]) for grn in specgroups])
        else:
            specnums = np.array(sorted(list(specnums)))
            specgroups = ['Spectrum{sn}'.format(sn=sn) for sn in specnums]
            sks1 = set(f1.keys())
            sks2 = set(f2.keys())
            
            retf = False
            if not set(specgroups).issubset(sks1):
                print('{f1} missing requested spectra: {k1}'.format(\
                      k1=sorted(list(specgroups - sks1)), **kwfmt))
                retf = True
            if not set(specgroups).issubset(sks2):
                print('{f2} missing requested spectra: {k2}'.format(\
                      k1=sorted(list(specgroups - sks2)), **kwfmt))
                retf = True
            if retf:
                return False
        
        # check if output spectra are the same
        for specgroup in specgroups:
            g1 = f1[specgroup]
            g2 = f2[specgroup]
            asame_right &= ts.checksame_attrs(f1[path], f2[path],\
                                        name1=name1 + ' ' + specgroup,\
                                        name2=name2 + ' ' + specgroup)
            # loop over subgroups (if different, should be reported in attributes already)
            lkeys = list(set(g1.keys()) & set(g2.keys()))
            for ion in lkeys:
                if ion == 'RealSpaceMassWeighted':
                    for sgrp in sgrps_rsmass_same:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path,
                                                                 **kwa)
                        didentical_right &= identical
                        dsimilar_right &= similar
                    for sgrp in sgrps_rsmass_diff:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path,
                                                                 **kwa)
                        dnonidentical_right &= not identical
                        dnonsimilar_right &= not similar
                        if identical:
                            print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                        elif similar:
                            print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))
                elif ion == 'Redshift_RealSpace':
                    path = '/'.join([specgroup, ion])
                    identical, similar = ts.checksame_arrays(f1, f2, path,
                                                             **kwa)
                    didentical_right &= identical
                    dsimilar_right &= similar
                elif ion == 'Flux':
                    path = '/'.join([specgroup, ion])
                    identical, similar = ts.checksame_arrays(f1, f2, path, 
                                                             **kwa)
                    dnonidentical_right &= not identical
                    dnonsimilar_right &= not similar
                    if identical:
                        print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                    elif similar:
                        print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))
                elif ion == 'ShortSpectraInfo':
                    for sgrp in sgrps_shortspecinfo_same:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path, **kwa)
                        didentical_right &= identical
                        dsimilar_right &= similar
                    for sgrp in sgrps_shortspecinfo_diff:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path, **kwa)
                        dnonidentical_right &= not identical
                        dnonsimilar_right &= not similar
                        if identical:
                            print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                        elif similar:
                            print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))
                else:
                    for sgrp in sgrps_ion_same:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path, 
                                                                 **kwa)
                        didentical_right &= identical
                        dsimilar_right &= similar
                    for sgrp in sgrps_ion_diff:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = ts.checksame_arrays(f1, f2, path,
                                                                 **kwa)
                        dnonidentical_right &= not identical
                        dnonsimilar_right &= not similar
                        if identical:
                            print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                        elif similar:
                            print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))

                    if ionw_st_same in g1[ion] and ionw_st_same in g2[ion]:
                        for s2grp in sgrps_ionw_same:
                            path = '/'.join([specgroup, ion, ionw_st_same, 
                                             s2grp])
                            identical, similar = ts.checksame_arrays(f1, f2, 
                                                        path, **kwa)
                            didentical_right &= identical
                            dsimilar_right &= similar
                        for s2grp in sgrps_ionw_diff:
                            path = '/'.join([specgroup, ion, ionw_st_same, 
                                             s2grp])
                            identical, similar = ts.checksame_arrays(f1, f2, 
                                                        path, **kwa)
                            dnonidentical_right &= not identical
                            dnonsimilar_right &= not similar
                            if identical:
                                print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                            elif similar:
                                print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))
                        for s2grp in gen_rsnionw_same:
                            path = '/'.join([specgroup, ion, s2grp])
                            identical, similar = ts.checksame_arrays(f1, f2, 
                                                        path, **kwa)
                            didentical_right &= identical
                            dsimilar_right &= similar
                        for s2grp in gen_rsnionw_diff:
                            path = '/'.join([specgroup, ion, s2grp])
                            identical, similar = ts.checksame_arrays(f1, f2, 
                                                        path, **kwa)
                            dnonidentical_right &= not identical
                            dnonsimilar_right &= not similar
                            if identical:
                                print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                            elif similar:
                                print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))

                    if tauw_st_diff in g1[ion] and tauw_st_diff in g2[ion]:
                        for s2grp in sgrps_ionw_same + sgrps_ionw_diff:
                            path = '/'.join([specgroup, ion, tauw_st_diff, 
                                             s2grp])
                            identical, similar = ts.checksame_arrays(f1, f2, 
                                                        path, **kwa)
                            dnonidentical_right &= not identical
                            dnonsimilar_right &= not similar
                            if identical:
                                print(msg_arrdiff_miss_id.format(n1=name1, 
                                        n2=name2, path=path, 
                                        specnum=specgroup))
                            elif similar:
                                print(msg_arrdiff_miss_sm.format(n1=name1, 
                                        n2=name2,path=aname, 
                                        specnum=specgroup))
    if asame_right:
        print('All attributes that should match, did')
    if adiff_right:
        print('All attributes that should not match, did not')
    if didentical_right:
        print('All arrays that should be identical, were')
    if dnonsimilar_right:
        print('All arrays that should not be similar, were not')
    if dsimilar_right and not didentical_right:
        print('All arrays that should be were similar, though at least some were not identical')
    if dnonidentical_right and not dnonsimilar_right:
       print('All arrays that should be were not identical, though at least some were similar')
    return (asame_right and didentical_right and dsimilar_right and \
            adiff_right and dnonidentical_right and dnonsimilar_right) 


def test_match_longspectra():
    '''
    Check whether the right arrays and attributes are the same
    and different between the old and new versions
    '''
    z0_new = 'L0100N1504_longspectrum_z0p01_test_new.hdf5'
    z0_old = 'L0100N1504_longspectrum_z0p01_test_old.hdf5'
    z0_older = 'L0100N1504_longspectrum_z0p01_test_old_before-recompile.hdf5'
    z3_new = 'L0100N1504_longspectrum_z3p03_test_new.hdf5'
    z3_old = 'L0100N1504_longspectrum_z3p03_test_old.hdf5'
    z3_older = 'L0100N1504_longspectrum_z3p03_test_old_before-recompile.hdf5'
    
    print('\n')
    print('Old version, before/after recompile')
    ts.testsame_longspectra(ddir + z0_old, ddir + z0_older, specnums='all',
                            name1='z0_after_recompile', 
                            name2='z0_before_recompile')
    ts.testsame_longspectra(ddir + z3_old, ddir + z3_older, specnums='all',
                            name1='z3_after_recompile', 
                            name2='z3_before_recompile')
    print('\n')
    print('z close to 0: differences should be small')
    ts.testsame_longspectra(ddir + z0_old, ddir + z0_new, specnums='all',
                            name1='z0p01_old', 
                            name2='z0p01_new')
    print('\n\n')
    print('z=0.01 expected differences')
    testrightdiffs_longspectra(ddir + z0_old, ddir + z0_new, specnums='all',
                            name1='z0p01_old', 
                            name2='z0p01_new')
    print('\n')
    print('z=3 expected differences')
    testrightdiffs_longspectra(ddir + z0_old, ddir + z0_new, specnums='all',
                            name1='z0p01_old', 
                            name2='z0p01_new')

def checkdiff(ao, an, oldovernew, name):
    nani_o = np.isnan(ao)
    nani_n = np.isnan(an)
    iscorrect = True
    if not np.all(nani_o == nani_n):
        iscorrect = False
        print('NaN locations mismatched for {}'.format(name))
        print('Could not check a scaling for {}'.format(name))
    elif not np.allclose(ao[nani_o], an[nani_n] * oldovernew):
        iscorrect = False
        msg = 'Velocities did not have the expected difference for {}'
        print(msg.format(name))
    return iscorrect

def test_vdiff_longspecfiles(fno, fnn, aexp, losn, spectra='all'):
    '''
    Parameters:
    -----------
    fno: str
        name of the old version file (directory is ddir)
    fnn: str:
        name of the new version file (directory is ddir)
    aexp: float
        expansion factor for the files
    losn: str
        name of the los file used (excl. directories)
    spectra: list-like of ints or 'all'
        which spectra to compare.
    Returns:
    --------
    correct: bool
        True if the real space velocity fields differ by the right factor
        (np.allclose), False if not
    '''
    # 'the correct conversion factor for the los output is 1/a not sqrt(a)'
    # old = v_raw * sqrt(a), new = v_raw * 1./a
    # old / new = a**1.5
    oldovernew = aexp**1.5
    correct = True
    nskeys = {'Units', 'Parameters', 'Constants', 'Header', 'Wavelength_Ang'}
    with h5py.File(ddir + fno, 'r') as fo,\
         h5py.File(ddir + fnn, 'r') as fn:
        if spectra == 'all':
            _siter = set(fo.keys())
            if _siter != set(fn.keys()):
                msg = 'The sets of spectra do not match between {} and {}'
                print(msg.format(fno, fnn))
                _siter = _siter.intersection(set(fn.keys()))
            _siter -= nskeys
            _siter = list(_siter)
            _siter.sort(key=lambda x: int(x[8:])) #Spectrum<number>
        else:
            _siter = ['Spectrum{}'.format(x) for x in spectra]
        for skey in _siter:
            losn_used_o = fo['{}/ShortSpectraInfo/FileUsed'.format(skey)]
            losn_used_o = [n.decode() for n in losn_used_o]
            losn_used_n = fn['{}/ShortSpectraInfo/FileUsed'.format(skey)]
            losn_used_n = [n.decode() for n in losn_used_n]
            rightlosfile = np.all([losn == losn_u for losn_u in losn_used_o])
            # wrong LOS files means wrong aexp -> stop
            if not rightlosfile:
                msg = 'file {} used los file {} instead of expected {}'
                raise RuntimeError(msg.format(fno, losn_used_o, losn))
            rightlosfile = np.all([losn == losn_u for losn_u in losn_used_n])
            if not rightlosfile:
                msg = 'file {} used los file {} instead of expected {}'
                raise RuntimeError(msg.format(fnn, losn_used_n, losn))

            # mass-weighted velocity
            masswkey = '/RealSpaceMassWeighted/LOSPeculiarVelocity_KMpS'
            key = skey + masswkey
            mv_o = fo[key][:]
            mv_n = fn[key][:]
            correct &= checkdiff(mv_o, mv_n, oldovernew, key)
            
            # ion-weighted velocity
            nionkeys = {'Flux', 'RealSpaceMassWeighted', 
                        'Redshift_RealSpace', 'ShortSpectraInfo',
                        }
            ionkeys = set(fo[skey].keys())
            ionkeys_n = set(fn[skey].keys())
            ionkeys = ionkeys.intersection(ionkeys_n)
            ionkeys -= nionkeys
            ionkeys =  list(ionkeys)
            ionkeys.sort() # consistent order between spectra
            rsionkey = 'RealSpaceNionWeighted/LOSPeculiarVelocity_KMpS'
            for ionkey in ionkeys:
                key = '/'.join([skey, ionkey, rsionkey])
                iv_o = fo[key][:]
                iv_n = fn[key][:]
                correct &= checkdiff(iv_o, iv_n, oldovernew, key)
    return correct

def test_vdiff_longspectra():
    '''
    compare mass-weighted and ion-weighted real space velocities,
    check if they differ by the constant factor expected
    '''

    # values taken from the los files used 
    a_z0 = 0.9901333652195401
    a_z3 = 0.24832638738644483
    filen_los_z0 = 'part_los_z0.010.hdf5'
    filen_los_z3 = 'part_los_z3.027.hdf5'
    
    # only old/new, old before/after recompile checked -> same spectra
    fn_z0_new = 'L0100N1504_longspectrum_z0p01_test_new.hdf5'
    fn_z0_old = 'L0100N1504_longspectrum_z0p01_test_old.hdf5'
    fn_z3_new = 'L0100N1504_longspectrum_z3p03_test_new.hdf5'
    fn_z3_old = 'L0100N1504_longspectrum_z3p03_test_old.hdf5'

    z0res = test_vdiff_longspecfiles(fn_z0_old, fn_z0_new, a_z0, filen_los_z0, 
                                     spectra='all')
    z3res = test_vdiff_longspecfiles(fn_z3_old, fn_z3_new, a_z3, filen_los_z3, 
                                     spectra='all')
    if z0res:
        print('z=0.01 spectra passed the v ratio test')
    else:
        print('z=0.01 spectra failed the v ratio test')
    if z3res:
        print('z=3.03 spectra passed the v ratio test')
    else:
        print('z=3.03 spectra failed the v ratio test')

def plothists(ax, array, binsize_kmps=10.):
    array = np.array(array)
    minv = np.min(array[np.logical_not(np.isnan(array))])
    maxv = np.max(array[np.logical_not(np.isnan(array))])
    minbin = np.floor(minv / binsize_kmps) * binsize_kmps
    maxbin = np.ceil(maxv / binsize_kmps) * binsize_kmps
    print(minv, maxv)
    print(minbin, maxbin)
    bins = np.arange(minbin, maxbin + 0.5 * binsize_kmps, binsize_kmps)
    ax.hist(array.flatten(), bins=bins, align='mid', histtype='step',
            color='black', linewidth=2.5, log=True)
    for i, row in enumerate(array):
        color = 'C{}'.format(i % 10)
        if np.all(np.isnan(row)):
            continue
        ax.hist(row, bins=bins, align='mid', histtype='step', 
                linestyle='dashed', linewidth=1.5, color=color,
                log=True)
    
def makeimg_vdist(fn_l_o, fn_l_n, fn_s, imgname):
    lo_hists = {'mass': [], 'h1_p': [], 'h1_z': []}
    ln_hists = {'mass': [], 'h1_p': [], 'h1_z': []}
    s_hists = {'mass': [], 'h1_p': [], 'h1_z': []}
    nskeys_l = {'Units', 'Parameters', 'Constants', 
                'Header', 'Wavelength_Ang'}
    masswkey = '/RealSpaceMassWeighted/LOSPeculiarVelocity_KMpS'
    ionwpkey = '/RealSpaceNionWeighted/LOSPeculiarVelocity_KMpS'
    ionwzkey = '/RedshiftSpaceOpticalDepthWeighted/LOSPeculiarVelocity_KMpS'
    nskeys_s = {'Units', 'Header', 'Parameters', 'Constants',
                'VHubble_KMpS', 'Projection'}
    
    with h5py.File(fn_l_o) as flo:
        _siter = set(flo.keys())
        _siter -= nskeys_l
        _siter = list(_siter)
        _siter.sort(key=lambda x: int(x[8:])) #Spectrum<number>
        for skey in _siter:
            key = skey + masswkey
            mv = flo[key][:]     
            lo_hists['mass'].append(mv)

            key = skey + '/h1' + ionwpkey
            ivp = flo[key][:]  
            lo_hists['h1_p'].append(ivp)

            key = skey + '/h1' + ionwzkey
            ivz = flo[key][:]  
            lo_hists['h1_z'].append(ivz)

    with h5py.File(fn_l_n) as fln:
        _siter = set(fln.keys())
        _siter -= nskeys_l
        _siter = list(_siter)
        _siter.sort(key=lambda x: int(x[8:])) #Spectrum<number>
        for skey in _siter:
            key = skey + masswkey
            mv = fln[key][:]     
            ln_hists['mass'].append(mv)

            key = skey + '/h1' + ionwpkey
            ivp = fln[key][:]  
            ln_hists['h1_p'].append(ivp)

            key = skey + '/h1' + ionwzkey
            ivz = fln[key][:]  
            ln_hists['h1_z'].append(ivz)

    with h5py.File(fn_s) as fs:
        _siter = set(fs.keys())
        _siter -= nskeys_s
        _siter = list(_siter)
        _siter.sort(key=lambda x: int(x[8:])) #Spectrum<number>
        for skey in _siter:
            key = skey + masswkey
            mv = fs[key][:]     
            s_hists['mass'].append(mv)

            key = skey + '/h1' + ionwpkey
            ivp = fs[key][:]  
            s_hists['h1_p'].append(ivp)

            key = skey + '/h1' + ionwzkey
            ivz = fs[key][:]  
            s_hists['h1_z'].append(ivz)
    
    fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(11., 11.),
                             sharex=True, sharey=True)
    fig.suptitle('pixel values for mass and H I weighted velocity', 
                 fontsize=12)
    ax = axes[0, 0]
    ax.set_title('real space mass weighted\nshort spectra')
    plothists(ax, s_hists['mass'], binsize_kmps=10.)
    ax.set_ylabel('Pixel count', fontsize=12)
    ax = axes[1, 0]
    ax.set_title('new long spectra')
    plothists(ax, ln_hists['mass'], binsize_kmps=10.)
    ax.set_ylabel('Pixel count', fontsize=12)
    ax = axes[2, 0]
    ax.set_title('old long spectra')
    plothists(ax, lo_hists['mass'], binsize_kmps=10.)
    ax.set_ylabel('Pixel count', fontsize=12)
    ax.set_xlabel('Mass-weighted velocity [km/s]')

    ax = axes[0, 1]
    ax.set_title('real space H I weighted\nshort spectra')
    plothists(ax, s_hists['h1_p'], binsize_kmps=10.)
    ax = axes[1, 1]
    ax.set_title('new long spectra')
    plothists(ax, ln_hists['h1_p'], binsize_kmps=10.)
    ax = axes[2, 1]
    ax.set_title('old long spectra')
    plothists(ax, lo_hists['h1_p'], binsize_kmps=10.)
    ax.set_xlabel('H I-weighted velocity [km/s]')

    ax = axes[0, 2]
    ax.set_title('z space H I weighted\nshort spectra')
    plothists(ax, s_hists['h1_z'], binsize_kmps=10.)
    ax = axes[1, 2]
    ax.set_title('new long spectra')
    plothists(ax, ln_hists['h1_z'], binsize_kmps=10.)
    ax = axes[2, 2]
    ax.set_title('old long spectra')
    plothists(ax, lo_hists['h1_z'], binsize_kmps=10.)
    ax.set_xlabel('H I-weighted velocity [km/s]')

    plt.savefig(imgname, bbox_inches='tight')

def comparevdists_longshort():
    '''
    compare ion- and mass-weighted velocity distributions along
    the lines of sight of the long and short spectra;
    qualitative measure of velocity agreement
    '''
    # snapshots old/new were the same -> only need one
    z0_s_o = ddir + 'spec.snap_028_z000p000.0_old.hdf5'
    z3_s_o = ddir + 'spec.snap_012_z003p017.0_old.hdf5'
    #z0_s_n = ddir + 'spec.{:03d}.snap_028_z000p000.0.hdf5'
    #z3_s_n = ddir + 'spec.{:03d}.snap_012_z003p017.0.hdf5'

    z0_l_n = ddir + 'L0100N1504_longspectrum_z0p01_test_new.hdf5'
    z0_l_o = ddir + 'L0100N1504_longspectrum_z0p01_test_old.hdf5'
    z3_l_n = ddir + 'L0100N1504_longspectrum_z3p03_test_new.hdf5'
    z3_l_o = ddir + 'L0100N1504_longspectrum_z3p03_test_old.hdf5'
    
    makeimg_vdist(z0_l_o, z0_l_n, z0_s_o, mdir + 'vdist_comps_z0.pdf')
    makeimg_vdist(z3_l_o, z3_l_n, z3_s_o, mdir + 'vdist_comps_z3.pdf')


def comparelongspectrum(fno, fnn, speckey, imgname):
    '''
    compare a specific long spectrum. set up for z=3
    '''
    co = 'C1'
    cn = 'C0'
    
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(11., 8.),
                             gridspec_kw={'hspace': 0.5})
    fig.suptitle('Comparing old/new {}'.format(speckey), fontsize=12)
    
    fax = axes[0]
    tax = axes[1]
    iax = axes[2]
    vax = axes[3] 

    wllim = (4850., 4940.) #z3
    a_z3 = 0.24832638738644483

    fax.set_ylabel('Flux', fontsize=12)
    fax.set_xlabel('$\\lambda \\; [\\AA]$', fontsize=12)
    fax.set_xlim(*wllim)
    tax.set_ylabel('$\\tau$', fontsize=12)
    tax.set_xlabel('$\\lambda \\; [\\AA]$', fontsize=12)
    tax.set_yscale('log')
    tax.set_xlim(*wllim)
    iax.set_ylabel('$\\mathrm{n}_{\\mathrm{H\\,I}} \\; [\\mathrm{cm}^{-3}]$',
                   fontsize=12)
    iax.set_xlabel('position [(1 + Hubble flow redshift) * rest wavelength]',
                   fontsize=12)
    iax.set_yscale('log')
    iax.set_xlim(*wllim)
    vax.set_ylabel('$\\mathrm{v}_{\\mathrm{H\\,I}} \\;' + \
                   ' [\\mathrm{km} \\,/\\, \\mathrm{s}]$',
                   fontsize=12)
    vax.set_xlabel('position [(1 + Hubble flow redshift) * rest wavelength]',
                   fontsize=12)
    vax.set_xlim(*wllim)

    with h5py.File(fno, 'r') as f:
        wl_rest = f['Header'].attrs['Transitions_Rest_Wavelength'][0]
    
        color = co
        ls = 'solid'
        label = 'old'
        kwargs = {'color': color, 'linestyle': ls, 'label': label,
                  'linewidth': 2.}

        wl = f['Wavelength_Ang'][:]
        flux = f[speckey + '/Flux'][:]
        fax.plot(wl, flux, **kwargs)
        
        tzkey = '/h1/RedshiftSpaceOpticalDepthOfStrongestTransition'
        tau = f[speckey + tzkey][:]
        tax.plot(wl, tau, **kwargs)
        
        zp = f[speckey + '/Redshift_RealSpace'][:]
        wl_p = (1. + zp) * wl_rest

        nion = f[speckey + '/h1/NIon_CM3'][:]
        iax.plot(wl_p, nion, **kwargs)
        
        vpkey = '/h1/RealSpaceNionWeighted/LOSPeculiarVelocity_KMpS'
        vion = f[speckey + vpkey][:]
        vax.plot(wl_p, vion, **kwargs)
    
    with h5py.File(fnn, 'r') as f:
        _wl_rest = f['Header'].attrs['Transitions_Rest_Wavelength'][0]
        if not _wl_rest == wl_rest:
            raise RuntimeError('Spectra use different wavelengths:', 
                               wl_rest, _wl_rest)
        color = cn
        ls = 'dashed'
        label = 'new'
        kwargs = {'color': color, 'linestyle': ls, 'label': label,
                  'linewidth': 2.}

        wl = f['Wavelength_Ang'][:]
        flux = f[speckey + '/Flux'][:]
        fax.plot(wl, flux, **kwargs)
        
        tzkey = '/h1/RedshiftSpaceOpticalDepthOfStrongestTransition'
        tau = f[speckey + tzkey][:]
        tax.plot(wl, tau, **kwargs)
        
        zp = f[speckey + '/Redshift_RealSpace'][:]
        wl_p = (1. + zp) * wl_rest

        nion = f[speckey + '/h1/NIon_CM3'][:]
        iax.plot(wl_p, nion, **kwargs)
        
        vpkey = '/h1/RealSpaceNionWeighted/LOSPeculiarVelocity_KMpS'
        vion = f[speckey + vpkey][:]
        vax.plot(wl_p, vion, **kwargs)
    
    fax.legend(fontsize=12)
    vax2 = vax.twinx()
    vlim = np.array(vax.get_ylim())
    lightspeed_kmps = 299792458. * 1e-3
    dz = vlim / lightspeed_kmps / a_z3 * wl_rest
    vax2.set_ylim(*tuple(dz))
    vax2.set_ylabel('$\\Delta \\lambda \\; [\\AA]$')

    plt.savefig(imgname, bbox_inches='tight')

def comparelongspectra_sample():
    '''
    Compare flux and optical depth for a few sample long spectra, 
    alongside their line of sight ion densities and velocities
    ''' 
    z3_l_n = ddir + 'L0100N1504_longspectrum_z3p03_test_new.hdf5'
    z3_l_o = ddir + 'L0100N1504_longspectrum_z3p03_test_old.hdf5'
    
    spectra = ['Spectrum1', 'Spectrum2', 'Spectrum3']
    outname = mdir + 'z3_comp_{}.pdf'
    for spec in spectra:
        comparelongspectrum(z3_l_o, z3_l_n, spec, outname.format(spec))



def comparespectra_sarah()
    '''
    quick visual comparison between the spectra for Sarah Walsh's
    project before and after the fix.
    '''


