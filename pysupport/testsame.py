#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 16:09:21 2020

@author: Nastasha
"""

import numpy as np
import h5py

def decodeif(x):
    if isinstance(x, type(b'')):
        out = x.decode()
    else:
        out = x
    return out

def checksame_attrs(group1, group2, name1='group1', name2='group2',):
    '''
    Check if the attributes of hdf5 two groups match exactly
    the exact differences are printed, but not returned
    
    input:
    ------
    group1: first hdf5 group
    group2: second hdf5 group
    name1:  name for the first hdf5 group (str; used in printing differences)
    name2:  name for the second hdf5 group (str; used in printing differences)
    
    returns:
    --------
    allsame (bool): True if all attributes present in one group are present in
                    the other, and their values are the same, otherwise False
    '''
    allsame = True
    kwfmt = {'g1': name1, 'g2': name2}
    
    attrs1 = {key: decodeif(val) for key, val in group1.attrs.items()}
    attrs2 = {key: decodeif(val) for key, val in group2.attrs.items()}
    
    keys1 = set(attrs1.keys())
    keys2 = set(attrs2.keys())
    keys_only1 = keys1 - keys2
    keys_only2 = keys2 - keys1
    keys_both = keys1 & keys2
    
    allsame = allsame and len(keys_only1) == 0 and len(keys_only2) == 0
    if len(keys_only1) > 0:
        print('{g1} contains keys {g2} does not:\n{keys}'.format(\
              keys=keys_only1, **kwfmt))
    if len(keys_only2) > 0:
        print('{g2} contains keys {g1} does not:\n{keys}'.format(\
              keys=keys_only2, **kwfmt))
    matchdct = {key: np.all(attrs1[key] == attrs2[key]) for key in keys_both}
    keymatch = np.all([matchdct[key] for key in matchdct])
    allsame &= keymatch
    if not keymatch:
        mismatch = {key if not matchdct[key] else None for key in matchdct}
        mismatch -= {None}
        print('{g1} and {g2} had different values for some attributes:'.format(\
              **kwfmt))
        for key in mismatch:
            print('{key}: \t {g1}: {at1} \t {g2}: {at2}'.format(\
                  key=key, at1=attrs1[key], at2=attrs2[key], **kwfmt))
    return allsame

def checkdiff_attrs(group1, group2, name1='group1', name2='group2',):
    '''
    Check if the attributes of hdf5 two groups all differ
    the exact differences are printed, but not returned
    
    input:
    ------
    group1: first hdf5 group
    group2: second hdf5 group
    name1:  name for the first hdf5 group (str; used in printing differences)
    name2:  name for the second hdf5 group (str; used in printing differences)
    
    returns:
    --------
    alldiff (bool): True if all attributes present in one group are present in
                    the other, and their values are the different, otherwise 
                    False
    '''
    alldiff = True
    kwfmt = {'g1': name1, 'g2': name2}
    
    attrs1 = {key: decodeif(val) for key, val in group1.attrs.items()}
    attrs2 = {key: decodeif(val) for key, val in group2.attrs.items()}
    
    keys1 = set(attrs1.keys())
    keys2 = set(attrs2.keys())
    keys_only1 = keys1 - keys2
    keys_only2 = keys2 - keys1
    keys_both = keys1 & keys2
    
    alldiff = alldiff and len(keys_only1) == 0 and len(keys_only2) == 0
    if len(keys_only1) > 0:
        print('{g1} contains keys {g2} does not:\n{keys}'.format(\
              keys=keys_only1, **kwfmt))
    if len(keys_only2) > 0:
        print('{g2} contains keys {g1} does not:\n{keys}'.format(\
              keys=keys_only2, **kwfmt))
    matchdct = {key: np.all(attrs1[key] != attrs2[key]) for key in keys_both}
    keymatch = np.all([matchdct[key] for key in matchdct])
    alldiff &= keymatch
    if not keymatch:
        mismatch = {key if not matchdct[key] else None for key in matchdct}
        mismatch -= {None}
        print('{g1} and {g2} had the same values for some attributes:'.format(\
              **kwfmt))
        for key in mismatch:
            print('{key}: \t {g1}: {at1} \t {g2}: {at2}'.format(\
                  key=key, at1=attrs1[key], at2=attrs2[key], **kwfmt))
    return alldiff

def checksame_arrays(f1, f2, path, name1='file1', name2='file2',\
                     kw_allclose=None):
    '''
    Check if the attributes of hdf5 two groups match exactly
    the exact differences are printed, but not returned
    
    input:
    ------
    f1:     first hdf5 File object
    f2:     second hdf5 File object
    path:   path to the array to compare
    name1:  name for the first hdf5 file (str; used in printing differences)
    name2:  name for the second hdf5 file (str; used in printing differences)
    kw_allclose: kwargs for numpy allclose, used to determine whether array
            values are similar
    
    returns:
    --------
    (identical, similar) (bool): whether the arrays are identical and whether 
                                 they are similar (numpy.allclose)
    '''
    similar = True
    identical = True
    kwfmt = {'f1': name1, 'f2': name2}
    
    a1 = np.array(f1[path])
    a2 = np.array(f2[path])
    try:
        n1 = np.isnan(a1)
        n2 = np.isnan(a2)
        if np.all(n1 == n2):
            nsel = np.logical_not(n1)
        else:
            print('{f1} and {f2} have differently positioned NaN values in {aname}'.format(\
                  aname=path, **kwfmt))
            return (False, False)
    except TypeError: # e.g. string array -> no NaNs to worry about
        nsel = ()
    identical = np.all(a1[nsel] == a2[nsel])
    if not identical:
        print('{f1} and {f2} do not have identical {aname}'.format(\
              aname=path, **kwfmt))
        similar = np.allclose(a1[nsel], a2[nsel], **kw_allclose)
        if not similar:
            print('{f1} and {f2} do not have similar {aname}'.format(\
              aname=path, **kwfmt))
            #print(a1)
            #print(a2)
    return (identical, similar)
    
def testsame_shortspectra(filen1, filen2, specnums='all',\
                          name1='file1', name2='file2'):
    '''
    test if the specwizard outputs file1 and file2 are the same (for the 
    selected spectra) 
    useful to check if debugging had unwanted/unexpected side effects
    
    Parameters:
    -----------
    filen1: str  
        name of the first hdf5 output file
    filen2: str   
        name of the second hdf5 output file
    specnums: list-like of integers or 'all'
        which spectrum numbers to compare. Numbers are assumed to match 
        between the files.
    name1: str   
        name for the first file (used in reporting differences)
    name2: str
        name for the second file (used in reporting differences)
    ignore_projection: bool
        ignore the 'Projection' group in the comparision. Useful in 
        quasar-style split-file outputs, where this data is only stored
        in the first file.
    
    returns:
    --------
    True/False: whether all outputs are identical 
    
    more detailed differences are printed
    '''
    
    checkpaths_attrs = ['Constants',
                        'Header',
                        'Header/ModifyMetallicityParameters',
                        'Parameters/ChemicalElements',
                        'Parameters/SpecWizardRuntimeParameters',
                        'Projection',
                        'Units']
    checksame_arns_sn = ['Projection/ncontr',
                         'Projection/x_fraction_array',
                         'Projection/y_fraction_array']
    checksame_arns = ['VHubble_KMpS']
    
    sgrps_mass = ['LOSPeculiarVelocity_KMpS',
                  'MetalMassFraction',
                  'OverDensity',
                  'Temperature_K']
    sgrps_ion = ['Flux',
                 'LogTotalIonColumnDensity',
                 'OpticalDepth']
    sgrps_ionw = ['LOSPeculiarVelocity_KMpS',
                  'NIon_CM3',
                  'OverDensity',
                  'Temperature_K']
    ionw_st = 'RealSpaceNionWeighted'
    tauw_st = 'RedshiftSpaceOpticalDepthWeighted'
    
    asame = True
    didentical = True
    dsimilar = True
    
    kwfmt = {'f1': name1, 'f2': name2}
    kw_allclose = {}
    kwa = {'name1': name1, 'name2': name2, 'kw_allclose': kw_allclose}
    
    with h5py.File(filen1, 'r') as f1, h5py.File(filen2, 'r') as f2:
        # check attributes
        for path in checkpaths_attrs:
            gn1 = name1 + ' ' + path
            gn2 = name2 + ' ' + path
            asame &= checksame_attrs(f1[path], f2[path], name1=gn1, name2=gn2)
        # check arrays not specific to spectra
        for aname in checksame_arns:
            identical, similar = checksame_arrays(f1, f2, aname, **kwa)
            didentical &= identical
            dsimilar &= similar
        
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
        
        # check general arrays pertaining to spectra
        for aname in checksame_arns_sn:
            a1 = f1[aname][:][specnums]
            a2 = f2[aname][:][specnums]
            
            identical = np.all(a1 == a2)
            didentical &= identical
            if not identical:
                print('{f1} and {f2} do not have identical {aname}'.format(\
                      aname=aname, **kwfmt))
                similar = np.allclose(a1, a2, **kw_allclose)
                dsimilar &= similar
                if not similar:
                    print('{f1} and {f2} do not have similar {aname}'.format(\
                      aname=aname, **kwfmt))
                    print(a1)
                    print(a2)
        
        # check if output spectra are the same
        for specgroup in specgroups:
            g1 = f1[specgroup]
            g2 = f2[specgroup]
            asame &= checksame_attrs(f1[path], f2[path],\
                                     name1=name1 + ' ' + specgroup,\
                                     name2=name2 + ' ' + specgroup)
            # loop over subgroups (if different, should be reported in attributes already)
            lkeys = list(set(g1.keys()) & set(g2.keys()))
            for ion in lkeys:
                if ion == 'RealSpaceMassWeighted':
                    sgrps = sgrps_mass
                    for sgrp in sgrps:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = checksame_arrays(f1, f2, path, **kwa)
                        didentical &= identical
                        dsimilar &= similar
                else:
                    sgrps = sgrps_ion
                    for sgrp in sgrps:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = checksame_arrays(f1, f2, path, **kwa)
                        didentical &= identical
                        dsimilar &= similar
                    if ionw_st in g1[ion] and ionw_st in g2[ion]:
                        for s2grp in sgrps_ionw:
                            path = '/'.join([specgroup, ion, ionw_st, s2grp])
                            identical, similar = checksame_arrays(f1, f2, path, **kwa)
                            didentical &= identical
                            dsimilar &= similar
                    if tauw_st in g1[ion] and tauw_st in g2[ion]:
                        for s2grp in sgrps_ionw:
                            path = '/'.join([specgroup, ion, tauw_st, s2grp])
                            identical, similar = checksame_arrays(f1, f2, path, **kwa)
                            didentical &= identical
                            dsimilar &= similar
    if asame:
        print('All attributes matched')
    if didentical:
        print('All arrays were identical')
    if dsimilar and not didentical:
        print('All arrays were similar, though at least some were not identical')
    return (asame and didentical and dsimilar)  
                
        
def testsame_longspectra(filen1, filen2, specnums='all',\
                          name1='file1', name2='file2'):
    '''
    test if the specwizard outputs file1 and file2 are the same (for the 
    selected spectra) 
    useful to check if debugging had unwanted/unexpected side effects
    
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
    True/False: whether all outputs are identical 
    
    more detailed differences are printed
    '''
    
    
    checkpaths_attrs = ['Constants',\
                        'Header',\
                        'Header/ModifyMetallicityParameters',\
                        'Parameters/ChemicalElements',\
                        'Parameters/SpecWizardRuntimeParameters',\
                        'Units']
    checksame_arns = ['Wavelength_Ang']
    
    sgrps_rsmass = ['LOSPeculiarVelocity_KMpS',\
                    'MetalMassFraction',\
                    'OverDensity',\
                    'Temperature_K']    
    sgrps_ion = ['LogTotalIonColumnDensity',\
                 'RedshiftSpaceOpticalDepthOfStrongestTransition']
    sgrps_ionw = ['LOSPeculiarVelocity_KMpS',\
                  'OverDensity',\
                  'Temperature_K']
    gen_rsnionw = ['NIon_CM3']
    
    sgrps_shortspecinfo =['FileUsed',\
                          'Ibfactor',\
                          'Icshift',\
                          'LosUsed',\
                          'NumberOfShortSpectra',\
                          'RandomSeeds',\
                          'x-axis',\
                          'x_simunits',\
                          'y-axis',\
                          'y_simunits',\
                          'z-axis',\
                          ]
    
    ionw_st = 'RealSpaceNionWeighted'
    tauw_st = 'RedshiftSpaceOpticalDepthWeighted'
    
    asame = True
    didentical = True
    dsimilar = True
    
    kwfmt = {'f1': name1, 'f2': name2}
    kw_allclose = {}
    kwa = {'name1': name1, 'name2': name2, 'kw_allclose': kw_allclose}
    
    with h5py.File(filen1, 'r') as f1, h5py.File(filen2, 'r') as f2:
        # check attributes
        for path in checkpaths_attrs:
            gn1 = name1 + ' ' + path
            gn2 = name2 + ' ' + path
            asame &= checksame_attrs(f1[path], f2[path], name1=gn1, name2=gn2)
        # check arrays not specific to spectra
        for aname in checksame_arns:
            identical, similar = checksame_arrays(f1, f2, aname, **kwa)
            didentical &= identical
            dsimilar &= similar
        
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
            asame &= checksame_attrs(f1[path], f2[path],\
                                     name1=name1 + ' ' + specgroup,\
                                     name2=name2 + ' ' + specgroup)
            # loop over subgroups (if different, should be reported in attributes already)
            lkeys = list(set(g1.keys()) & set(g2.keys()))
            for ion in lkeys:
                if ion == 'RealSpaceMassWeighted':
                    sgrps = sgrps_rsmass
                    for sgrp in sgrps:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = checksame_arrays(f1, f2, path, **kwa)
                        didentical &= identical
                        dsimilar &= similar
                elif ion == 'Redshift_RealSpace':
                    path = '/'.join([specgroup, ion])
                    identical, similar = checksame_arrays(f1, f2, path, **kwa)
                    didentical &= identical
                    dsimilar &= similar
                elif ion == 'Flux':
                    path = '/'.join([specgroup, ion])
                    identical, similar = checksame_arrays(f1, f2, path, **kwa)
                    didentical &= identical
                    dsimilar &= similar
                elif ion == 'ShortSpectraInfo':
                    for sgrp in sgrps_shortspecinfo:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = checksame_arrays(f1, f2, path, **kwa)
                        didentical &= identical
                        dsimilar &= similar
                else:
                    sgrps = sgrps_ion
                    for sgrp in sgrps:
                        path = '/'.join([specgroup, ion, sgrp])
                        identical, similar = checksame_arrays(f1, f2, path, **kwa)
                        didentical &= identical
                        dsimilar &= similar
                    if ionw_st in g1[ion] and ionw_st in g2[ion]:
                        for s2grp in sgrps_ionw:
                            path = '/'.join([specgroup, ion, ionw_st, s2grp])
                            identical, similar = checksame_arrays(f1, f2, path, **kwa)
                            didentical &= identical
                            dsimilar &= similar
                        for s2grp in gen_rsnionw:
                            path = '/'.join([specgroup, ion, s2grp])
                            identical, similar = checksame_arrays(f1, f2, path, **kwa)
                            didentical &= identical
                            dsimilar &= similar
                    if tauw_st in g1[ion] and tauw_st in g2[ion]:
                        for s2grp in sgrps_ionw:
                            path = '/'.join([specgroup, ion, tauw_st, s2grp])
                            identical, similar = checksame_arrays(f1, f2, path, **kwa)
                            didentical &= identical
                            dsimilar &= similar
    if asame:
        print('All attributes matched')
    if didentical:
        print('All arrays were identical')
    if dsimilar and not didentical:
        print('All arrays were similar, though at least some were not identical')
    return (asame and didentical and dsimilar) 

def recursivecomp_attrs(g1, g2, exclkeys):
    '''
    recursively check if all groups and datasets in an hdf5 file 
    have the same attributes. Any differences are printed.

    Parameters:
    -----------
    g1: str
        first hdf5 file
    g2: str
        second hdf5 file
    exclkeys: set of strings
        groups (at any single level) to exclude from the search

    Returns:
    --------
    None
    '''
    at1 = dict(g1.attrs.items())
    at2 = dict(g2.attrs.items())
    ak1 = set(at1.keys())
    ak2 = set(at2.keys())
    if ak1 == ak1:
        _ikeys = list(ak1)
    else:
        m1 = ak2 - ak1
        m2 = ak1 - ak1
        print('{} is missing attributes: {}'.format(g1, m1))
        print('{} is missing attributes: {}'.format(g2, m2))     
        _ikeys = list(ak1.intersection(ak2))
    for ak in _ikeys:
        if not np.all(at1[ak] == at2[ak]):
            print('{}, {} have different {}: {}, {}'.format(g1, g2, ak, at1[ak], at2[ak])) 
    if not hasattr(g1, 'keys'):
        if not hasattr(g2, 'keys'):
            return None
        else:
            print('{} has keys, {} does not'.format(g2, g1))
            return None
    elif not hasattr(g2, 'keys'):
        print('{} has keys, {} does not'.format(g1, g2))
        return None
    k1 = set(g1.keys())
    k2 = set(g2.keys())
    if k1 == k2:
        _ikeys = list(k1 - set(exclkeys))
    else:
        print('Keys in {}, {}, mismatch: {}, {}'.format(g1, g2, k1, k2))
        _ikeys = list(k1.intersection(k2) - set(exclkeys))
    if len(_ikeys) == 0:
        return None
    for key in _ikeys:
        recursivecomp_attrs(g1[key], g2[key], exclkeys)
    