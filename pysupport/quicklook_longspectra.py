#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 17:42:13 2020

@author: Nastasha
"""

import numpy as np
import h5py

import matplotlib.pyplot as plt

def plotlongspectrum(fi, specnum, savename=None):
    '''
    plot the long spectrum in file fi, spectrum specnum 
    The figure object is not closed and can be viewed with plt.show() .
    
    input:
    ------
    fi:       open hdf5 file (specwizard output; h5py.File object)
    specnum:  spectrum number (int)
    savename: name of the file to save the plot to 
              (string or None -> not saved)
    
    returns:
    --------
    nothing
    '''
    plt.figure()
    if 'Wavelength_Ang' in fi:
        xv = fi['Wavelength_Ang'][:] 
        plt.xlabel('$\\lambda \\; [\\mathrm{\\AA}]$')
    else:
        xv = fi['Frequency_MHz'][:] 
        plt.xlabel('$\\nu \\; [\\mathrm{MHz}]$')
    yv = fi['Spectrum{sn}/Flux'.format(sn=specnum)][:]
    plt.ylabel('Flux')
    plt.plot(xv, yv)
    if savename is not None:
        plt.savefig(savename)
        
def plotlongspectra(filename, specnums='all', savename=None):
    '''
    plot long spectra in the specwizard output filename
    
    input:
    ------
    filename: name of the specwizard output file (string)
    specnums: 'all' or list-like of integers; which spectra in the file to 
              plot
    savename: name of the file(s) to save the spectra to
              None -> not saved
              string -> [/path/to/filename/]mainname.extension format is 
                        assumed. Spectra are saved as 
                        [/path/to/filename/]mainname_<specnum>.extension
              list of strings -> spectra are saved to these file names in 
                        order of spectrum number (sorted if 'all', 
                        otherwise matched by index to the list)
    returns:
    --------
    nothing                    
    '''
    fi = h5py.File(filename, 'r')
    if specnums == 'all':
        nspec = fi['Header'].attrs['NumberOfSpectra']
        specnums = range(1, nspec + 1) # specnums for long spectra start at 1
    elif not hasattr(specnums, '__len__'): # assume single number
        specnums = [specnums]
    else:
        specnums = list(specnums)
    if savename is None:
        savename = [savename] * len(specnums)
    elif not hasattr(savename[0], '__len__'): # one string
        parts = [name.split('.') for name in savename]
        savename = ['.'.join(parts[i][:-2] +\
                             [parts[i][-2] + '_{sn}'.format(sn=specnums[i])] +\
                             parts[-1]) \
                    for i in range(len(parts))]
    for i in range(len(specnums)):
        plotlongspectrum(fi, specnums[i], savename=savename[i])
        
def plotweightedspectrum(fi, specnum, savename=None,\
                         kind='RealSpaceMassWeighted', ion=None):
    '''
    input:
    ------
    fi:       open hdf5 file (specwizard output; h5py.File object)
    specnum:  spectrum number (int)
    savename: name of the file to save the plot to 
              (string or None -> not saved)
    kind:     what to plot --
              'RealSpaceMassWeighted', 'RealSpaceNionWeighted', or
              'RedshiftSpaceOpticalDepthWeighted'
    ion:      for ion- or optical-depth-weighted kinds: which ion to use
              (string)
    '''
    
    if 'RealSpace' in kind:
        xv = fi['Spectrum{sn}/Redshift_RealSpace'.format(sn=specnum)]
        xlabel = 'redshift (no peculiar velocities)'
    else:
        if 'Wavelength_Ang' in fi:
            xv = fi['Wavelength_Ang'][:] 
            xlabel = '$\\lambda \\; [\\mathrm{\\AA}]$'
        else:
            xv = fi['Frequency_MHz'][:] 
            xlabel = '$\\nu \\; [\\mathrm{MHz}]$'
    
    if kind == 'RealSpaceMassWeighted':
        title = 'Spectrum {sn}, real space mass-weighted'.format(sn=specnum)
        sgrp = fi['Spectrum{sn}/RealSpaceMassWeighted'.format(sn=specnum)]
        ykeys = ['rho', 'T', 'vpec', 'Z']
        wkey = 'rho'
        ylabels = {'rho': '$\\log_{10}(1 + \\delta)$',\
                   'T'  : '$\\log_{10}\\ \\mathrm{T} \\; [\\mathrm{K}]$',\
                   'vpec': '$\\mathrm{v}_{\\mathrm{p}} \\, [\\mathrm{km} \\, \\mathrm{s}^{-1}]$',\
                   'Z': '$\\log_{10} \\, \\mathrm{Z}$'}
        ytext = {'rho': 'weight',\
                 'T': None,\
                 'vpec': None,\
                 'Z':   'mass fraction'}
        data = {'rho': sgrp['OverDensity'][:],\
                'T':   sgrp['Temperature_K'][:],\
                'Z':   sgrp['MetalMassFraction'][:],\
                'vpec': sgrp['LOSPeculiarVelocity_KMpS'][:],\
                }
    elif kind == 'RealSpaceNionWeighted':
        title = 'Spectrum {sn}, real space {ion}-weighted'.format(sn=specnum,\
                                                                  ion=ion)
        sgrp = fi['Spectrum{sn}/{ion}'.format(sn=specnum, ion=ion)]
        ykeys = ['nion', 'rho', 'T', 'vpec']
        wkey = 'nion'
        ylabels = {'rho': '$\\log_{10}(1 + \\delta)$',\
                   'T'  : '$\\log_{10}\\ \\mathrm{T} \\; [\\mathrm{K}]$',\
                   'vpec': '$\\mathrm{v}_{\\mathrm{p}} \\, [\\mathrm{km} \\, \\mathrm{s}^{-1}]$',\
                   'nion': '$\\log_{{10}} \\, \\mathrm{{n}}(\\mathrm{{{ion}}}) \\; [\\mathrm{{cm}}^{{-3}}]$'.format(ion=ion)}
        ytext = {'nion': 'weight',\
                 'T': None,\
                 'vpec': None,\
                 'rho': None}
        data = {'rho': sgrp['RealSpaceNionWeighted/OverDensity'][:],\
                'T':   sgrp['RealSpaceNionWeighted/Temperature_K'][:],\
                'vpec': sgrp['RealSpaceNionWeighted/LOSPeculiarVelocity_KMpS'][:],\
                'nion': sgrp['NIon_CM3'][:],\
                }
    elif kind == 'RedshiftSpaceOpticalDepthWeighted':
        title = 'Spectrum {sn}, redshift space {ion}-weighted'.format(sn=specnum,\
                                                                  ion=ion)
        sgrp = fi['Spectrum{sn}/{ion}'.format(sn=specnum, ion=ion)]
        ykeys = ['tau_first', 'rho', 'T', 'vpec']
        wkey = 'tau_first'
        ylabels = {'rho': '$\\log_{10}(1 + \\delta)$',\
                   'T'  : '$\\log_{10}\\ \\mathrm{T} \\; [\\mathrm{K}]$',\
                   'vpec': '$\\mathrm{v}_{\\mathrm{p}} \\, [\\mathrm{km} \\, \\mathrm{s}^{-1}]$',\
                   'tau_first': '$\\log_{{10}} \\, \\tau(\\mathrm{{{ion}}})$'.format(ion=ion)}
        ytext = {'tau_first': 'weight',\
                 'T': None,\
                 'vpec': None,\
                 'rho': None}
        data = {'rho': sgrp['RedshiftSpaceOpticalDepthWeighted/OverDensity'][:],\
                'T':   sgrp['RedshiftSpaceOpticalDepthWeighted/Temperature_K'][:],\
                'vpec': sgrp['RedshiftSpaceOpticalDepthWeighted/LOSPeculiarVelocity_KMpS'][:],\
                'tau_first': sgrp['RedshiftSpaceOpticalDepthOfStrongestTransition'][:],\
                }
    
    fig, axes = plt.subplots(figsize=(5.5, 5.5), ncols=1, nrows=len(ykeys))
    fontsize = 12
    takelogs = ['rho', 'T', 'Z', 'nion', 'tau_first']
    minshows = {'T': 3.5,\
                'rho': -2.,\
                'Z':   -4.,\
                'tau_first': -3.,\
                'nion': -20.}
    
    fig.suptitle(title, fontsize=fontsize)
    for i in range(len(ykeys)):
        ax = axes[i]
        ykey = ykeys[i]
        doxlabel = i == len(ykeys) - 1
        if doxlabel:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabels[ykey], fontsize=fontsize)
        if ytext[ykey] is not None:
            ax.text(0.98, 0.98, ytext[ykey], fontsize=fontsize,\
                    transform=ax.transAxes,\
                    horizontalalignment='right', verticalalignment='top')
        yv = data[ykey]
        if ykey in takelogs:
            yv = np.log10(yv)
        ax.plot(xv, yv)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,\
                       labelbottom=doxlabel, labelleft=True,\
                       right=True, top=True)
        if ykey != wkey:
            # sanity check missing values
            assert np.all(np.isnan(data[ykey]) == (data[wkey] == 0.))
        if ykey in minshows:
            ylim = ax.get_ylim()
            if ylim[0] < minshows[ykey]:
                ymax = np.max(yv[np.isfinite(yv)])
                ymin = minshows[ykey]
                ax.set_ylim(ymin, ymax + 0.02 * (ymax - ymin))
        
    # sync x ranges
    xlims = [list(ax.get_xlim()) for ax in axes]
    xmin = min([xlim[0] for xlim in xlims])
    xmax = max([xlim[1] for xlim in xlims])
    [ax.set_xlim(xmin, xmax) for ax in axes]
    
    if savename is not None:
        plt.savefig(savename)