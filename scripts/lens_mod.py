"""
Date: 24 Nov 2022
"""

from scripts.imports import *

def get_num_dens_gal(z, survey):
    """
    -- z = redshift value(s).
    -- survey = 'SKA2' or 'Euclid'
    
    output: dNdz, integrated dNdz, normalization.
    """
    if survey == 'SKA2':
        # SKA2 specification (from https://arxiv.org/abs/1601.03947)
        alpha, zm = np.sqrt(2), 1.3 # respectively scale parameter and the median redshift of sources
        z0 = zm/alpha 
        gamma = 1.25 
        beta = 2
        tot_ngal = 10 # total number of galaxies
        
    elif survey == 'Euclid':
        # Euclid specification (from https://arxiv.org/abs/1601.03947)
        alpha, zm = np.sqrt(2), 0.9 
        z0 = zm/alpha 
        gamma = 1.5 
        beta = 2
        tot_ngal = 30 
    
    # print('survey:',survey, ' alpha:',alpha,' zm:', zm, ' gamma:',gamma, ' beta:',beta, ' num gal:',tot_ngal)
    #redshift number density distribution 
    num_dens_gal = z**beta * np.exp(-(z/z0)**gamma) 
    num_gal = integrate.trapz(num_dens_gal, x = z) 
    K = tot_ngal/num_gal
    
    return num_dens_gal, num_gal, K

def get_dndz(z, alpha, zm, gamma, beta, tot_ngal):
    """
    -- z = redshift value(s).
    -- alpha, zm, gamma, beta, tot_ngal = see eq. (18) of https://arxiv.org/pdf/1601.03947.pdf
    
    output: dNdz, integrated dNdz, normalization.
    """
    z0 = zm/alpha 
    #redshift number density distribution 
    num_dens_gal = z**beta * np.exp(-(z/z0)**gamma) 
    num_gal = integrate.trapz(num_dens_gal, x = z) 
    K = tot_ngal/num_gal
    
    return num_dens_gal, num_gal, K
