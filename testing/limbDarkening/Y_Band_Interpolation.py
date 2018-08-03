# Function to get y-band LDCs for any Teff, logg, M_H
# To be put into Aaron's code (EBLSST.py)
# written by Andrew Bowen, Northwestern undergraduate, funded by LSSTC grant (summer 2018)

def get_y_LDC(Teff, logg, M_H):
    
#     All filters/wavelength arrays
    SDSSfilters = ['u_','g_','r_','i_','z_', "J", 'H', "K" ]  #Only 2MASS/SDSS filters (8 in total)
    SDSSwavelength = np.array([354, 464, 621.5, 754.5, 870, 1220, 1630, 2190])
    y_wavelength = np.array(1004)
    
#     Getting coefficients from ELLC and appending them to specific coeff arrays
    SDSSfiltVals = np.array([])
    a1_array = np.array([])
    a2_array = np.array([])
    a3_array = np.array([])
    a4_array = np.array([])
#     Gets LDCs for all filters
    for w,f in zip(SDSSwavelength, SDSSfilters):
        ldy_filt = ellc.ldy.LimbGravityDarkeningCoeffs(f)
        a1, a2, a3, a4, y = ldy_filt(Teff, logg, M_H)
        a1_array = np.append(a1_array, a1)
        a2_array = np.append(a2_array, a2)
        a3_array = np.append(a3_array, a3)
        a4_array = np.append(a4_array, a4)

#     Sets up interpolation for y-band for each coeff
    find_y_a1 = np.interp(y_wavelength, SDSSwavelength, a1_array)
    find_y_a2 = np.interp(y_wavelength, SDSSwavelength, a2_array)
    find_y_a3 = np.interp(y_wavelength, SDSSwavelength, a3_array)
    find_y_a4 = np.interp(y_wavelength, SDSSwavelength, a4_array)
    
    return find_y_a1, find_y_a2, find_y_a3, find_y_a4
