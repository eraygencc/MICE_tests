from scipy.optimize import curve_fit
import pylab as py 
import numpy as np
import os
import glob
from astropy.io import fits
from astropy.wcs import WCS
import astropy.table
import scipy
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

import treecorr

from matplotlib import rc, rcParams
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)
rcParams['figure.figsize'] = (8., 6.)
rcParams['axes.linewidth'] = 2
rcParams['axes.labelsize'] = 24
rcParams['axes.titlesize'] = 30
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2
rcParams['lines.markersize'] = 5
rcParams['lines.markeredgewidth'] = 2
rcParams['xtick.major.size'] = 8
rcParams['xtick.minor.size'] = 4
rcParams['xtick.major.width'] = 2
rcParams['xtick.labelsize'] = 16
rcParams['ytick.major.size'] = 8
rcParams['ytick.minor.size'] = 4
rcParams['ytick.major.width'] = 2
rcParams['ytick.labelsize'] = 16

hdul=fits.open("MICE_cats/mice2_cat_1.fits")
data = hdul[1].data

x1 = data["ra_gal_mag"]
y1 = data["dec_gal_mag"]
gamma_1 = data["gamma1"]
gamma_2 = data["gamma2"]
z_true = data["z_cgal_v"]
z_B = data["Z_B"]
kappa = data["kappa"]
weight = data["weight"]

m_evo_u = data["m_evo_u"]
m_evo_g = data["m_evo_g"]
m_evo_r = data["m_evo_r"]
m_evo_i = data["m_evo_i"]

m_obs_u = data["m_obs_u"]
m_obs_g = data["m_obs_g"]
m_obs_r = data["m_obs_r"]
m_obs_i = data["m_obs_i"]

hdul.close()

###############################################
### needed for reddening calculations ###
lambda_u = 365
lambda_g = 464
lambda_r = 658
lambda_i = 806

total_sep = np.loadtxt('add_dust_final.txt')
###############################################
## Functions

thetas_ggl = np.array([ 1.197,1.7149 ,2.457,3.5202 ,5.0434 ,7.2258,10.353,14.832, 21.251,30.446, 43.621,62.497,89.54,128.29,183.8])
thetas_magn = np.array([0.16177,0.42335,1.1079 ,2.8994,7.5877,19.857,51.966,135.99])
mice_nonlin_x, mice_nonlin = np.loadtxt("Txt_files/mice_magn_nonlin_values.txt",usecols=(0,1),unpack=True)

def measure_ggl(ra_s,dec_s,ra_d,dec_d,redg_1s, redg_2s):
    cat_2 = treecorr.Catalog(x=ra_s*60,x_units='arcmin', y=dec_s*60, y_units='arcmin', g1=-redg_1s, g2=-redg_2s)
    cat_1 = treecorr.Catalog(x=ra_d*60,x_units='arcmin', y=dec_d*60, y_units='arcmin')

    fn_output = "Txt_files/GGL/ggl_"+str(i)+'_'+str(j)+'_'+str(k)+".txt"

    ng = treecorr.NGCorrelation(min_sep=1,max_sep=220,nbins=15,sep_units='arcmin')
    ng.process(cat_1,cat_2)
    ng.write(fn_output)

    ### PLOTTINGG ###
    (fig , ax) = plt.subplots()

    #ax.set_ylim(0.5e-4,5e-3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\theta$ in arcmin")
    plt.ylabel(r"$\gamma_\mathrm{t}$")

    plt.errorbar(thetas_ggl, ng.xi, yerr=ng.varxi, fmt="^--", color="black" ,label=r"$0.75<z<1.4$")
    
    fig_ggl = 'Figures/MICE_test/GGL/ggl_'+str(i)+'_'+str(j)+'_'+str(k)+'.png'

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(fig_ggl, bbox_inches='tight')
    return

def measure_magnification(ra_s,dec_s,ra_d,dec_d,magnitudes):
    mag_u = magnitudes[0]
    mag_g = magnitudes[1]
    mag_r = magnitudes[2]
    mag_i = magnitudes[3]
                
    mag_sample_u = mag_u[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_g = mag_g[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_r = mag_r[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_i = mag_i[np.where((z>z_min_s) & (z<z_max_s))]
                
    mag_input_u = mag_sample_u - np.mean(mag_sample_u)
    mag_input_g = mag_sample_g - np.mean(mag_sample_g)
    mag_input_r = mag_sample_r - np.mean(mag_sample_r)
    mag_input_i = mag_sample_i - np.mean(mag_sample_i)
    
    cat_2_u = treecorr.Catalog(x=ra_s,x_units='deg', y=dec_s, y_units='deg', k=mag_input_u)
    cat_2_g = treecorr.Catalog(x=ra_s,x_units='deg', y=dec_s, y_units='deg', k=mag_input_g)
    cat_2_r = treecorr.Catalog(x=ra_s,x_units='deg', y=dec_s, y_units='deg', k=mag_input_r)
    cat_2_i = treecorr.Catalog(x=ra_s,x_units='deg', y=dec_s, y_units='deg', k=mag_input_i)
    
    cat_1 = treecorr.Catalog(x=ra_d,x_units='deg', y=dec_d, y_units='deg')

    fn_output_u = "Txt_files/Magnification/magn_u_"+str(i)+'_'+str(j)+'_'+str(k)+".txt"
    nk_u = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_u.process(cat_1,cat_2_u)
    nk_u.write(fn_output_u)
    
    fn_output_g = "Txt_files/Magnification/magn_r_"+str(i)+'_'+str(j)+'_'+str(k)+".txt"
    nk_g = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_g.process(cat_1,cat_2_u)
    nk_g.write(fn_output_g)
    
    fn_output_r = "Txt_files/Magnification/magn_r_"+str(i)+'_'+str(j)+'_'+str(k)+".txt"
    nk_r = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_r.process(cat_1,cat_2_r)
    nk_r.write(fn_output_r)
    
    fn_output_i = "Txt_files/Magnification/magn_i_"+str(i)+'_'+str(j)+'_'+str(k)+".txt"
    nk_i = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_i.process(cat_1,cat_2_i)
    nk_i.write(fn_output_i)

    (fig, ax) = plt.subplots()

    plt.errorbar(thetas_magn,nk_u.xi,yerr=nk_u.varxi,fmt="^--",label=r'$u$') #[2:]
    plt.errorbar(thetas_magn,nk_g.xi,yerr=nk_g.varxi,fmt="^--",label=r'$g$')
    plt.errorbar(thetas_magn,nk_r.xi,yerr=nk_r.varxi,fmt="^--",label=r'$r$')
    plt.errorbar(thetas_magn,nk_i.xi,yerr=nk_i.varxi,fmt="^--",label=r'$i$')

    plt.plot(mice_nonlin_x, mice_nonlin,label="MICE")

    plt.xscale("log")
    #ax.set_ylim(-0.02,0.001)

    plt.xlabel(r"$\theta$ in arcmin")
    plt.ylabel(r"$\langle\delta_\mathrm{g}\Delta m \rangle$")

    plt.axhline(y = 0, color="black",linestyle = '--')
    plt.axvline(x = 1, color="black",linestyle = '--')
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig_magn = 'Figures/MICE_test/Magnification/magnification_'+str(i)+'_'+str(j)+'_'+str(k)+'.png'
    plt.savefig(fig_magn, bbox_inches='tight')
    return
    
def measure_reddening(ra_s,dec_s,ra_d,dec_d,magnitudes):
    mag_u = magnitudes[0]
    mag_g = magnitudes[1]
    mag_r = magnitudes[2]
    mag_i = magnitudes[3]
                
    mag_sample_u = mag_u[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_g = mag_g[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_r = mag_r[np.where((z>z_min_s) & (z<z_max_s))]
    mag_sample_i = mag_i[np.where((z>z_min_s) & (z<z_max_s))]
    
    color_ug = mag_sample_u - mag_sample_g
    color_gr = mag_sample_g - mag_sample_r
    color_ri = mag_sample_r - mag_sample_i
    color_gi = mag_sample_g - mag_sample_i
    
    color_input_ug = color_ug - np.mean(color_ug)
    color_input_gr = color_gr - np.mean(color_gr)
    color_input_ri = color_ri - np.mean(color_ri)
    color_input_gi = color_gi - np.mean(color_gi)

    cat_2_ug = treecorr.Catalog(x=ra_s,x_units='deg', y=ra_s, y_units='deg', k=color_input_ug)
    cat_2_gr = treecorr.Catalog(x=ra_s,x_units='deg', y=ra_s, y_units='deg', k=color_input_gr)
    cat_2_ri = treecorr.Catalog(x=ra_s,x_units='deg', y=ra_s, y_units='deg', k=color_input_ri)
    cat_2_gi = treecorr.Catalog(x=ra_s,x_units='deg', y=ra_s, y_units='deg', k=color_input_gi)
    
    cat_1 = treecorr.Catalog(x=ra_d,x_units='deg', y=ra_d, y_units='deg')

    fn_output_ug = "Txt_files/Reddening/redd_"+str(i)+'_'+str(j)+'_'+str(k)+"_ug.txt"
    nk_ug = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_ug.process(cat_1,cat_2_ug)
    nk_ug.write(fn_output_ug)
    
    fn_output_gr = "Txt_files/Reddening/redd_"+str(i)+'_'+str(j)+'_'+str(k)+"_gr.txt"
    nk_gr = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_gr.process(cat_1,cat_2_gr)
    nk_gr.write(fn_output_gr)
    
    fn_output_ri = "Txt_files/Reddening/redd_"+str(i)+'_'+str(j)+'_'+str(k)+"_ri.txt"
    nk_ri = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_ri.process(cat_1,cat_2_ri)
    nk_ri.write(fn_output_ri)
    
    fn_output_gi = "Txt_files/Reddening/redd_"+str(i)+'_'+str(j)+'_'+str(k)+"_gi.txt"
    nk_gi = treecorr.NKCorrelation(min_sep=0.1,max_sep=220,nbins=8,sep_units='arcmin')
    nk_gi.process(cat_1,cat_2_gi)
    nk_gi.write(fn_output_gi)

    (fig, ax) = plt.subplots()

    plt.errorbar(thetas_magn,nk_ug.xi,yerr=nk_ug.varxi,fmt="^--",label=r"$u-g$")
    plt.errorbar(thetas_magn,nk_gr.xi,yerr=nk_gr.varxi,fmt="^--",label=r"$g-r$")
    plt.errorbar(thetas_magn,nk_ri.xi,yerr=nk_ri.varxi,fmt="^--",label=r"$r-i$")
    plt.errorbar(thetas_magn,nk_gi.xi,yerr=nk_gi.varxi,fmt="^--",label=r"$g-i$")

    plt.xscale("log")
    #plt.yscale("log")
    #ax.set_ylim(-0.02,0.001)

    plt.xlabel(r"$\theta$ in arcmin")
    plt.ylabel(r"$\langle\delta_\mathrm{g}\delta E_{ij} \rangle$")

    plt.axhline(y = 0, color="black",linestyle = '--')
    plt.axvline(x = 1, color="black",linestyle = '--')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig_redd = 'Figures/MICE_test/Reddening/redd_'+str(i)+'_'+str(j)+'_'+str(k)+'.png'
    plt.savefig(figg_redd, bbox_inches='tight')
    return

###############################################

### Redshift: true/photo-z
z_list = [z_true, z_B]

### Lensing : on/off
#on:
redg_1_on = gamma_1/(1-kappa)
redg_2_on = gamma_2/(1-kappa)
#off:
gamma_1_random = np.array([np.random.uniform(-0.22,0.22) for i in range(len(x1))])
gamma_2_random = np.array([np.random.uniform(-0.22,0.22) for i in range(len(x1))])
redg_1_off = gamma_1_random/(1-kappa)
redg_2_off = gamma_2_random/(1-kappa)

redg1_list = [redg_1_on, redg_1_off]
redg2_list = [redg_2_on, redg_2_off]
redg_list = [redg1_list, redg2_list]

### MAGNITUDE SELECTION: Magnification, Extinction, Photometric Noise:
#no magnification, no extinction,no noise
mag_u_1 = m_evo_u + 2.5*np.log10(1+(2*kappa))
mag_g_1 = m_evo_g + 2.5*np.log10(1+(2*kappa))
mag_r_1 = m_evo_r + 2.5*np.log10(1+(2*kappa))
mag_i_1 = m_evo_i + 2.5*np.log10(1+(2*kappa))
mag_1_list = [mag_u_1, mag_g_1, mag_r_1, mag_i_1]

#magnified, no extinction, no noise,
mag_u_2 = m_evo_u 
mag_g_2 = m_evo_g
mag_r_2 = m_evo_r
mag_i_2 = m_evo_i
mag_2_list = [mag_u_2, mag_g_2, mag_r_2, mag_i_2]

#magnified, extinction included, no noise
mag_u_3 = m_evo_u + 1.09 * 2.3*(10**(-3)) * (551/lambda_u) * total_sep
mag_g_3 = m_evo_g + 1.09 * 2.3*(10**(-3)) * (551/lambda_g) * total_sep
mag_r_3 = m_evo_r + 1.09 * 2.3*(10**(-3)) * (551/lambda_r) * total_sep
mag_i_3 = m_evo_i + 1.09 * 2.3*(10**(-3)) * (551/lambda_i) * total_sep
mag_3_list = [mag_u_3, mag_g_3, mag_r_3, mag_i_3]

#magnified, extinction included, photometric noise included
mag_u_4 = m_obs_u + 1.09 * 2.3*(10**(-3)) * (551/lambda_u) * total_sep
mag_g_4 = m_obs_g + 1.09 * 2.3*(10**(-3)) * (551/lambda_g) * total_sep
mag_r_4 = m_obs_r + 1.09 * 2.3*(10**(-3)) * (551/lambda_r) * total_sep
mag_i_4 = m_obs_i + 1.09 * 2.3*(10**(-3)) * (551/lambda_i) * total_sep
mag_4_list = [mag_u_4, mag_g_4, mag_r_4, mag_i_4]

#no magnification, extinction included, no noise
mag_u_5 = mag_u_1 + 1.09 * 2.3*(10**(-3)) * (551/lambda_u) * total_sep
mag_g_5 = mag_g_1 + 1.09 * 2.3*(10**(-3)) * (551/lambda_g) * total_sep
mag_r_5 = mag_r_1 + 1.09 * 2.3*(10**(-3)) * (551/lambda_r) * total_sep
mag_i_5 = mag_i_1 + 1.09 * 2.3*(10**(-3)) * (551/lambda_i) * total_sep
mag_5_list = [mag_u_5, mag_g_5, mag_r_5, mag_i_5]

mag_list= [mag_1_list, mag_2_list, mag_3_list, mag_4_list, mag_5_list] 

### Shape noise: on/off


### Source selection: r-band magnitude limit: 21,22,23& source redshift selection

##### REDSHIFTS #####
z_min_s = 0.75
z_max_s = 1.4
z_min_l = 0.45
z_max_l = 0.55

############# MEASUREMENTS ###########

for i,z in enumerate(z_list):
    for j,redg in enumerate(redg_list):
        for k,mag in enumerate(mag_list):
            #SOURCES
            x1_s = x1[np.where((z>z_min_s) & (z<z_max_s))] #&(m_mag_r<25)
            y1_s = y1[np.where((z>z_min_s) & (z<z_max_s))]
            
            redg_1 = redg[0]
            redg_2 = redg[1]
            redg_1s = redg_1[np.where((z>z_min_s)& (z<z_max_s))]
            redg_2s = redg_2[np.where((z>z_min_s) & (z<z_max_s))]
            w_s = weight[np.where((z>z_min_s) & (z<z_max_s))]

            #LENSES
            x1_d = x1[np.where((z>z_min_l) & (z<z_max_l))] #& (M_halo < 15) & (M_halo>12)
            y1_d = y1[np.where((z>z_min_l) & (z<z_max_l))]

            ### GALAXY-GALAXY LENSING ###
            measure_ggl(x1_s,y1_s,x1_d,y1_d,redg_1s,redg_2s)

            ### MAGNIFICATION ###
            measure_magnification(x1_s,y1_s,x1_d,y1_d,mag)

            ### REDDENING ###
            measure_reddening(x1_s,y1_s,x1_d,y1_d,mag)