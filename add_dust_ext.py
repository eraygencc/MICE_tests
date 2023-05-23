import numpy as np
from astropy.io import fits
import astropy.table
import scipy
import scipy.spatial
from scipy.spatial import kdtree
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import multiprocessing
from pathlib import Path
import time
import itertools

hdul=fits.open("MICE_cats/mice2_cat_1.fits")
data = hdul[1].data

x1 = data["ra_gal_mag"]
y1 = data["dec_gal_mag"]
z = data["z_cgal_v"]

hdul.close()

#Define the functions needed

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
   
zs = np.linspace(0.01,1.7,100000)
d_ang_func = interpolate.interp1d(zs, ((cosmo.angular_diameter_distance(zs) / u.Mpc) * 0.7) * 1000)

def sep_calculator(z_input):
    D_ang = ((cosmo.angular_diameter_distance(z_input) / u.Mpc) * 0.7) 
    return ((1 / D_ang) / 0.00029)

theta_max_func = interpolate.interp1d(zs, sep_calculator(zs))

def dist(x1,y1,x2,y2,z_lens,z_source): 
    return np.where(z_lens<z_source, ((np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )*0.00029)*(d_ang_func(z_lens)))**(-0.8), 0)

def do_kdtree(combined_x_y_arrays,points,theta):
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays) 
    results = mytree.query_ball_point(points,theta)
    return results

#Redshift boundaries of lenses
z_list = np.linspace(0.1,1.3,13)


for i in range(50):
    x_s = np.zeros(len(x1))
    x_s[np.where(z>z_list[i])] = x1[np.where(z>z_list[i])] * 60
    y_s = np.zeros(len(y1))
    y_s[np.where(z>z_list[i])] = y1[np.where(z>z_list[i])] * 60
    z_s = np.zeros(len(z))
    z_s[np.where(z>z_list[i])] = z[np.where(z>z_list[i])]
    x_l = x1[np.where((z>z_list[i])& (z<z_list[i+1]))] * 60
    y_l = y1[np.where((z>z_list[i])& (z<z_list[i+1]))] * 60
    z_l = z[np.where((z>z_list[i])& (z<z_list[i+1]))]

    comb_xy = np.dstack([x_s.ravel(),y_s.ravel()])[0]
    pts = np.dstack((x_l,y_l))
    points_list = list(pts[0])

    #calculate the apparent size of the maximum separation (1 Mpc h^-1) for the corresponding redshift
    theta_max = theta_max_func(np.mean(z_l)) 

    #find the nearest neighbours for the given maximum separation
    inds = do_kdtree(comb_xy,points_list,theta_max)
    
    add_seps = np.zeros(len(x_s))

    #calculate the scaled separations for each source and write into the corresponding index
    for l in range(len(x_l)):
        separations = dist_1(x_s[inds[l]], y_s[inds[l]], x_l[l], y_l[l], z_l[l], z_s[inds[l]])
        separations[separations==np.NaN] = 0
        separations[separations==np.inf] = 0
        separations[separations==-np.inf] = 0
        add_seps[inds[l]] = add_seps[inds[l]] + separations
        
    fn = "Txt_files/add_dust_"+ str(i) + ".txt"
    f = open(fn, 'w')

    for k in range(len(add_seps)):
         f.write(str(add_seps[k]) + "\n")
        
    f.close()
