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
gamma_1 = data["gamma1"]
gamma_2 = data["gamma2"]
z = data["z_cgal_v"]
z_B = data["Z_B"]
kappa = data["kappa"]
weight = data["weight"]

m_mag_u = data["m_evo_u"]
m_mag_g = data["m_evo_g"]
m_mag_r = data["m_evo_r"]
m_mag_i = data["m_evo_i"]

m_obs_u = data["m_obs_u"]
m_obs_g = data["m_obs_g"]
m_obs_r = data["m_obs_r"]
m_obs_i = data["m_obs_i"]

hdul.close()

m_evo_u = (m_mag_u + 2.5*np.log10(1+(2*kappa)))
m_evo_g = (m_mag_g + 2.5*np.log10(1+(2*kappa)))
m_evo_r = (m_mag_r + 2.5*np.log10(1+(2*kappa))) 
m_evo_i = (m_mag_i + 2.5*np.log10(1+(2*kappa)))

#ADDING EXTINCTION DUE TO INTERGALACTIC DUST

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    
#wavelengths [nm]
lambda_u = 365
lambda_g = 464
lambda_r = 658
lambda_i = 806


def dist(x1,y1,x2,y2,redshift): 
    theta=np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 ) * 0.00029 #from arcmin to radian
    D_ang = ((cosmo.angular_diameter_distance(redshift) / u.Mpc) * 0.7) * 1000  #from Mpc to kpc as Menard suggested in their paper
    r_p = theta * D_ang
    q=r_p**(-0.8)
    return q

def do_kdtree(combined_x_y_arrays,points,theta):
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays) 
    results = mytree.query_ball_point(points,theta)
    return results

#Middle redshifts of lenses and their corresponding maximum halo scales (1 Mpc h^-1)
z_list = np.linspace(0.1,1.3,13)
thetas_max_list = np.array([9.13504891,6.10805403,4.83488206,4.14683927,3.72506269,3.44682619,3.25479989,3.11859465,3.02059203,2.94986953, 2.89929382, 2.86400811, 2.84059113])


#for i in range(12):
def test_func(i):
    x_s = np.zeros(len(x1))
    x_s[np.where(z>z_list[i+1])] = x1[np.where(z>z_list[i+1])] * 60
    y_s = np.zeros(len(y1))
    y_s[np.where(z>z_list[i+1])] = y1[np.where(z>z_list[i+1])] * 60
    x_l = x1[np.where((z>z_list[i])& (z<z_list[i+1]))] * 60
    y_l = y1[np.where((z>z_list[i])& (z<z_list[i+1]))] * 60
    z_l = z[np.where((z>z_list[i])& (z<z_list[i+1]))]

    comb_xy = np.dstack([x_s.ravel(),y_s.ravel()])[0]
    pts = np.dstack((x_l,y_l))
    points_list = list(pts[0])

    theta_max = thetas_max_list[i]

    inds = do_kdtree(comb_xy,points_list,theta_max)
    
    add_seps = np.zeros(len(x_s))

    for l in range(len(x_l)):
        separations = dist(x_s[inds[l]], y_s[inds[l]], x_l[l], y_l[l], z_l[l])
        add_seps[inds[l]] = add_seps[inds[l]] + separations
            
        #add_mag_u = 1.09 * 2.3*(10**(-3)) * (551/lambda_u) * separation
        #add_mag_g = 1.09 * 2.3*(10**(-3)) * (551/lambda_g) * separation
        #add_mag_r = 1.09 * 2.3*(10**(-3)) * (551/lambda_r) * separation
        #add_mag_i = 1.09 * 2.3*(10**(-3)) * (551/lambda_i) * separation
        
        #add_ext_u[inds[l][j]]+=add_mag_u
        #add_ext_g[inds[l][j]]+=add_mag_g
        #add_ext_r[inds[l][j]]+=add_mag_r
        #add_ext_i[inds[l][j]]+=add_mag_i
        
    fn = "Txt_files/add_dust_"+ str(i) + ".txt"
    f = open(fn, 'w')

    for k in range(len(add_seps)):
         f.write(str(add_seps[k]) + "\n")
        
    f.close()
    return
               
#los_2 = range(10)
#arg=list(itertools.product(los_1,los_2))
arg = zip(range(len(z_list)-1))
pool = multiprocessing.Pool(processes=12)
pool.starmap(test_func, arg)
pool.close()
    
