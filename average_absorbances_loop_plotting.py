# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 19:28:54 2025

@author: Laura SC
"""

import numpy as np
from normal_func import normal
from refract_normal_func import refract_normal
import matplotlib.pyplot as plt
import math

#convert degrees to radians
def degree2rad(deg):
    return((deg/360.0)*2.0*np.pi) #converts degrees to radians

#convert radians to degrees
def rad2degree(rad):
    return(rad*(180.0/np.pi)) #converts radians to degrees

#polar equations to calculate the coordinates of the ellipse, taking into account rotation around angle phi
def ellipse_polar(t, ax, ay, phi):
    x_prime = (ax * np.cos(t) * np.cos(phi)) - (ay * np.sin(t) * np.sin(phi)) #Rotated x coordinates
    y_prime = (ax * np.cos(t) * np.sin(phi)) + (ay * np.sin(t) * np.cos(phi)) #Rotated y coordinates
    return(x_prime, y_prime)

#first instance of when the incident light ray intersects the ellipse with magnitude and direction
def first_instance(inside_ri, outside_ri, incidence_flux_p, incidence_flux_s):
    # #Find normal
    ni, Ti, xni, yni = normal(ray, ellipse_params, normal_params, control_params)
    mag_ray = np.sqrt((ray[0][0]**2) + (ray[0][1]**2)) 
    unit_ray = ray[0]/mag_ray
    mag_ni = np.sqrt((ni[0]**2) + (ni[1]**2))
    unit_ni = ni/mag_ni
    ni_point = xni, yni #x and y coordinates of intersection of ray on ellipse

    #Regulate the direction of the normal according to direction of ray
    dot_uray_uni = unit_ni[0]*unit_ray[0] + unit_ni[1]*unit_ray[1]
    if dot_uray_uni > 0.0:
        unit_ni = -1.0*unit_ni
    else:
        unit_ni = unit_ni

    # #Initial refraction
    n1 = outside_ri
    n2 = inside_ri
    snell = n1/n2
    dot_uray_uni = unit_ni[0]*unit_ray[0] + unit_ni[1]*unit_ray[1]
    beta = unit_ray - dot_uray_uni*unit_ni
    gamma = (snell**2)*(1-(dot_uray_uni**2))
    unit_vt = snell*beta - unit_ni*np.abs(np.sqrt(1.0-gamma)) #unit vector for refracted ray transmitted into ellipse
    vt_origin = xni, yni #origin coordinates

    # #Transmitted reflection
    dot_ray_unit_ni = ray[0][0]*unit_ni[0] + ray[0][1]*unit_ni[1]
    vrx = ray[0][0] - (2*dot_ray_unit_ni*unit_ni[0])
    vry = ray[0][1] - (2*dot_ray_unit_ni*unit_ni[1])
    vr = vrx, vry #vector for the reflected ray inside ellipse
    mag_vr = np.sqrt((vr[0]*vr[0]) + (vr[1]*vr[1]))
    unit_vr = vr/mag_vr #unit vector of vr above

    # #Fresnell laws to calculate the flux intensities for ray encoding
    dot_init_incidence = unit_ray[0]*unit_ni[0] + unit_ray[1]*unit_ni[1]
    dot_init_transmittance = unit_ni[0]*unit_vt[0] + unit_ni[1]*unit_vt[1]
    n1 = outside_ri
    n2 = inside_ri
    init_rp = (abs(((1.0/n2 * dot_init_transmittance) - (1.0/n1 * dot_init_incidence)) \
                   / ((1.0/n2 * dot_init_transmittance) + (1.0/n1 * dot_init_incidence))))**2
    init_rs = (abs(((1.0/n2 * dot_init_incidence) - (1.0/n1 * dot_init_transmittance)) \
                   / ((1.0/n2 * dot_init_incidence) + (1.0/n1 * dot_init_transmittance))))**2
        #intensity flux of ray reflected away from ellipse
    init_flux_p_ref = init_rp * incidence_flux_p
    init_flux_s_ref = init_rs * incidence_flux_s
        #intensity flux of ray transmitted into ellipse
    init_flux_p_trans = (1.0-init_rp) * incidence_flux_p
    init_flux_s_trans = (1.0-init_rs) * incidence_flux_s 
    return(ni_point, unit_ni, unit_ray, unit_vt, vt_origin, vr, mag_vr, unit_vr, \
           init_flux_p_ref, init_flux_s_ref, init_flux_p_trans, init_flux_s_trans)

#following instances of intersection between the ray and the ellipse internal surface when the ray is
    #now originating from the inside of the ellipse
def subs_instances(inside_ri, outside_ri, unit_vt, vt_origin, prev_flux_p_trans, prev_flux_s_trans):
    # #INTERNAL REFLECTION
    #Find normal
    vt = unit_vt*mag_vt #vector for initial internal ray to be reflected
    vt_origin = xni, yni #intersection point on ellipse to be used as origin of ray
    vrf = (vt, vt_origin) #ray vectors format to be used by funtion below
    nvti, Tvti, xvti, yvti = refract_normal(vrf, ellipse_params, normal_params, control_params) 
    mag_nvti = np.sqrt((nvti[0]**2) + (nvti[1]**2))
    unit_nvti = nvti/mag_nvti #unit vector for the normal at the following intersection after origin of internal ray
    vti_points = xvti, yvti #intersection points for the second intersection of internal ray of interest

    #Regulate normal direction according to direction of ray
    dot_uvt_unvti = unit_nvti[0]*unit_vt[0] + unit_nvti[1]*unit_vt[1]
    if dot_uvt_unvti > 0.0:
        unit_nvti = -1.0*unit_nvti
    else:
        unit_nvti = unit_nvti
        
    # #Internal reflection
    dot = unit_vt[0]*unit_nvti[0] + unit_vt[1]*unit_nvti[1]
    nrix = unit_vt[0] - (2*dot*unit_nvti[0])
    nriy = unit_vt[1] - (2*dot*unit_nvti[1])
    nri = nrix, nriy
    mag_nri = np.sqrt((nri[0]*nri[0]) + (nri[1]*nri[1]))
    unit_nri = nri/mag_nri #unit vector for reflected ray

    # #Transmitted refraction out of ellipse
    n1 = inside_ri
    n2 = outside_ri
    snell = n1/n2
    dot_unri_unvti = unit_nvti[0]*unit_nri[0] + unit_nvti[1]*unit_nri[1]
    beta = unit_nri - dot_unri_unvti*unit_nvti
    gamma = (snell**2)*(1.0-(dot_unri_unvti**2))
    if gamma > 1.0: #checks for total internal reflection
        unit_vt2 = (0.0, 0.0) #no vector for the transmitted out ray since there is no ray transmitted out (all flux kept inside ellipse)
        vt_origin2 = xvti, yvti
    else:
        unit_vt2 = snell*beta - unit_nvti*np.abs(np.sqrt(1.0-gamma)) #unit vector for ray transmitted out of ellipse 
        vt_origin2 = xvti, yvti #intersection vector for location of where transmitted out ray originates

    # #Fresnell
    target = (0.0, 0.0)
    if tuple(unit_vt2) == target: #error handling for total internal reflection -> all flux is kept inside
        sequ_flux_p_ref = prev_flux_p_trans #all flux for p polarization is reflected back into ellipse
        sequ_flux_s_ref = prev_flux_s_trans #all flux for s polarization is reflected back into ellipse
        sequ_flux_p_trans = 0.0 #no p polarization is transmitted out
        sequ_flux_s_trans = 0.0 #no s polarization is transmitted out
    else:
        dot_incidence = unit_nvti[0]*unit_vt[0] + unit_nvti[1]*unit_vt[1]
        dot_transmittance = unit_nvti[0]*unit_vt2[0] + unit_nvti[1]*unit_vt2[1]
        n1 = inside_ri
        n2 = outside_ri
        rp = (abs(((1.0/n2 * dot_transmittance) - (1.0/n1 * dot_incidence)) \
                  / ((1.0/n2 * dot_transmittance) + (1.0/n1 * dot_incidence))))**2
        rs = (abs(((1.0/n2 * dot_incidence) - (1.0/n1 * dot_transmittance)) \
                  / ((1.0/n2 * dot_incidence) + (1.0/n1 * dot_transmittance))))**2
    
        sequ_flux_p_ref = rp * prev_flux_p_trans #p polarization that is reflected and stays within the ellipse
        sequ_flux_s_ref = rs * prev_flux_s_trans #s polarization that is reflected and stays within the ellipse
        
        sequ_flux_p_trans = (1.0-rp) * prev_flux_p_trans #p polarization that is reftracted and transmitted out of the ellipse
        sequ_flux_s_trans = (1.0-rs) * prev_flux_s_trans #s polarization that is reftracted and transmitted out of the ellipse
    return(vti_points, unit_nvti, nri, mag_nri, unit_nri, unit_vt2, vt_origin2, \
           sequ_flux_p_ref, sequ_flux_s_ref, sequ_flux_p_trans, sequ_flux_s_trans)

#################################################################################################################################################################################################################################################################
# PARAMETERS
#################################################################################################################################################################################################################################################################

#TANGENT & NORMAL PARAMETERS:
delta_t = 0.001 #Small change in t to calculate the tangent in normal calculation
tol = 0.0001 #Angle tolerance
ang_res = 10000.0 #Angle resolution
normal_params = (delta_t, tol, int(ang_res))

#CONTROL PARAMETERS:
t_range_param = 10000.0
plot_range = 10000.0 #Avoids errors with intersect calculations
control_params = (int(t_range_param), plot_range)

#REFRACTION PARAMETERS:
outside_ri = 1.0 #Refraction index of media outside ellipse
inside_ri = 1.4 ##Refraction index of media inside ellipse
mag_vt = 100000.0 #Gives magnitude for calculating vector of initial ray inside ellipse -> amplifies it

#RAY CODING:
incidence_flux_p = 500.0 #p polarization flux of incoming ray
incidence_flux_s = 500.0 #s polarization flux of incoming ray
mu = 1.0 #Absorbance parameter by the ellipse: higher = more opaque thus more absorbant, and vice-versa
minimum_absorption = 10.0 #Cut-off for loop of absorbance inside cell, lower means more loops occur thus higher absorbance

#ABSORBANCE LOOPING:
area = 1.0 #control for area to be constant

#################################################################################################################################################################################################################################################################
## WORKSPACE
#################################################################################################################################################################################################################################################################

avg_scan_abs = [] #stores the average absorbance of all the heights scanned per ellipse
avg_dimension_abs = [] #stores the average absorbance of all the dimensions scanned per ellipse,
                        #^averaged in each ellise from the height and rotation
dimensions = [] #stores dimensions for x-axis plotting
hit_counter = [] #counts total number of hits per dimensions

k = 1.0
while k <= 10.0:
    print('\n', k)
    #Find dimensions of ellipse according to constant area
    a = ((k*area)/np.pi)**0.5
    b = a/k

    ax = a
    ay = b
    
    dimen = a/b
    dimensions.append(dimen)
    no_hits = 0.0
    
    rotate_abs = [] #stores the percentage absorbances in each "rotation" in steps in accordance with deg
    deg = 0.0
    while deg <= 360.0:
        
        scan_abs = [] #stores the percentage absorbances in each "height" in steps in accordance with h
        
        h = -2.0
        while h <= 2.0:
            
            ##VARIABLE PARAMETERS
            ray = ((10.0, h) , (-10.0, h)) #Ray = ((ray vector) , (origin vector))
            phi = degree2rad(deg)
            ellipse_params = (ax, ay, phi)
                
            ##ARRAYS
            normal_vecs = [] #Stores all vectors for the normals
            trans_rays = [] #Stores all vectors for rays trasnmitted and scattered out by the ellipse
            refl_intensities = [] #Stores all intensities in tuples of p, s polarization for rays reflected within the ellipse
            trans_intensities = [] #Stores all intensities for rays trasnmitted and scattered out by the ellipse
            refl_ray_vecs = [] #Stores all vectors for rays reflected within the ellipse
            unit_refl_ray_vecs = [] #Stores all unit vectors for rays reflected within the ellipse
            intersects = [] #Stores all intersects by the ray as it reflects within the ellipse
            refract = True
            try: #checks to see if error occurs due to ray not hitting the ellipse
                ##1st INSTANCE:
                ni_point, unit_ni, unit_ray, unit_vt, vt_origin, vr, mag_vr, unit_vr, \
                    init_flux_p_ref, init_flux_s_ref, init_flux_p_trans, init_flux_s_trans \
                        = first_instance(inside_ri, outside_ri, incidence_flux_p, incidence_flux_s) #first instance function
            except:
                ni_point = 0.0
            
            if ni_point != 0.0: #error handling to skip this section of code if the ray does not hit ellipse
                xni, yni = ni_point #intersection points for first instance
                init_intensity = (init_flux_p_trans, init_flux_s_trans) #flux intensity of rays refracted into ellipse
                trans_intensity = (init_flux_p_ref + init_flux_s_ref) 
                    #^sum of p,s polarizations of rays trasnmitted and scattered out by the ellipse
                
                intersects.append(vt_origin) #appends first intersect to array
                unit_refl_ray_vecs.append(unit_vt) 
                normal_vecs.append(unit_ni) #appends unit vector of normal at first intersection
                trans_rays.append(unit_vr) #appends unit vector of first ray refracted into ellipse
                trans_intensities.append(trans_intensity)
                
                ##SUBSEQUENT INSTANCES:
                new_refl_vec = unit_vt * mag_vt #vector for first ray refracted into ellipse for use within loop
                new_refl_origin = vt_origin #first origin/point of intersection on ellipse
                ray_intensity = init_flux_p_trans + init_flux_s_trans #sets up initial intensity refracted into ellipse for while loop
                
                while ray_intensity > minimum_absorption:
                    #Set-up ray
                    new_refl = (new_refl_vec, new_refl_origin)
                    
                    try:
                        #Calculate new intersects, refraction, and transmission vectors
                        vti_points, unit_nvti, nri, mag_nri, unit_nri, unit_vt2, vt_origin2,\
                            sequ_flux_p_ref, sequ_flux_s_ref, sequ_flux_p_trans, sequ_flux_s_trans\
                                = subs_instances(inside_ri, outside_ri, new_refl[0], new_refl_origin,\
                                                 init_flux_p_trans, init_flux_s_trans)
                    except:
                        vti_points = 0.0
                    
                    if vti_points != 0.0:
                        xvti, yvti = vt_origin2
                    
                        #Calculate ray intensities and attenuation
                        trans_intensity = (sequ_flux_p_trans + sequ_flux_s_trans)
                        mag_newray = np.sqrt(((abs(xvti-xni))**2) + ((abs(yvti-yni)**2)))
                        sequ_intensity = (sequ_flux_p_ref, sequ_flux_s_ref)
                        refl_p_intensity = sequ_flux_p_ref * (np.exp(-mu * mag_newray))
                        refl_s_intensity = sequ_flux_s_ref * (np.exp(-mu * mag_newray))
                        refl_intensity = refl_p_intensity, refl_s_intensity
                        
                        #Add vectors to arrays:
                        normal_vecs.append(unit_nvti)
                        trans_rays.append(unit_vt2)
                        refl_intensities.append(refl_intensity)
                        trans_intensities.append(trans_intensity) #Transmitted out intensity
                        refl_ray_vecs.append(new_refl[0]) #Internal reflected ray vectors
                        unit_refl_ray_vecs.append(unit_nri) #Internal normal vectors
                        intersects.append(vt_origin2) #Intersection points
                        
                        #Reset variables:
                        unit_new_refl = unit_nri
                        new_refl_vec = unit_new_refl * mag_newray
                        xni = vti_points [0]
                        yni = vti_points [1]
                        new_refl_origin = xni, yni
                        init_flux_p_trans = refl_p_intensity
                        init_flux_s_trans = refl_s_intensity
                        
                        ray_intensity = sum(refl_intensity) #reset for loop to new intensity 
                        
                    else:
                        refract = False
                        break
                    
                if refract == True:
                    ## CALCULATE TOTAL ABSORBANCE 
                    intensity_absorbance = (incidence_flux_p + incidence_flux_s) - sum(trans_intensities) - ray_intensity
                    percent_absorbance = (intensity_absorbance / (incidence_flux_p + incidence_flux_s)) *100.0 
                        #^absorbance as a percentage of the original ray in intensity
                    if math.isnan(percent_absorbance):
                        percent_absorbance = 0.0
                        scan_abs.append(percent_absorbance)
                        print(k, deg, h, percent_absorbance, 'nan')
                    else:
                        scan_abs.append(percent_absorbance)
                        print(k, deg, h, percent_absorbance)
                else:
                    percent_absorbance = 0.0 #if ray does not refract into ellipse then no absorbance can occur
                    scan_abs.append(percent_absorbance)
                    print(k, deg, h, percent_absorbance, 'REFRACT ERROR')                
            else:
                percent_absorbance = 0.0 #if ray does not hit ellipse then no absorbance can occur
                scan_abs.append(percent_absorbance)
                no_hits += 1.0
                print(k, deg, h, percent_absorbance, 'NO HIT')
        
            h += 0.1 #steps of 1 in 100 from y value of -2.0 to 2.0
        
        avg_scan_abs.append(scan_abs)
        average_scan_abs = sum(scan_abs) / len(scan_abs)
        avg_scan_abs.append(average_scan_abs)
    
        rotate_abs.append(average_scan_abs)
        deg += 1.0 #steps of 1 in 360 so for every degree in a full rotation
    
    hit_counter.append(no_hits)
    average_abs = sum(rotate_abs) / len(rotate_abs)
    avg_dimension_abs.append(average_abs)
    print(hit_counter)
    print(avg_dimension_abs)
    k += 1.0

#################################################################################################################################################################################################################################################################
# PLOTS
#################################################################################################################################################################################################################################################################

# Plot the scatter graph
plt.figure(figsize=(10, 6))
plt.scatter(dimensions, avg_dimension_abs, color='black', marker = 'o', s=17)
plt.xlabel('a/b', fontname = 'serif')
plt.ylabel('Abs (%)', fontname = 'serif')
plt.title('Absorption for a/b at 1.5 RIU (averaged from height and rotation)', fontname = 'serif')
plt.grid(True)
plt.show()






