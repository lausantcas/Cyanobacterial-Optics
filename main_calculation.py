# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:08:49 2024

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

#ELLIPSE PARAMETERS:
area = 1.0 #internal area of ellispe
k = 4/3 #proportional parameter for how much larger a (ax) is larger than b (ay)
a = ((k*area)/np.pi)**0.5
b = a/k

ax = a #1/2 total width of ellipse along semi-major axis irrespective of rotation
ay = b #1/2 total height of ellipse semi-major axis irrespective of rotation

#RAY PARAMETERS:
h = b #height of ellipse
ray = ((10.0, h) , (-10.0, h)) #Ray = ((ray vector) , (origin vector))
    
phi = degree2rad(135.0) #Angle of rotation of ellipse
ellipse_params = (ax, ay, phi)

#TANGENT & NORMAL PARAMETERS:
delta_t = 0.001 #Small change in t to calculate the tangent in normal calculation
tol = 0.0001 #Angle tolerance
ang_res = 10000.0 #Angle resolution
normal_params = (delta_t, tol, int(ang_res))

#CONTROL PARAMETERS:
t_range_param = 10000.0
plot_range = 10000.0 #Avoids errors with refraction intersection calculations
control_params = (int(t_range_param), plot_range)

#REFRACTION PARAMETERS:
outside_ri = 1.0 #Refraction index of media outside ellipse
inside_ri = 1.4 ##Refraction index of media inside ellipse
mag_vt = 100000.0 #Gives magnitude for calculating to initial ray inside ellipse

#RAY CODING:
incidence_flux_p = 500.0 #p polarization flux
incidence_flux_s = 500.0 #s polarization flux
mu = 1.0 #Absorbance parameter by the ellipse: higher = more opaque thus more absorbant, and vice-versa
minimum_absorption = 0.001 #Cut-off for loop of absorbance inside cell, lower means more loops occur thus higher absorbance

#################################################################################################################################################################################################################################################################
## WORKSPACE & DEBUGGING
#################################################################################################################################################################################################################################################################

##ARRAYS
normal_vecs = [] #Stores all vectors for the normals
trans_rays = [] #Stores all vectors for rays trasnmitted and scattered out by the ellipse
refl_intensities = [] #Stores all intensities in tuples of p, s polarization for rays reflected within the ellipse
trans_intensities = [] #Stores all intensities for rays trasnmitted and scattered out by the ellipse
refl_ray_vecs = [] #Stores all vectors for rays reflected within the ellipse
unit_refl_ray_vecs = [] #Stores all unit vectors for rays reflected within the ellipse
intersects = [] #Stores all intersects by the ray as it reflects within the ellipse
refract = True
try:
    ##1st INSTANCE:
    ni_point, unit_ni, unit_ray, unit_vt, vt_origin, vr, mag_vr, unit_vr, \
        init_flux_p_ref, init_flux_s_ref, init_flux_p_trans, init_flux_s_trans \
            = first_instance(inside_ri, outside_ri, incidence_flux_p, incidence_flux_s) #first instance function
except:
    ni_point = 0.0

if ni_point != 0.0:
    xni, yni = ni_point #intersection points for first instance
    init_intensity = (init_flux_p_trans, init_flux_s_trans) #flux intensity of rays refracted into ellipse
    trans_intensity = (init_flux_p_ref + init_flux_s_ref) #sum of p,s polarizations of rays trasnmitted and scattered out by the ellipse
    
    intersects.append(vt_origin) #appends first intersect to array
    unit_refl_ray_vecs.append(unit_vt) #appends 
    normal_vecs.append(unit_ni) #appends unit vector of normal at first intersection
    trans_rays.append(unit_vr) #appends unit vector of first ray refracted into ellipse
    trans_intensities.append(trans_intensity)
    
    ##SUBSEQUENT INSTANCES:
    new_refl_vec = unit_vt * mag_vt #vector for first ray refracted into ellipse for use within loop
    new_refl_origin = vt_origin #first origin/point of intersection on ellipse
    ray_intensity = init_flux_p_trans + init_flux_s_trans #sets up initial intensity refracted into ellipse for while loop
    
    # refract = True
    
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
        print('Total absorbance by cell:', round(percent_absorbance, 3),'%')
        if math.isnan(percent_absorbance):
            percent_absorbance = 0.0
            print('Total absorbance by cell:', round(percent_absorbance, 3),'%')
    else:
        percent_absorbance = 0.0 #if ray does not refract into ellipse then no absorbance can occur
        print('Total absorbance by cell:', round(percent_absorbance, 3),'%')
        print('The ray does not refract into the ellipse despite intercepting within error margin.\n')
else:
    percent_absorbance = 0.0 #if ray does not hit ellipse then no absorbance can occur
    print('Total absorbance by cell:', round(percent_absorbance, 3),'%')
    print('The ray does not intercept the ellipse.')


#################################################################################################################################################################################################################################################################
# PLOTS
#################################################################################################################################################################################################################################################################

#PLOTTING ELLIPSE
t_range = np.linspace(0.0, 2.0*np.pi, control_params[0]) #t is the angle in radians of the ellipse
x_ellipse_points = []
y_ellipse_points = []
for t in t_range:
    point = ellipse_polar(t, ax, ay, phi)
    x_ellipse_points.append(point[0]) #Creates list of x coordinates on ellipse
    y_ellipse_points.append(point[1]) #Creates list of y coordinates on ellipse
plt.plot(x_ellipse_points, y_ellipse_points, c="g", alpha = 0.55, label = 'Ellipse/Circle')
mag_ray = np.sqrt((ray[0][0]**2) + (ray[0][1]**2)) 

if percent_absorbance == 0.0 :
    #PLOTTING RAY:
    ray_origin = ray[1]
    ray_vector = ray[0]
    x_ray_points = [ray_origin[0], ray_vector[0]]
    y_ray_points = [ray_origin[1], ray_vector[1]]
    plt.plot(x_ray_points, y_ray_points, c="blue", alpha = 0.75, label = 'Ray')
elif refract == False or math.isnan(percent_absorbance):
    #PLOTTING RAY:
    ray_origin = ray[1]
    ray_vector = ray[0]
    x_ray_points = [ray_origin[0], intersects[0][0]]
    y_ray_points = [ray_origin[1], intersects[0][1]]
    plt.plot(x_ray_points, y_ray_points, c="blue", alpha = 0.75, label = 'Ray')
    
    #PLOTTING NORMAL:
    x_normal_points = [intersects[0][0], normal_vecs[0][0]]
    y_normal_points = [intersects[0][1], normal_vecs[0][1]]
    plt.plot(x_normal_points, y_normal_points, '--', c='black', alpha = 0.5, label = 'Normal')
  
    # #PLOTTING INTERNAL REFLECTION VECTORS:
    i = 0
    x_vti_points = [intersects[i][0], unit_refl_ray_vecs[0][0]]
    y_vti_points = [intersects[i][1], unit_refl_ray_vecs[0][1]]
    plt.plot(x_vti_points, y_vti_points, c="orange", alpha = 0.75, label = "Reflected" if i == 0 else "")
    
    #PLOTTING TRANSMITTANCE VECTORS:
    k = 0
    target = (0.0, 0.0)
    if tuple(trans_rays[k]) == target: #checks for total internal reflection as this means no transmittance ray plotted
        a = 1.0
    else:
        x_vt2_points = [intersects[k][0], mag_ray*trans_rays[k][0]]
        y_vt2_points = [intersects[k][1], mag_ray*trans_rays[k][1]]
    plt.plot(x_vt2_points, y_vt2_points, c="red", alpha = 0.55, label = "Transmitted" if k == 0 else "")
else:
    #PLOTTING RAY:
    ray_origin = ray[1]
    ray_vector = ray[0]
    x_ray_points = [ray_origin[0], intersects[0][0]]
    y_ray_points = [ray_origin[1], intersects[0][1]]
    plt.plot(x_ray_points, y_ray_points, c="blue", alpha = 0.75, label = 'Ray')
    
    #PLOTTING NORMAL:
    x_normal_points = [intersects[0][0], plot_range*unit_ni[0]]
    y_normal_points = [intersects[0][1], plot_range*unit_ni[1]]
    plt.plot(x_normal_points, y_normal_points, '--', c='black', alpha = 0.5, label = 'Normal')
    
    # #PLOTTING INTERNAL REFLECTION VECTORS:
    plot_loop = len(intersects)
    i = 0
    while plot_loop > 1:
            #Reflected vectors plotted
        x_vti_points = [intersects[i][0], intersects[i+1][0]]
        y_vti_points = [intersects[i][1], intersects[i+1][1]]
        plt.plot(x_vti_points, y_vti_points, c="orange", alpha = 0.75, label = "Reflected" if i == 0 else "")
            #Normal vectors plotted
        x_vtn_points = [intersects[i+1][0], plot_range*normal_vecs[i+1][0]]
        y_vtn_points = [intersects[i+1][1], plot_range*normal_vecs[i+1][1]]
        plt.plot(x_vtn_points, y_vtn_points, '--', c='black' , alpha = 0.5)
        i += 1
        plot_loop -= 1
        
    #PLOTTING TRANSMITTANCE VECTORS:
    loop1 = len(intersects)+1
    k = 0
    target = (0.0, 0.0)
    while loop1 > 1:
        if tuple(trans_rays[k]) == target: #checks for total internal reflection as this means no transmittance ray plotted
            a = 1.0
        else:
            x_vt2_points = [intersects[k][0], mag_ray*trans_rays[k][0]]
            y_vt2_points = [intersects[k][1], mag_ray*trans_rays[k][1]]
        plt.plot(x_vt2_points, y_vt2_points, c="red", alpha = 0.55, label = "Transmitted" if k == 0 else "")
        k += 1
        loop1 -= 1
    

#OTHER PLOTTING:
plt.xlim(-max(ax*1.15,ay*1.15),max(ax*1.15,ay*1.15))
plt.ylim(-max(ax*1.15,ay*1.15),max(ax*1.15,ay*1.15))
plt.legend(loc='lower left',
          prop=dict(family='serif', size ='8'),
          bbox_to_anchor=(1, 0))
plt.gca().set_aspect('equal')
plt.xlabel('Cell length (micrometers)', fontname = 'serif')
plt.ylabel('Cell height (micrometers)', fontname = 'serif')
plt.title('Ellipse', fontname = 'serif')
plt.show()


