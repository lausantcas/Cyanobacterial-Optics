# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:45:23 2024

@author: Laura SC
"""

import numpy as np
from intersect import intersection

def degree2rad(deg):
    return((deg/360.0)*2.0*np.pi)

def rad2degree(rad):
    return(rad*(180.0/np.pi))

def ellipse_polar(t, ax, ay, phi):
    x_prime = (ax * np.cos(t) * np.cos(phi)) - (ay * np.sin(t) * np.sin(phi)) #Rotated x coordinates
    y_prime = (ax * np.cos(t) * np.sin(phi)) + (ay * np.sin(t) * np.cos(phi)) #Rotated y coordinates
    return(x_prime, y_prime)

"""############################################################################
Function to calculate the normal unit vector, ni, at the point of intersection 
between an optical ray and an arbitrary rotated ellipse centred on the origin.

Inputs: 
    (1) ray = (vec, origin), the optical ray specified as a unit vector and 
                             point of origin
               vec = (vec_x, vec_y)
               origin = (orig_x, orig_y)
            
    (2) ellipse = (ax, ay, phi)
        ax = 1/2 width of ellipse
        ay = 1/2 height of ellipse
        phi = angle between centre of ellipse and intersection point
    
    (3) normal = (delta_t, tol, int(ang_res))
        delta_t = small change in t
        ang_res = fraction of 2pi
        tol = tolerance        

    (4) params = (int(t_range_param), plot_range)
        t_range_param = number of points in ellipse
        plot_range = Distance between points for generation of line for intersect
                     and plotting calculation

Function:
    (1) Unpack params
    (2) Generate ellipse
    (3) Generate ray
    (4) Find intersection 
    (5) Select for intersection coordinates depending on number of intersection
        points
    (6) Calculate vectors for tangent and normal
        a. Case 1: if a circle
        b. Case 2: if intersect lying upon semi-major or semi-minor axes
        c. Case 3: if none of the above

Outputs: 
    (1) The normal unit vector, ni
    (2) The tangent unit vector, Ti
    (3) The x-coordinate of the intersection point of interest, xi
    (4) The y-coordinate of the intersection point of interest, yi

############################################################################"""

def refract_normal(ray, ellipse, normal, params):
    
    #(1) UNPACK PARAMS:
    ax, ay, phi = ellipse
    delta_t, tolerance, ang_res = normal
    t_range_param, plot_range = params
    
    #(2) GENERATE ELLIPSE
    t_range = np.linspace(0.0, 2.0*np.pi, t_range_param) #t is the angle in radians of the ellipse
    x_ellipse_points = []
    y_ellipse_points = []
    for t in t_range:
        point = ellipse_polar(t, ax, ay, phi)
        x_ellipse_points.append(point[0]) #Creates list of x coordinates on ellipse
        y_ellipse_points.append(point[1]) #Creates list of y coordinates on ellipse

    #(3) GENERATE RAY:
    ray_vec = ray[0]
    origin = ray[1]
    x_ray_points = [origin[0], ray_vec[0]]
    y_ray_points = [origin[1], ray_vec[1]]
    
    #(4) FIND INTERSECTION COORDINATES:    
    xi, yi = intersection(x_ray_points, y_ray_points, x_ellipse_points, y_ellipse_points)
    if len(xi) == 1: #Selects depending on number of intersection points
        i1 = [xi[0], yi[0]]
        xi, yi = i1
    elif len(xi) == 2:
        i1 = [xi[0], yi[0]]
        i2 = [xi[1], yi[1]]
        mag_i1 = np.sqrt((i1[0]-origin[0]) * (i1[0]-origin[0]) + (i1[1]-origin[1]) * (i1[1]-origin[1]))
        mag_i2 = np.sqrt((i2[0]-origin[0]) * (i2[0]-origin[0]) + (i2[1]-origin[1]) * (i2[1]-origin[1]))
        #(5) SELECT FOR CORRECT INTERSECTION COORDINATES
        if mag_i1 < mag_i2:
            xi, yi = i2
        else:
            xi, yi = i1
    elif len(xi) == 0:
        print('ERROR! No refraction intersections found - ray is tangential to surface, no angle of refraction found.')
        return(0)
    else:
        print('ERROR! Normal refraction intersection error!')  
        return(0)

    #(6) CALCULATE TANGENT & NORMAL:       
    mag_ri = np.sqrt((xi*xi) + (yi*yi))
    #(6.a) Case 1
    if ax == ay: #ie. If a circle
        ni = (xi/mag_ri , yi/mag_ri)
        Ti = (-(yi/mag_ri) , (xi/mag_ri))
    #(6.b)Case 2
    elif ax != ay and (abs(mag_ri - ax) < tolerance or abs(mag_ri - ay) < tolerance):
        ni = (xi/mag_ri , yi/mag_ri)
        Ti = (-(yi/mag_ri) , (xi/mag_ri))
    #(6.c) Case 3 
    else:
        cos2t = ((2.0*mag_ri*mag_ri) - ((ax*ax) + (ay*ay))) / ((ax*ax) - (ay*ay))
        t1 = 0.5 * np.arccos(cos2t)
        t2 = np.pi - t1
        t3 = np.pi + t1
        t4 = (2.0*np.pi) - t1
        #For non-rotated coords
        x_non = (xi * np.cos(phi)) + (yi * np.sin(phi))
        y_non = (yi * np.cos(phi)) - (xi * np.sin(phi))
        if x_non > 0 and y_non >= 0:
            ti = t1
        elif x_non <= 0 and y_non > 0:
            ti = t2
        elif x_non < 0 and y_non <= 0:
            ti = t3
        elif x_non >= 0.0 and y_non < 0.0:
            ti = t4
        else:
            print('ERROR! n_i is undefined for (x_i,y_i) = (0,0) (REFRACTION)\n')
            return('ERROR \n')
        #Finding the normal
        plus_xt, plus_yt = ellipse_polar((ti + (delta_t/2.0)), ax, ay, phi)
        minus_xt, minus_yt = ellipse_polar((ti - (delta_t/2.0)), ax, ay, phi)
        mag_T = np.sqrt(((plus_xt - minus_xt)**2) + ((plus_yt - minus_yt)**2))
        Tx = (plus_xt - minus_xt)/mag_T
        Ty = (plus_yt - minus_yt)/mag_T
        Ti = (Tx, Ty)
        if x_non > 0 and y_non >= 0:
            ni = (-Ty , Tx)
        elif x_non <= 0 and y_non > 0:
            ni = (Ty , -Tx)
        elif x_non < 0 and y_non <= 0:
            ni = (-Ty , Tx)
        elif x_non >= 0.0 and y_non < 0.0:
            ni = (Ty , -Tx)
        else:
           print('ERROR! n_i is undefined for (xi , yi) = (0,0) (REFRACTION)\n')
           return('ERROR \n')
       
    return(ni, Ti, xi, yi)


