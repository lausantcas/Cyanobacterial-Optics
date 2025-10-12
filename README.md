# Cyanobacterial Optics
Cyanobacteria have been found to elongate under low irradiance conditions. The model found in this repository treats a single cyanobacterial cell as a 2D circle or ellipse, of which the elongation can be controlled in a set parameter (k), and calculates the total absorbance by a cell when a collimated light ray intersects it. A second script also averages results across all rotation (0–360°) and vertical offsets to approximate diffuse-light conditions and produces comparison plots versus ellipticity (ratio ax/ay).

The model tracks a collimated incidence light ray as it intersects the boundary of the cell, reflecting and refracting according to Fresnel’s equations and Snell’s Law, in vector form. The ray is also attenuated inside the cell via the Beer–Lambert law until it once again intersects the internal surface of the cell. The process of reflection and refraction then repeats, taking into account the change in refractive indices of each medium, and this new point of intersection is the new point of origin for the next iteration. This is then overall looped until the minimum absorption condition is met, which is when the intensity flux encoded on the ray falls below this set minimum absorption parameter. Total absorption by the cell is computed as the fraction of the incident intensity that is lost inside the cell across multiple internal passes. 

## What the code does
### 1. Define cell geometry
The cell is a rotated ellipse with semi-axes (ax, ay), and rotation phi. The parametric form is used to locate point on the boundary and to handle rotation cleanly.

### 2. Cast a light ray
A directed ray with origin r0​ and unit direction vi is traced to the first boundary intersection. The correct intersection is selected by magnitude and direction filters.

### 3. Compute the surface normal
The unit tangent at the hit point is estimated, the the outward normal is chosen based on quadrant logic to ensure it points away from the cell.

### 4. Apply interface optics
Using the Snell's Law in vector form to compute the transmitted vector; using Fresnel equations to split the incident ray intensity into reflected and transmitted s and p components, then summed for total flux.

### 5. Propagate inside the cell
Between hits, intensity attenuates by Beer-Lambert Law with the path length L. The loop continues: at each interior boundary, computes new reflection and transmission vectors to update the intensity fluxes. This loop ends when the internal flux falls below a minimum flux control.

### 6. Aggregate absorption
Absorption = (incident intensity − sum of all intensities transmitted to the outside) / (incident intensity). The single-run script plots the ray paths and reports % absorption.

### 7. Diffuse-light averaging
Repeat steps 1–6 across all rotations and ray heights, average per geometry, and plot absorption vs ellipticity; this reproduces the rising trend with increasing ellipticity ratio of ax/ay.

## Results found
Circles vs ellipses: circles do not exhibit total internal reflection (TIR) no matter the height of the ray due to the geometry; ellipses do at certain rotations and ray heights, which helps retain light internally and increases absorption.

<img width="500" height="222" alt="image" src="https://github.com/user-attachments/assets/4772cbcc-91a8-4601-a729-c34c0606ee36" />


Ellipticity trend: holding area constant (1 μm²), average absorption increases with ellipticity under diffuse-light conditions. Below can be seen some preliminary results supporting this.

<img width="500" height="324" alt="image" src="https://github.com/user-attachments/assets/e1c154c3-8955-43c4-93a2-ce4f2c168fcf" />


## What each code file correlates to
### normal_function.py and refract_normal_function.py
The normal_function.py script carries out step 3 from the above section, which computes the surface normal at the first intersection point where the ray intersects the cell and the point of intersection of interest is the first.

However, once the ray is looped and bounced around the inside of the cell, it is in interest to maintain the first intersection point as the origin coordinates, and so the normal calculation for the following intersection point means that the second intersection of this ray is the one of interest now. This is calculated in the script found within refract_normal_function.py

### main_calculation.py
This script carries out steps 1-6, calling upon the normal_function.py and refract_normal_function.py scripts for step 4.

### average_absorbances_loop_plotting.py
This script loops and averages using the main_calculation.py script to simulate diffuse-light conditions. It does this by looping the main_calculation.py file and changing the following parameters:
- The ellipticity ratio in steps predetermined in parameter k.
- The height at which a horizontal incident collimated light ray intersects the cell in constant steps within the loop, this is determined from a minimum start to a maximum, with height h.
- The rotation angle of the cell determined as phi, in constant steps determined by the user.

The absorption of each height and rotation are appended into arrays and averaged overall to find the average absorption for a set control area and ellipticity dimensions of each cell.

## More information
This project was developed for a dissertation, please contact for further details.
