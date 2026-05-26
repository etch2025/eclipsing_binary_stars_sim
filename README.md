**Eclipsing Binary Stars Light Curve Simulator**  

**MATLAB Version: LC4_5.m  
Python Version: LC4_5.py  
Desmos Version: https://www.desmos.com/calculator/hvjfjnlnyz**

<img width="6321" height="4925" alt="image" src="https://github.com/user-attachments/assets/fb5880f9-1a8d-4bfc-9f83-3fb4001799f7" />

Major Assumptions:
- Uniform Luminosity Distribution (No limb darkening, luminosity/projected area = const)
- Circular Orbit (e = 0)
- Celestial bodies m1 and m2 is perfectly spherical with projected areas being perfectly circular

Outputs:
- Light Curve and Orbital Diagram
- Eclipse Luminosity and Duration
- Points of Contact
- Impact Parameter
- L1 Point (Roche Lobe Overflow)

Input Variable Units:
- Mass (m1, m2): Solar Masses
- Radii (r1, r2): Solar Radii
- Luminosities (L1, L2): Solar Luminosities
- Orbital Period (P): Days
- Orbital Inclination (i): Degrees

Boundary Conditions:
- r1 > r2
- P > P_minimum = sqrt((((r1+r2)*R_Sol)/AU)^3/(m1+m2)) * yr;
-   At P_minimum, the orbital semi-major axis becomes equal to r1 + r2, forming a contact binary star
- i > acos(((r1 + r2)*R_Sol)/sma) * (180/pi);
-   Minimum Orbital Inclination for an eclipse in Deg
