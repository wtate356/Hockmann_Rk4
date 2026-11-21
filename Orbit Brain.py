import numpy as np

# parameters (using example satellite from AI)
ve = 3100                  #fuel exit velocity in m/s
mi = 3000                   #kg
massburn = 0.16             #kg/s rate at which mass is burnt
mu = 3.986004418e14         #GM in m^3 / s^2

delta_t = 10                #time step in seconds
initialRadius = 6671000     #initial circular orbit radius in m
finalRadius = 7000000       #destination circular orbit radius in m
initialVelocity = 7730      #initial velocity in m/s


def iterative_calculator(ri, vi, mi, a, delta_t):        #The brain, calculates radius change and tells when to turn on or off
    if a < finalRadius and ri < finalRadius:
        delta_v = ve * np.log(mi / (mi - massburn * delta_t))       #if spacecraft is below turns on for time interval, otherwise coasts (below if apogee distance < finalRadius)
        mi = mi - massburn * delta_t            #updates mass
    else:
        delta_v = 0

    a = abs((2/ri - ((vi + delta_v) ** 2)/mu) ** -1)         #semi-major axis of new ellipse
    ecc = 1- ri/a      #eccentricity
    M = np.sqrt(mu / a**3) * delta_t    #mean anomaly

    # solves for E which is an IVP
    Ei = M
    change = 10
    while abs(change) >= 1e-6:          #finds E (eccentric anomaly)
        change = (Ei - ecc*np.sin(Ei) - M) / (1 - ecc * np.cos(Ei))
        Ef = Ei - change
        Ei = Ef

    costheta = (np.cos(Ef) - ecc) / (1 - ecc * np.cos(Ef))      #finds cos(theta) where theta = true anomaly
    theta = np.arccos(costheta) * 180 / np.pi                   #finds angle
    rf = a * (1 - ecc**2) / (1 + ecc * costheta)                #finds radius after the time step where travels across "new ellipse"
    vf = np.sqrt(mu * (2 / rf - 1 / a))                         #finds velocity after time step at that radius (because velocity slows after travelling across ellipse)

    return float(rf), float(vf), float(mi), float(theta), float(a)  #outputs rf, vf, and mi to iterate upon, theta to graph it, a for conditioning

# sets up the while loop to have variables to call
inTransfer = True
radius = initialRadius
velocity  = initialVelocity
mass = mi
a = initialRadius
n = 0

listTime = [0]
listRadius = [initialRadius]
listVelocity = [initialVelocity]
listMass = [mi]
listAngles = [0]


while inTransfer:
    x = iterative_calculator(radius, velocity, mass, a, delta_t)
    radius = x[0]
    velocity = x[1]
    mass = x[2]
    a = x[4]
    n += 1

    if np.isnan(radius):    #stops loop if mass all burnt
        inTransfer = False
    if abs(radius - finalRadius) < 10 and abs(radius - listRadius[-1]) < 10:    #tolerances
        inTransfer = False
    
    listTime.append(n*delta_t)
    listRadius.append(radius)
    listVelocity.append(velocity)
    listMass.append(mass)
    listAngles.append(x[3] + listAngles[-1])

#print(f"Time array: {listTime}")
#print(f"Radii: {listRadius}")
#print(f"Velocity: {listVelocity}")
#print(f"Mass: {listMass}")
#print(f"Angles: {listAngles}")

