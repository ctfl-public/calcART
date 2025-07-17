from sparta import sparta 
import os
from constants import SIGMA
import sys
import numpy as np

# output directory name
# mydir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
mydir = "dump"
if not os.path.exists(mydir):
    os.makedirs(mydir)

def calc_dqrad(kappa, sigma_sca, T, outputName, limits=[0, 1], \
               size=1 , nRays=1000, SF='LA', g1=0, \
                nonhomogeneous=False, \
                    in_rad=0):
    """
    Calculate radiative heat flux (Dqrad).

    The medium is confined between ylo and yhi in y direction;
    ylo is cold black surface, yhi is black surface with {in_rad} radiation input.

    Args:
        kappa (float): Absorption coefficient.
        sigma_sca (float): Scattering coefficient.
        T (float or str): Temperature of the medium, can be a constant value or a file.
            file should contain x, T if nonhomogeneous is False,
            or x, T, kappa, sigma_sca if nonhomogeneous is True.
        outputName (str): Name of output file (includes extension).
        limits (list): [ylo, yhi] limits of the grid in y direction.
        size (int): Size of the grid along y axis (default is 1).
        nRays (int): Number of rays to be emitted from each cell (default is 1000).
        SF (str): Scattering function, 'LA' for Linear Anisotropic, 'HG' for Henyey-Greenstein.
        g1 (float): Anisotropy factor for LA or HG.
        nonhomogeneous (bool): If True, T_profile is a file with x, T, kappa, and sigma_sca values. 
        in_rad (float): Radiation input value.

    Returns:
        None. Writes output to a file in the mydir directory.
        The output file contains yc and Dqrad values.
    """
    
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile):
        print(f"Case already exists: {outfile}")
        return
    
    T_profile = T
    if (type(T) is str):
        T = 0

    # change to your machine name used to compile sparta: "serial", "serial_debug", ..
    spa = sparta("serial") 

    spa.command("seed 8887435")
    spa.command("units si")
    spa.command("dimension 3")
    spa.command("boundary ss ss ss")

    if (SF == 'LA'): # Linear Anisotropic
        spa.command("global radiation RMCRT pathlength 1 \
        sigma_sca " + str(sigma_sca) + " g1 " + str(g1) +" kappa " + str(kappa) +  " T " + str(T))
    elif (SF == 'HG'):  # Henyey-Greenstein
        spa.command("global radiation RMCRT pathlength 1 \
        sigma_sca " + str(sigma_sca) + " HG g1 "+ str(g1) +" kappa " + str(kappa) +  " T " + str(T))
    else:
        sys.exit(SF+" is not recognized.")

    spa.command("photon_continuous none constant")
    spa.command("create_box -1 1 " + int2str(limits) + " -1 1")
    spa.command("create_grid 1 " + str(size) + " 1 block * * *")

    if (type(T_profile) is str):
        # variable kappa and sigma, file contains x, T, kappa, and sigma
        if nonhomogeneous: 
            spa.command("global MR " + T_profile)
        # constant kappa and sigma, file contains only x and T
        else: 
            spa.command("global T MR " + T_profile)

    spa.command("timestep 0.5")
    spa.command("balance_grid rcb part")

    # black surface
    spa.command("surf_collide 1 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp "+ str((in_rad/SIGMA)**0.25))
    # cold black surface
    spa.command("surf_collide 2 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp 0")
    # cold mirror surface
    spa.command("surf_collide 3 radiationboundary epsilon 0 alpha 0 rho_s 1 rho_d 0 temp 0")

    # hack to delete photon
    spa.command("surf_react 1 global 1.0 0.0")

    # assign collide modules
    spa.command("bound_modify yhi collide 1 react 1")
    spa.command("bound_modify ylo collide 2 react 1")
    spa.command("bound_modify xlo xhi zlo zhi collide 3 react 1")

    spa.command("create_photons npc " + str(nRays) + " xzCenter")

    spa.command("stats 5000")
    spa.command("run 1")
    spa.command("dump 1 grid all 999999999 " + outfile + " yc Dqrad")
    spa.command("run 200000000")

    spa.close()

# convert list of integers/floats to space separated string
def int2str(intList):
    return str(intList).replace(',','').replace('[','').replace(']','')

# return emissivity value from sparta output file
def read_dqrad(fileName):
    filePath = os.path.join(mydir,fileName)
    yc = []
    dqrad = []
    with open(filePath, 'r') as f:
        lines = f.readlines()

        # skip first 9 lines
        lines = lines[9:]

        # read corresponding data
        for line in lines:
            vals = [float(i) for i in line.split(" ") if i.strip()]
            
            yc.append(vals[0])
            dqrad.append(-vals[1])

    return yc, dqrad


def calc_ref(thickness, ext, omega,
                nRays=200000, SF='LA', A1=0, 
                g1=0, g2=0, f1=0, f2=0,
                theta=None, theta_in=None, div_angle=None,
                machine = "serial", 
                cmdargs = ["-screen","none"]):
    
    # ensure old version using A1 is working if not replaced with g1
    g1 = A1 if (A1 != 0) else g1

    if (div_angle and not is_num(theta_in)): raise ValueError("Incident direction is not defined.")

    sigma_sca = ext * omega
    kappa = ext * (1-omega)

    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.0f}.ref")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.0f}.ref")

    if os.path.exists(outfile):
        print(f"Case already exists: {outfile}")
    else:
        spa = sparta(machine, cmdargs) 

        spa.command("seed 8887435")
        spa.command("units si")
        spa.command("dimension 3")
        spa.command("boundary ss ss ss")

        if (SF == 'LA'): # Linear Anisotropic
            spa.command("global radiation RMCRT pathlength 1 \
            sigma_sca " + str(sigma_sca) + " g1 " + str(g1) +" kappa " + str(kappa) +  " T 0")
        elif (SF == 'HG'):  # Henyey-Greenstein
            spa.command("global radiation RMCRT pathlength 1 \
            sigma_sca " + str(sigma_sca) + " HG g1 "+ str(g1) +" kappa " + str(kappa) +  " T 0")
        elif (SF == 'NC'):  # combination of HG, Nicolau et al. (1994)
            spa.command(f"global radiation RMCRT pathlength 1 \
            sigma_sca {sigma_sca} NC g1 {g1} g2 {g2} f1 {f1} f2 {f2} kappa {kappa} T 0")
        else:
            sys.exit(SF+" is not recognized.")

        spa.command("photon_continuous none constant")
        spa.command("create_box -1 1 0 "+ str(thickness)+ " -1 1")
        spa.command("create_grid 1 1 1 block * * *")

        spa.command("timestep 0.5")
        spa.command("balance_grid rcb part")

        # cold mirror surface
        spa.command("surf_collide 1 radiationboundary epsilon 0 alpha 0 rho_s 1 rho_d 0 temp 0")
        # cold black surface
        spa.command("surf_collide 2 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp 0")
        # black surface (unit emission)
        if (div_angle):
            a1 = np.sin(theta_in)
            a2 = -np.cos(theta_in)
            a3 = 0
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25} "+
                        f" collimated {a1} {a2} {a3} div {div_angle}")
        elif (is_num(theta_in)):
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}"+ \
                        f" theta_in {theta_in}")
        else:
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 alpha 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}")

        # hack to delete photon
        spa.command("surf_react 1 global 1.0 0.0")

        # assign collide modules
        spa.command("bound_modify xlo xhi zlo zhi collide 1 react 1")
        spa.command("bound_modify ylo collide 2 react 1")
        spa.command("bound_modify yhi collide 3 react 1")

        #create_photons n 1000
        if (is_num(theta)):
            spa.command(f"fix 1 emit/face/photon yhi n {nRays} nevery 999999999 limit xzCenter theta {theta}")
        else:
            spa.command(f"fix 1 emit/face/photon yhi n {nRays} nevery 999999999 limit xzCenter")

        spa.command("variable epsilon equal 1")
        spa.command("variable T equal 0")

        spa.command("stats 5000")
        spa.command("run 1")

        spa.command(f"dump 1 grid all 999999999 {outfile} id xc yc zc I countEmitted")
        spa.command("run 200000000")

        spa.close()

    # dir-dir (div_angle method)
    if (div_angle and is_num(theta)): 
        return read_prop(outfile) / (2*np.pi*(1-np.cos(div_angle)))
    # dir-dir (cos method)
    elif (not div_angle and is_num(theta_in) and is_num(theta)): 
        return read_prop(outfile) #* np.cos(theta_in) 
    # hem-dir
    elif (not div_angle and not is_num(theta_in) and is_num(theta)): 
        return read_prop(outfile) 
    # dir-hem (div_angle method)
    elif (div_angle and not is_num(theta)): 
        # print("Warning: Consider integrating over all out directions to get dir-hem property.")
        return read_prop(outfile) / (2*(1-np.cos(div_angle)))
    # dir-hem (cos method)
    elif (not div_angle and is_num(theta_in) and not is_num(theta)): 
        # print("Warning: Consider integrating over all out directions to get dir-hem property.")
        return read_prop(outfile) * np.pi
    # hem-hem
    elif (not div_angle and not is_num(theta_in) and not is_num(theta)): 
        return read_prop(outfile) * np.pi


def is_num(var):
    return isinstance(var, (int, float))

def read_prop(fileName):
    with open(fileName, 'r') as f:
        # skip first 9 lines
        for _ in range(9):
            next(f)

        # read corresponding data
        line = f.readline()
        vals = [float(i) for i in line.split(" ") if i.strip()]
        
        id = vals[0]
        xc = vals[1]
        yc = vals[2]
        zc = vals[3]
        out = vals[4]
        count = vals[5]

    return out


# # rum main application
# if __name__ == "__main__":
#     main()