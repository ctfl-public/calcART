from sparta import sparta 
import os
from constants import SIGMA
import sys
import numpy as np
import warnings

# output directory name
# mydir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
mydir = os.path.join(os.getcwd(), "dump")
if not os.path.exists(mydir):
    os.makedirs(mydir)

def calc_dqrad(kappa, sigma_sca, T, outputName=None, limits=[0, 1], \
                size=1 , nRays=1000, SF='LA', g1=0, \
                nonhomogeneous=False, \
                in_rad=0, \
                machine = "serial", \
                cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux (Dqrad, W/m3).

    The medium is confined between ylo and yhi in y direction;
    ylo is cold black surface, yhi is black surface with {in_rad} radiation input.

    Args:
        kappa (float): Absorption coefficient.
        sigma_sca (float): Scattering coefficient.
        T (float or str): Temperature of the medium, can be a constant value or a file.
            file should contain x, T if nonhomogeneous is False,
            or x, T, kappa, sigma_sca if nonhomogeneous is True.
        outputName (str): Name of output file (includes extension). 
            if None, output is renamed to default name. (default is None).
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
    
    D = abs(limits[1] - limits[0])
    if not outputName:
        outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D:0.3f}-size{size}-nRays{nRays}-inrad{in_rad:0.1e}.dqrad"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"Case already exists and is not empty: {outfile}")
        return
    
    T_profile = T
    if (type(T) is str):
        T = 0

    # change to your machine name used to compile sparta: "serial", "serial_debug", ..
    spa = sparta(machine, cmdargs) 

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
    spa.command("surf_collide 1 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str((in_rad/SIGMA)**0.25))
    # cold black surface
    spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")
    # cold mirror surface
    spa.command("surf_collide 3 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")

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

def int2str(intList):
    """
    Convert list of integers/floats to space separated string.
    
    Args:
        intList (list): List of integers or floats to convert.
        
    Returns:
        str: Space-separated string representation of the input list.
    """
    return str(intList).replace(',','').replace('[','').replace(']','')

# return emissivity value from sparta output file
def read_dqrad(fileName):
    """   
    Parses the output file generated by calc_dqrad function to extract
    y-coordinates and corresponding radiative heat flux (Dqrad) values.
    
    Args:
        fileName (str): Name of the output file in the dump directory.
        
    Returns:
        tuple: A tuple containing two lists:
            - yc (list): Y-coordinate values
            - dqrad (list): Radiative heat flux values (negated from file)
                positive values indicate incoming radiation (heat gain),
                negative values indicate outgoing radiation (heat loss).
    """
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


def calc_trans(thickness, ext, omega,        
                nRays=200000, SF='LA', A1=0, 
                g1=0, g2=0, f1=0, f2=0,
                theta=None, theta_in=None, div_angle=None,
                machine = "serial",
                cmdargs = ["-screen","none"]):
    """
    Calculate transmisivity of a medium.
    
    for various angular configurations (directional-directional,
    hemispherical-directional, directional-hemispherical, hemispherical-hemispherical).
    
    Args:
        thickness (float): Thickness of the medium slab.
        ext (float): Extinction coefficient of the medium.
        omega (float): Scattering albedo (ratio of scattering to extinction).
        nRays (int): Number of rays from surface to simulate (default: 200000).
        SF (str): Scattering function type - 'LA' (Linear Anisotropic), 
                 'HG' (Henyey-Greenstein), or 'NC' (Nicolau combination).
        A1 (float): Legacy parameter, use g1 instead (default: 0).
        g1 (float): First anisotropy parameter (default: 0).
        g2 (float): Second anisotropy parameter for NC scattering (default: 0).
        f1 (float): First fraction parameter for NC scattering (default: 0).
        f2 (float): Second fraction parameter for NC scattering (default: 0).
        theta (float, optional): Detection angle in radians.
        theta_in (float, optional): Incident/Source/Emitting angle in radians.
        div_angle (float, optional): Divergence angle for collimated beam.
        machine (str): SPARTA machine configuration (default: "serial").
        cmdargs (list): Command line arguments for SPARTA (default: ["-screen","none"]).
        
    Returns:
        float: Transmission coefficient corrected for the specific angular configuration.
    """
    
    # ensure old version using A1 is working if not replaced with g1
    if (A1 != 0):
        warnings.warn("A1 is deprecated, use g1 instead.", DeprecationWarning)
    g1 = A1 if (A1 != 0) else g1

    if (div_angle and not is_num(theta_in)): raise ValueError("Incident direction is not defined.")
    
    sigma_sca = ext * omega
    kappa = ext * (1-omega)

    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.3f}.trans")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.3f}.trans")

    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"Case already exists and is not empty: {outfile}")
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
        spa.command("surf_collide 1 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")
        # cold black surface
        spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")
        # black surface (unit emission)
        if (div_angle):
            a1 = np.sin(theta_in)
            a2 = -np.cos(theta_in)
            a3 = 0
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25} "+
                        f" collimated {a1} {a2} {a3} div {div_angle}")
        elif (is_num(theta_in)):
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}"+ \
                        f" theta_in {theta_in}")
        else:
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}")

        # hack to delete photon
        spa.command("surf_react 1 global 1.0 0.0")

        # assign collide modules
        spa.command("bound_modify xlo xhi zlo zhi collide 1 react 1")
        spa.command("bound_modify ylo collide 2 react 1")
        spa.command("bound_modify yhi collide 3 react 1")

        if (is_num(theta)):
            spa.command(f"fix 1 emit/face/photon ylo n {nRays} nevery 999999999 limit xzCenter theta {theta}")
        else:
            spa.command(f"fix 1 emit/face/photon ylo n {nRays} nevery 999999999 limit xzCenter")

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


def calc_ref(thickness, ext, omega,
                nRays=200000, SF='LA', A1=0, 
                g1=0, g2=0, f1=0, f2=0,
                theta=None, theta_in=None, div_angle=None,
                machine = "serial", 
                cmdargs = ["-screen","none"]):
    """
    Calculate reflectivity of a medium.
    
    Args:
        thickness (float): Thickness of the medium slab.
        ext (float): Extinction coefficient of the medium.
        omega (float): Single scattering albedo (ratio of scattering to extinction).
        nRays (int): Number of photon rays to simulate (default: 200000).
        SF (str): Scattering function type - 'LA' (Linear Anisotropic), 
                 'HG' (Henyey-Greenstein), or 'NC' (Nicolau combination).
        A1 (float): Legacy parameter, use g1 instead (default: 0).
        g1 (float): First anisotropy parameter (default: 0).
        g2 (float): Second anisotropy parameter for NC scattering (default: 0).
        f1 (float): First fraction parameter for NC scattering (default: 0).
        f2 (float): Second fraction parameter for NC scattering (default: 0).
        theta (float, optional): Detection angle in radians.
        theta_in (float, optional): Incident/Source/Emitting angle in radians.
        div_angle (float, optional): Divergence angle for collimated beam.
        machine (str): SPARTA machine configuration (default: "serial").
        cmdargs (list): Command line arguments for SPARTA (default: ["-screen","none"]).
        
    Returns:
        float: Reflection coefficient corrected for the specific angular configuration.
    """
    
    # ensure old version using A1 is working if not replaced with g1
    if (A1 != 0):
        warnings.warn("A1 is deprecated, use g1 instead.", DeprecationWarning)
    g1 = A1 if (A1 != 0) else g1

    if (div_angle and not is_num(theta_in)): raise ValueError("Incident direction is not defined.")

    sigma_sca = ext * omega
    kappa = ext * (1-omega)

    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.3f}.ref")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.3f}.ref")

    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"Case already exists and is not empty: {outfile}")
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
        spa.command("surf_collide 1 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")
        # cold black surface
        spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")
        # black surface (unit emission)
        if (div_angle):
            a1 = np.sin(theta_in)
            a2 = -np.cos(theta_in)
            a3 = 0
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25} "+
                        f" collimated {a1} {a2} {a3} div {div_angle}")
        elif (is_num(theta_in)):
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}"+ \
                        f" theta_in {theta_in}")
        else:
            spa.command(f"surf_collide 3 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp {(1/SIGMA)**0.25}")

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


def calc_abs(thickness, ext, omega,
                nRays=200000, SF='LA', A1=0, 
                g1=0,
                theta=-1,
                machine = "serial",
                cmdargs = ["-screen","none"]):
    """
    Calculate absorptivity (==emissivity) of a medium.
    
    Args:
        thickness (float): Thickness of the medium slab.
        ext (float): Extinction coefficient of the medium.
        omega (float): Single scattering albedo (ratio of scattering to extinction).
        nRays (int): Number of photon rays to simulate (default: 200000).
        SF (str): Scattering function type - 'LA' (Linear Anisotropic) or 
                 'HG' (Henyey-Greenstein).
        A1 (float): Legacy parameter, use g1 instead (default: 0).
        g1 (float): Anisotropy parameter (default: 0).
        theta (float): Detection angle in radians. (default: -1 for hemispherical).
        machine (str): SPARTA machine configuration (default: "serial").
        cmdargs (list): Command line arguments for SPARTA (default: ["-screen","none"]).
        
    Returns:
        float: Absorption coefficient.
    """
    
    # ensure old version using A1 is working if not replaced with g1
    if (A1 != 0):
        warnings.warn("A1 is deprecated, use g1 instead.", DeprecationWarning)
    g1 = A1 if (A1 != 0) else g1

    sigma_sca = ext * omega
    kappa = ext * (1-omega)

    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.3f}.abs")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.3f}.abs")

    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"Case already exists and is not empty: {outfile}")
    else:
        spa = sparta(machine, cmdargs) 

        spa.command("seed 8887435")
        spa.command("units si")
        spa.command("dimension 3")
        spa.command("boundary ss ss ss")

        if (SF == 'LA'): # Linear Anisotropic
            spa.command("global radiation RMCRT pathlength 1 \
            sigma_sca " + str(sigma_sca) + " g1 " + str(g1) +" kappa " + str(kappa) +  " T " + str((1/SIGMA)**0.25))
        elif (SF == 'HG'):  # Henyey-Greenstein
            spa.command("global radiation RMCRT pathlength 1 \
            sigma_sca " + str(sigma_sca) + " HG g1 "+ str(g1) +" kappa " + str(kappa) +  " T " + str((1/SIGMA)**0.25))
        else:
            sys.exit(SF+" is not recognized.")

        spa.command("photon_continuous none constant")
        spa.command("create_box -1 1 0 "+ str(thickness)+ " -1 1")
        spa.command("create_grid 1 1 1 block * * *")

        spa.command("timestep 0.5")
        spa.command("balance_grid rcb part")

        # cold mirror surface
        spa.command("surf_collide 1 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")
        # cold black surface
        spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")

        # hack to delete photon
        spa.command("surf_react 1 global 1.0 0.0")

        # assign collide modules
        spa.command("bound_modify xlo xhi zlo zhi collide 1 react 1")
        spa.command("bound_modify ylo yhi collide 2 react 1")

        #create_photons n 1000
        if (theta == -1):
            spa.command("fix 1 emit/face/photon yhi n " + str(nRays) + " nevery 999999999 limit xzCenter")
        else:
            spa.command("fix 1 emit/face/photon yhi n " + str(nRays) + " nevery 999999999 limit xzCenter theta " + str(theta))

        spa.command("variable epsilon equal 1")
        spa.command("variable T equal 0")

        spa.command("stats 5000")
        spa.command("run 1")

        spa.command("dump 1 grid all 999999999 " + outfile + " id xc yc zc qrad countEmitted")
        spa.command("run 200000000")

        spa.close()

    return -read_prop(outfile)


def calc_emission(ext, omega, T, limits=[0, 1], \
                size=1 , nRays=200000, SF='LA', g1=0, \
                nonhomogeneous=False, \
                machine = "serial", \
                cmdargs = ["-screen","none"]):
    """
    Calculate radiative heat flux (qrad, W/m2).

    The medium is confined between 2 cold black plates in ylo and yhi directions. 
    The emission (qrad) is calculated at yhi surface.

    Args:
        kappa (float): Absorption coefficient.
        sigma_sca (float): Scattering coefficient.
        T (float or str): Temperature of the medium, can be a constant value or a file.
            file should contain x, T if nonhomogeneous is False,
            or x, T, kappa, sigma_sca if nonhomogeneous is True.
        limits (list): [ylo, yhi] limits of the grid in y direction.
        size (int): Size of the grid along y axis (default is 1).
        nRays (int): Number of rays to be emitted from each cell (default is 200000).
        SF (str): Scattering function, 'LA' for Linear Anisotropic, 'HG' for Henyey-Greenstein.
        g1 (float): Anisotropy factor for LA or HG.
        nonhomogeneous (bool): If True, T_profile is a file with x, T, kappa, and sigma_sca values. (default: False).
        machine (str): SPARTA machine configuration (default: "serial").
        cmdargs (list): Command line arguments for SPARTA (default: ["-screen","none"]).

    Returns:
        float: radiative emission at yhi, W/m2.
    """

    sigma_sca = ext * omega
    kappa = ext * (1-omega)
    D = abs(limits[1] - limits[0])

    outfile = os.path.join(mydir,f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}-size{size}-nRays{nRays}.emi")

    T_profile = T
    if (type(T) is str):
        T = 0

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"Case already exists and is not empty: {outfile}")
    else:

        spa = sparta(machine, cmdargs) 

        spa.command("seed 8887435")
        spa.command("units si")
        spa.command("dimension 3")
        spa.command("boundary ss ss ss")

        if (SF in ['LA', 'HG']):
            spa.command("global radiation RMCRT pathlength 1 \
            sigma_sca " + str(sigma_sca) + " " + SF + " g1 " + str(g1) +" kappa " + str(kappa) +  " T " + str(T))
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

        # cold mirror surface
        spa.command("surf_collide 1 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")
        # cold black surface
        spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")

        # hack to delete photon
        spa.command("surf_react 1 global 1.0 0.0")

        # assign collide modules
        spa.command("bound_modify xlo xhi zlo zhi collide 1 react 1")
        spa.command("bound_modify ylo yhi collide 2 react 1")

        spa.command("fix 1 emit/face/photon yhi n " + str(nRays) + " nevery 999999999 limit xzCenter")

        spa.command("variable epsilon equal 1")
        spa.command("variable T equal 0")

        spa.command("stats 5000")
        spa.command("run 1")

        spa.command("dump 1 grid all 999999999 " + outfile + " id xc yc zc qrad countEmitted")
        spa.command("run 200000000")

        spa.close()

    return -read_prop(outfile)


def is_num(var):
    """
    Check if a variable is a numeric type (int or float).
    
    Utility function to determine if a variable contains a numeric value,
    used throughout the code to validate input parameters.
    
    Args:
        var: Variable to check for numeric type.
        
    Returns:
        bool: True if variable is int or float, False otherwise.
    """
    return isinstance(var, (int, float))


def read_prop(fileName):
    """
    Read a single property value from SPARTA output file.
    
    Parses SPARTA output files to extract a single property value from
    the first data line. Used by calc_trans, calc_ref, and calc_abs functions.
    
    Args:
        fileName (str): Path to the SPARTA output file.
        
    Returns:
        float: The property value (5th column) from the first data line.
    """
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


def calc_EWET(t:np.ndarray, T:np.ndarray, D:float, beta:float, omega:float, SF:str, g1:float):
	"""
	Calculate emission from infinite slab with nonuniform temperature 
	using the exponential weighted effective temperature (EWET) emission model.

	Args:
		t (np.ndarray): Indepth distance from radiating surface.
		T (np.ndarray): Temperature profile.
		D (float): Thickness of the slab.
		beta (float): Absorption coefficient (1/m).
		omega (float): Scattering coefficient (1/m).
		SF (str): Scattering function type, e.g., "HG" for Henyey-Greenstein.
		g1 (float): Asymmetry parameter for scattering.

	Returns:
		float: Emission from the slab (W/m^2).
	"""
	if SF == "HG":
		p = [0.27, -0.96, 1.35, -0.4098]

		rhos = np.exp(-(p[0]*(1-omega)**2+p[1]*omega**2+p[2]+p[3]*g1*omega)*beta*t)
		T4 = np.average(np.power(T, 4), weights=rhos)
		epsilon = calc_abs(D, beta, omega, SF="HG", g1=g1) # to be replaced with emissivity model
		return epsilon * SIGMA * T4
	else:
		raise NotImplementedError(f"SF={SF} not implemented")
	

def calc_ED(q:float, t:np.ndarray, D:float, beta:float, omega:float, SF:str, g1:float):
	"""
	Calculate exponential decay of div.q within an infinite slab using the 
	exponential decay (ED) model.

	Args:
		q (float): Incident heat flux.
		t (np.ndarray): Thickness of the slab.
		D (float): Thickness of the slab.
		beta (float): Absorption coefficient (1/m).
		omega (float): Scattering coefficient (1/m).
		SF (str): Scattering function type, e.g., "HG" for Henyey-Greenstein.
		g1 (float): Asymmetry parameter for scattering.

	Returns:
		float: Divergence of q profile along the slab (W/m^3). 
		Negated so that positive is for gain, negative is for loss.
	"""
	if SF == "HG":
		# optimization with wieghts (favors indepth fitting)
		c0, c1, c2, p1 = [ 1.2850711,  -0.48658091, -0.51118171,  4.87108912]
		c11, c21, c3 = [-0.70981822,  0.46435636,  3.10418282]
		# if beta*(1-omega) < 500:
		# 	# without using weights (favors near surface fitting)
		# 	c0, c1, c2, p1 = [ 1.74447467, -0.9899896,  -0.46715914,  4.30635443]
		# 	c11, c21, c3 = [-0.7545377,   0.57699982,  1.38960422]


		z = (c0 + (c1+c11*g1)*omega + (c2+c21*g1)*omega**(p1+c3*g1))*beta
		rho = calc_ref(D, beta, omega, SF="HG", g1=g1) # to be replaced with reflectivity model
		dq = z * np.exp(-z*t) * (1-rho) * q
		return dq
	else:
		raise NotImplementedError(f"SF={SF} not implemented")
	