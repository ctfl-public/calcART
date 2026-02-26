from sparta import sparta 
import os
from constants import SIGMA
import sys
import numpy as np
import warnings
from data_management import check_version_file

# initialize data management
mydir = check_version_file()
print(f"Data directory: {mydir}")


def calc_dq(kappa:float, sigma_sca:float, T:str|float, outputName=None, limits=[0, 1], \
            size=1 , nRays=1000, SF='LA', g1=0, \
            nonhomogeneous=False, \
            in_rad=0.0, \
            machine = "serial", \
            cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux (Dqrad, W/m3) of a medium slab 
    confined between: 
    - ylo: black surfaces.
    - yhi: black surface with radiation input (in_rad). 

    Temperatures of ylo and yhi black surfaces are T, if T is float, or T[0], if T is MR_profile path.

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
        machine (str): Name of the machine used to run sparta, e.g., "serial", "serial_debug".
        cmdargs (list): Additional command line arguments for sparta.

    Returns:
        tuple: A tuple containing two lists:
            - yc (list): Y-coordinate values
            - dqrad (list): Radiative heat flux values (negated from file)
                positive values indicate incoming radiation (heat gain),
                negative values indicate outgoing radiation (heat loss).
    """
    
    D = abs(limits[1] - limits[0])
    if (D*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
        outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}mm-size{size:0.0f}-nRays{nRays:0.0f}-inrad{in_rad:0.1e}.dq"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        # print(f"Case already exists and is not empty: {outfile}")
        return read_dqrad(outputName)

    if isinstance(T, str):
        MR_profile = T
        T = 0 
        # read MRProfile and assign Tlo and Thi
        xMR, T_xMR = readMRProfile(MR_profile)
        Tlo = np.interp(limits[0], xMR, T_xMR)
        Thi = np.interp(limits[1], xMR, T_xMR)
    else:
        MR_profile = None
        Tlo = T
        Thi = T

    if (in_rad):
        Thi = (in_rad/SIGMA)**0.25 

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

    if (type(MR_profile) is str):
        # variable kappa and sigma, file contains x, T, kappa, and sigma
        if nonhomogeneous: 
            spa.command("global MR " + MR_profile)
        # constant kappa and sigma, file contains only x and T
        else: 
            spa.command("global T MR " + MR_profile)

    spa.command("timestep 0.5")
    spa.command("balance_grid rcb part")

    # black surface
    spa.command("surf_collide 1 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str(Thi))
    # cold black surface
    spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str(Tlo))
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

    return read_dqrad(outputName)


def calc_dq_medium(kappa:float, sigma_sca:float, T:str|float, outputName=None, limits=[0, 1], \
            size=1 , nRays=1000, SF='LA', g1=0, \
            nonhomogeneous=False, \
            machine = "serial", \
            cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux (Dqrad, W/m3) from neighboring cells
    of a medium slab confined between 2 cold black plates.

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
        machine (str): Name of the machine used to run sparta, e.g., "serial", "serial_debug".
        cmdargs (list): Additional command line arguments for sparta.

    Returns:
		list: divergence of radiative heat flux values
			positive values indicate incoming radiation (heat gain),
			negative values indicate outgoing radiation (heat loss).
    """
    
    D = abs(limits[1] - limits[0])
    
    if (D*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
        outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}mm-size{size:0.0f}-nRays{nRays:0.0f}.dq_medium"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        # print(f"Case already exists and is not empty: {outfile}")
        _, I = read_dqrad(outputName)
        dq_medium = -4*np.pi*kappa*np.array(I)  
        return dq_medium

    if isinstance(T, str):
        MR_profile = T
        T = 0 
        # # read MRProfile and assign Tlo and Thi
        # xMR, T_xMR = readMRProfile(MR_profile)
        # Tlo = np.interp(limits[0], xMR, T_xMR)
        # Thi = np.interp(limits[1], xMR, T_xMR)
    else:
        MR_profile = None
        # Tlo = T
        # Thi = T
        
    Tlo = 0
    Thi = 0

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

    if (type(MR_profile) is str):
        # variable kappa and sigma, file contains x, T, kappa, and sigma
        if nonhomogeneous: 
            spa.command("global MR " + MR_profile)
        # constant kappa and sigma, file contains only x and T
        else: 
            spa.command("global T MR " + MR_profile)

    spa.command("timestep 0.5")
    spa.command("balance_grid rcb part")

    # black surface
    spa.command("surf_collide 1 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str(Thi))
    # black surface
    spa.command("surf_collide 2 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str(Tlo))
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
    spa.command("dump 1 grid all 999999999 " + outfile + " yc I")
    spa.command("run 200000000")

    spa.close()

    _, I = read_dqrad(outputName)
    dq_medium = -4*np.pi*kappa*np.array(I)  
    return dq_medium


def calc_dq_cooling(kappa:float, sigma_sca:float, T:float, D:float, \
			nRays=200000, SF='LA', g1=0, \
            machine = "serial", \
            cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux of infinite slab of thickness D at temperature T surrounded by cold black medium.
    
    Args:
		kappa (float): Absorption coefficient.
		sigma_sca (float): Scattering coefficient.
		T (float): Temperature of the medium. (constant)
		D (float): Thickness of the slab.
		nRays (int): Number of rays to be emitted from each cell (default is 1000).
		SF (str): Scattering function, 'LA' for Linear Anisotropic, 'HG' for Henyey-Greenstein.
		g1 (float): Anisotropy factor for LA or HG.
          
    Returns:
		float: Dqrad (W/m3) value. Positive values indicate heat loss.
    """

    
    if (D*1e6 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    outputName = f"T{T:0.4f}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1e6:0.3f}microns-nRays{nRays:0.0f}.dq_cooling"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        # print(f"Case already exists and is not empty: {outfile}")
        return read_prop(outfile)

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
    spa.command("create_box -1 1 0 " + int2str(D) + " -1 1")
    spa.command("create_grid 1 1 1 block * * *")

    spa.command("timestep 0.5")
    spa.command("balance_grid rcb part")

    # old black surface
    spa.command("surf_collide 1 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp 0")
    # cold mirror surface
    spa.command("surf_collide 2 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")

    # hack to delete photon
    spa.command("surf_react 1 global 1.0 0.0")

    # assign collide modules
    spa.command("bound_modify ylo yhi collide 1 react 1")
    spa.command("bound_modify xlo xhi zlo zhi collide 2 react 1")

    spa.command("create_photons npc " + str(nRays) + " xzCenter")

    spa.command("stats 5000")
    spa.command("run 1")
    spa.command("dump 1 grid all 999999999 " + outfile + " id xc yc zc Dqrad countEmitted")
    spa.command("run 200000000")

    spa.close()

    return read_prop(outfile)


def readMRProfile(file):
		with open(file, 'r') as f:
			lines = f.readlines()
		# Initialize lists for storing the data
		x_values = []
		y_values = []
		# Process the lines and extract the data from the first and second columns
		for line in lines:
			if not (line.startswith(' #') or line.startswith('#')):
				columns = line.split()
				x_values.append(float(columns[0]))
				y_values.append(float(columns[1]))
		xMR = np.array(x_values)
		T_xMR = np.array(y_values)

		return xMR, T_xMR


def calc_dqrad(kappa, sigma_sca, T:str|float, outputName=None, limits=[0, 1], \
                size=1 , nRays=1000, SF='LA', g1=0, \
                nonhomogeneous=False, \
                in_rad=0, \
                machine = "serial", \
                cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux (Dqrad, W/m3) of a medium slab 
    confined between: 
    - ylo: cold black surfaces
    - yhi: black surface with radiation input (in_rad). If in_rad = 0, 
    the top surface is also a cold black surface.

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
    if (D <= 0.001):
        high_precision = True
    elif (D*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
        if high_precision:
            outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}mm-size{size:0.0f}-nRays{nRays:0.0f}-inrad{in_rad:0.1e}.dqrad"
        else:
            outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D:0.3f}-size{size:0.0f}-nRays{nRays:0.0f}-inrad{in_rad:0.1e}.dqrad"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        #print(f"Case already exists and is not empty: {outfile}")
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


def calc_dq_equilibrium(kappa:float, sigma_sca:float, T:float, outputName=None, limits=[0, 1], \
                size=1 , nRays=1000, SF='LA', g1=0, \
                machine = "serial", \
                cmdargs = ["-screen","none"]):
    """
    Calculate divergence of radiative heat flux (Dqrad, W/m3). 
    The slab medium is isothermal and in equilibrium with the surrounding.

    The medium is confined between ylo and yhi in y direction;
    ylo is cold black surface, yhi is black surface with {in_rad} radiation input.

    Args:
        kappa (float): Absorption coefficient.
        sigma_sca (float): Scattering coefficient.
        T (float): Temperature of the medium. (constant)
        outputName (str): Name of output file (includes extension). 
            if None, output is renamed to default name. (default is None).
        limits (list): [ylo, yhi] limits of the grid in y direction.
        size (int): Size of the grid along y axis (default is 1).
        nRays (int): Number of rays to be emitted from each cell (default is 1000).
        SF (str): Scattering function, 'LA' for Linear Anisotropic, 'HG' for Henyey-Greenstein.
        g1 (float): Anisotropy factor for LA or HG.

    Returns:
        tuple: A tuple containing two lists:
            - yc (list): Y-coordinate values
            - dqrad (list): Radiative heat flux values (negated from file)
                positive values indicate incoming radiation (heat gain),
                negative values indicate outgoing radiation (heat loss).
    """
    
    D = abs(limits[1] - limits[0])
    if (D <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
        outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D:0.3f}-size{size}-nRays{nRays}.dqequilibrium"
    outfile = os.path.join(mydir,outputName)

    # skip runing of file exists
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        #print(f"Case already exists and is not empty: {outfile}")
        return read_dqrad(outputName)

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

    spa.command("timestep 0.5")
    spa.command("balance_grid rcb part")

    # black surface
    spa.command("surf_collide 1 radiationboundary epsilon 1 rho_s 0 rho_d 0 temp "+ str(T))
    # cold mirror surface
    spa.command("surf_collide 2 radiationboundary epsilon 0 rho_s 1 rho_d 0 temp 0")

    # hack to delete photon
    spa.command("surf_react 1 global 1.0 0.0")

    # assign collide modules
    spa.command("bound_modify ylo yhi collide 1 react 1")
    spa.command("bound_modify xlo xhi zlo zhi collide 2 react 1")

    spa.command("create_photons npc " + str(nRays) + " xzCenter")

    spa.command("stats 5000")
    spa.command("run 1")
    spa.command("dump 1 grid all 999999999 " + outfile + " yc Dqrad")
    spa.command("run 200000000")

    spa.close()

    return read_dqrad(outputName)


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
    
    if (thickness*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.3f}.trans")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.3f}.trans")

    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
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

    if (thickness*1e6 <= 0.001):
        raise ValueError(f"Error: Thickness {thickness*1e6:.3f} microns is too small, results may conflict.")
    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1e6:.3f}microns.ref")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1e6:.3f}microns.ref")

    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
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


def calc_emissivity_bezier(omega, SF='HG', g1=0):
    """
    Calculate emissivity using Bezier curve fit based on ext, omega, and scattering function.
    
    Args:
        omega (float or array-like): Single scattering albedo values. Must be within [0, 1].
        SF (str): Scattering function type - 'LA' (Linear Anisotropic) or 'HG' (Henyey-Greenstein).
        g1 (float): Anisotropy parameter (default: 0).
        
    Returns:
        float or numpy.ndarray: Emissivity value calculated from the Bezier curve fit with the same shape as `omega`.
    """

    if (SF == 'HG'):
        return bezier_HG(omega, g1)
    else:
        raise ValueError(SF+" is not recognized.")


def bezier_HG(omega, g1, n_t=101):
    """
    Approximate emissivity using the cubic Bezier HG model from tabulated control points.

    Parameters
    ----------
    omega : float or array-like
        Single scattering albedo values. Must be within [0, 1].
    g1 : float
        Asymmetry factor. Must be within [-0.8, 0.8].
    n_t : int, optional
        Number of parametric samples used to invert omega(t) -> epsilon(t).

    Returns
    -------
    float or numpy.ndarray
        Emissivity values with the same shape as `omega`.
    """
    # HG control points from the paper's Table 1.
    g_table = np.array([-0.8, -0.6, -0.3, 0.0, 0.3, 0.6, 0.8], dtype=float)
    p1x_table = np.array([0.856, 0.819, 0.846, 0.849, 0.859, 0.894, 0.942], dtype=float)
    p1y_table = np.array([0.648, 0.706, 0.764, 0.825, 0.880, 0.930, 0.948], dtype=float)
    p2x_table = np.array([0.997, 0.998, 0.998, 0.998, 0.998, 0.997, 0.995], dtype=float)
    p2y_table = np.array([0.412, 0.465, 0.483, 0.543, 0.636, 0.778, 0.954], dtype=float)

    # Validate g1 range
    if g1 < g_table[0] or g1 > g_table[-1]:
        raise ValueError("g1 must be within [-0.8, 0.8] for the HG Bezier model.")

    # Interpolate control points for the given g1.
    p1x = np.interp(g1, g_table, p1x_table)
    p1y = np.interp(g1, g_table, p1y_table)
    p2x = np.interp(g1, g_table, p2x_table)
    p2y = np.interp(g1, g_table, p2y_table)

    # Validate omega range and convert to numpy array for processing.
    omega_arr = np.asarray(omega, dtype=float)
    omega_shape = omega_arr.shape
    omega_flat = np.atleast_1d(omega_arr).ravel()
    if np.any((omega_flat < 0.0) | (omega_flat > 1.0)):
        raise ValueError("omega must be within [0, 1] for the HG Bezier model.")

    # Cubic Bezier with P0=(0,1), P3=(1,0), and interpolated P1/P2.
    # omega(t) = bx(t) = Bezier curve in x direction
    # epsilon(t) = by(t) = Bezier curve in y direction.
    t = np.linspace(0.0, 1.0, int(n_t), dtype=float)
    one_minus_t = 1.0 - t
    bx = (
        3.0 * (one_minus_t ** 2) * t * p1x
        + 3.0 * one_minus_t * (t ** 2) * p2x
        + (t ** 3)
    )
    by = (
        (one_minus_t ** 3)
        + 3.0 * (one_minus_t ** 2) * t * p1y
        + 3.0 * one_minus_t * (t ** 2) * p2y
    )

    # Invert the Bezier curve to find epsilon for each omega.
    # (instead of both being functions of the curve parameter t)
    sort_idx = np.argsort(bx)
    bx_sorted = bx[sort_idx]
    by_sorted = by[sort_idx]
    bx_unique, unique_idx = np.unique(bx_sorted, return_index=True)
    by_unique = by_sorted[unique_idx]
    if bx_unique.size < 2:
        raise ValueError("Bezier curve inversion failed: omega(t) is not usable for interpolation.")

    # Interpolate to find epsilon values corresponding to the input omega values.
    eps_flat = np.interp(omega_flat, bx_unique, by_unique)
    # Reshape the output to match the input shape of omega.
    eps = eps_flat.reshape(omega_shape)

    if omega_shape == ():
        return float(eps_flat[0])
    return eps


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

    if (thickness*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if (SF == 'LA'):
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-th{thickness*1000:.3f}.abs")
    else:
        outfile = os.path.join(mydir,f"ext{ext:.0f}-sigma{sigma_sca:.0f}-HG-g{g1:0.3f}-th{thickness*1000:.3f}.abs")

    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
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
                outputName=None, \
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

    if (D*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
        outputName = f"T{T}-abs{kappa:0.0f}-sca{sigma_sca:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}-size{size}-nRays{nRays}.emi"
    outfile = os.path.join(mydir, outputName)

    T_profile = T
    if (type(T) is str):
        T = 0

    # skip runing if file exists
    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
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


def calc_EWET(t:np.ndarray, T:np.ndarray, beta:float, omega:float, SF:str, g1:float, D:float=None, eps:float=None):
	"""
	Calculate emission from infinite slab with nonuniform temperature 
	using the exponential weighted effective temperature (EWET) emission model.

	Args:
		t (np.ndarray): Indepth distance from radiating surface.
		T (np.ndarray): Temperature profile.
		beta (float): Absorption coefficient (1/m).
		omega (float): Scattering coefficient (1/m).
		SF (str): Scattering function type, e.g., "HG" for Henyey-Greenstein.
		g1 (float): Asymmetry parameter for scattering.
		D (float, None if eps is provided): Thickness of the slab (m).
        eps (float, None if D is provided): Emissivity of the slab. If None, it will be calculated using calc_abs.

	Returns:
		float: Emission from the slab (W/m^2).
	"""
	if D is None and eps is None:
		raise ValueError("At least one of 'D' or 'eps' must be provided.")

	if SF == "HG":
		p = [0.27, -0.96, 1.35, -0.4098]

		rhos = np.exp(-(p[0]*(1-omega)**2+p[1]*omega**2+p[2]+p[3]*g1*omega)*beta*t)
		T4 = np.average(np.power(T, 4), weights=rhos)
		if eps is None:
			eps = calc_abs(D, beta, omega, SF="HG", g1=g1) # to be replaced with emissivity model
		# print(f"\t\t\tEWET: eps={eps}, T={T4**0.25:.2f} K")
		return eps * SIGMA * T4
	else:
		raise NotImplementedError(f"SF={SF} not implemented")
	

def calc_ED(q:float, t:np.ndarray, beta:float, omega:float, SF:str, g1:float, D:float=None, rho:float=None):
	"""
	The exponential decay (ED) model of div.q within an infinite slab, independent of dx.

	Args:
		q (float): Incident heat flux.
		t (np.ndarray): indepth distance from radiating surface.
		beta (float): Absorption coefficient (1/m).
		omega (float): Scattering coefficient (1/m).
		SF (str): Scattering function type, e.g., "HG" for Henyey-Greenstein.
		g1 (float): Asymmetry parameter for scattering.
		D (float, None if rho is provided): Thickness of the slab (m).
		rho (float, None if D is provided): Reflectivity of the slab. If None, it will be calculated using calc_ref.

	Returns:
		ndarray (float): Divergence of q profile along the slab (W/m^3). 
		Negated so that positive is for gain, negative is for loss.
	"""
	if rho is None and D is None:
		raise ValueError("At least one of 'rho' or 'D' must be provided.")

	c = []
	if isinstance(t, (int, float)):
		t = np.array([t])
	if SF == "HG":
		ps = [
               [1.7491399843309074, -1.1552594930057578, -0.4318258386332926, 9.25024515431138, -0.5761350422793642, 0.4998417390151376, 4.150870636913617],    # tau_t < 1
               [1.578784151010635, -0.6761764004143443, -0.7219166423129564, 4.665250384559308, -0.8561787105800065, 0.7139544714826039, 2.134979292479359],    # tau_t < 2
               [1.4415746726535175, -0.5522202698008893, -0.594943039246587, 3.9763416567871137, -0.6642272289495321, 0.3342743371826334, 1.804684535855158],    # tau_t < 3
               [1.3651711029736548, -0.481249796441904, -0.6037599935279878, 4.096469833658471, -0.6748722255372567, 0.37406994955475137, 2.0727482949134086],   # tau_t < 4
               [1.3229258185051147, -0.45903291820330927, -0.6027260159585395, 4.364585823775098, -0.6645554467694439, 0.35258738146975727, 2.4180688229303366], # tau_t < 5
               [1.289489168059346, -0.4177240924094028, -0.6171045241663791, 4.401463729819236, -0.6808194394661212, 0.33953672020492376, 2.659895035734948], 
               [1.2639583479263201, -0.39451456701717014, -0.6224591148816251, 4.541433669113842, -0.6770165843928688, 0.32161190684697943, 2.8880162869284596], 
               [1.2443726135814845, -0.37595462417131176, -0.6244164406520526, 4.607958830137856, -0.6752681604287076, 0.31279419547367326, 3.0159062794384113], 
               [1.2266970852416734, -0.35750342716665073, -0.6265037133321759, 4.642930732241344, -0.674616139632077, 0.3167836876971489, 3.0404604999464766], 
               [1.213979626639181, -0.33996247153073766, -0.6314472427445118, 4.632819491431386, -0.6789104912870932, 0.32185876117673495, 3.077489453763049], 
               [1.201389330322055, -0.31900168795793266, -0.6406970951116406, 4.6055760296443715, -0.6847230613741695, 0.3392725614748317, 3.03601277803], 
               [1.1906053878924654, -0.312414364262092, -0.6367143776176585, 4.66148011297818, -0.6780910425301713, 0.33444358421543346, 3.043963011306982], 
               [1.1808922764719154, -0.2910265929340613, -0.6498175351272703, 4.626470417109754, -0.69056586099523, 0.3591622217014206, 2.9880700867543575], 
               [1.172712560913378, -0.27924156079224843, -0.6558555187835624, 4.6531383269181, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
               [1.165035577505646, -0.26540870315383525, -0.6632878955929292, 4.674213986337265, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
               [1.157858423330447, -0.24727043908854712, -0.6773358208078232, 4.642706234882837, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
               [1.1541810656604963, -0.23139644393823666, -0.6897639660440309, 4.586850973503088, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
               [1.1449179545338313, -0.19695353455614692, -0.707709827590225, 4.434736751279188, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
               [1.1383146341996437, -0.16914693047171248, -0.7235731819347645, 4.3664519666863555, -0.6980403604238626, 0.37995173344074995, 2.943079450625958], 
            #    [1.1210566165075444, -0.04915201316161835, -0.8223186216189828, 4.264249847720826, -0.6980403604238626, 0.37995173344074995, 2.943079450625958]
               ]

		ps = np.array(ps)
		for depth in t:
			tau_t = depth * beta
			if tau_t <= 1:
				c.append(ps[0])
			elif tau_t <= 2:
				c.append(ps[1])
			elif tau_t <= 3:
				c.append(ps[2])
			elif tau_t <= 4:
				c.append(ps[3])
			elif tau_t <= 5:
				c.append(ps[4])
			elif tau_t <= 6:
				c.append(ps[5])
			elif tau_t <= 7:
				c.append(ps[6])
			elif tau_t <= 8:
				c.append(ps[7])
			elif tau_t <= 9:
				c.append(ps[8])
			elif tau_t <= 10:
				c.append(ps[9])
			elif tau_t <= 11:
				c.append(ps[10])
			elif tau_t <= 12:
				c.append(ps[11])
			elif tau_t <= 13:
				c.append(ps[12])
			elif tau_t <= 14:
				c.append(ps[13])
			elif tau_t <= 15:
				c.append(ps[14])
			# elif tau_t <= 16:	
			# 	c.append(ps[15])
			# elif tau_t <= 17:
			# 	c.append(ps[16])
			# elif tau_t <= 18:
			# 	c.append(ps[17])
			# elif tau_t <= 19:
			# 	c.append(ps[18])
			# elif tau_t <= 20:
			# 	c.append(ps[19])
			else: # tau_t > 20
				# print(f"Warning: tau_t={tau_t:0.1f} > 20 at depth of {depth*1000:0.3f}mm, using coefficients for tau [19:20]")
				c.append(ps[14])
			
               
		c = np.array(c)
		c0 = c[:, 0]
		c1 = c[:, 1]
		c2 = c[:, 2]
		p1 = c[:, 3]
		c11 = c[:, 4]
		c21 = c[:, 5]
		c3 = c[:, 6]
		z = (c0 + (c1+c11*g1)*omega + (c2+c21*g1)*omega**(p1+c3*g1))*beta
		if rho is None:
			rho = calc_ref(D, beta, omega, SF="HG", g1=g1) # to be replaced with reflectivity model
		dq = z * np.exp(-z*t) * (1-rho) * q
		return dq
	else:
		raise NotImplementedError(f"SF={SF} not implemented")
	

def calc_dq_ED(q:float, t:np.ndarray, beta:float, omega:float, SF:str, g1:float, rho:float=None):
	"""
	Calculate exponential decay of div.q within an infinite slab using the 
	exponential decay (ED) model, dependent of dx.

	Args:
		q (float): Incident heat flux.
		t (np.ndarray): indepth distance from radiating surface.
		beta (float): Absorption coefficient (1/m).
		omega (float): Scattering coefficient (1/m).
		SF (str): Scattering function type, e.g., "HG" for Henyey-Greenstein.
		g1 (float): Asymmetry parameter for scattering.
		rho (float, optional): Reflectivity of the slab. If None, it will be calculated using calc_ref.

	Returns:
		ndarray (float): Divergence of q profile along the slab (W/m^3). 
		Negated so that positive is for gain, negative is for loss.
	"""
	if isinstance(t, (int, float)):
		t = np.array([t])

	size_l = len(t)
	dx_l = t[1]-t[0] if size_l > 1 else 2*t[0] # +/- spacing in t
	D = max(t) + abs(dx_l)/2
	if rho is None:
		rho = calc_ref(D, beta, omega, SF=SF, g1=g1) # to be replaced with reflectivity model
     
	model_size = 600
	model_tau = 30
	size_h = int(D*beta * model_size / model_tau)
	dq_interpolated = np.zeros_like(t)
	N = int(size_h/size_l) # number of high-res points per low-res point
	size_h = N * size_l
	if size_l >= size_h:
		# print(f"Warning: size_l={size_l} >= size_h={size_h}, Using ED model directly.")
		return calc_ED(q, t, beta, omega, SF, g1, rho=rho)

	t_0 = t[0] - dx_l/2
	for i in range(size_l):
		t_N = t_0 + dx_l
		t_model = np.linspace(t_0, t_N, N, endpoint=True)
		dq_model = calc_ED(q, t_model, beta, omega, SF, g1, rho=rho)
		dq_interpolated[i] = np.trapz(dq_model, t_model) / dx_l
		t_0 = t_N
	return dq_interpolated


def calc_first_cell_dq(kappa:float, sigma:float, nRays:int, SF:str, g1:float, q:float, D:float, size:int):
    """ 
    Wrapper for calc_dq() to calculate dq for a given cell due to neighbors emission only at the yhi side.
    ylo is treated as cold black.
    
    Args:
        q (float): incident radiation
        D (float): slab thickness
        size (int): number of cells
    """
    print(f"Calculating first_cell_dq using RMCRT for D={D*1000:0.3f} mm, size={size}...", end=' ')
    if (D*1000 <= 0.001):
        raise ValueError("Error: Thickness is too small, results may conflict.")
    outputName = f"T0-abs{kappa:0.0f}-sca{sigma:0.0f}-{SF}-g1{g1:0.3f}-D{D*1000:0.3f}mm-size{size:0.0f}-nRays{nRays:0.0f}-inrad{q:0.1e}.dq_ED_medium"
    _, dq = calc_dq(
        kappa=kappa,
        sigma_sca=sigma,
        T=0.0,
        limits=[0, D],
        size=size,
        nRays=nRays,
        SF=SF,
        g1=g1,
        in_rad=q,
        outputName=outputName
    )
    print("\tDone.")
    return dq[-1]