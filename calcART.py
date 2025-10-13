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
			nRays=1000, SF='LA', g1=0, \
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
        raise ValueError("Error: Thickness is too small, results may conflict.")
    if not outputName:
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
	Calculate exponential decay of div.q within an infinite slab using the 
	exponential decay (ED) model.

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
				[ 1.77884405, -1.13756974, -0.36217746,  5.55132145], # tau_t <= 1
				[ 1.57960717, -0.6777506 , -0.74774612,  4.79422544], # tau_t <= 2
				[ 1.44528636, -0.54006825, -0.55390418,  3.52680565], # tau_t <= 3
				[ 1.36353177, -0.46856152, -0.56847874,  3.6571423 ], # tau_t <= 4
				[ 1.3200944 , -0.42218924, -0.5772813 ,  3.74555988], # tau_t <= 5
				[ 1.28631732, -0.39120893, -0.57786769,  3.78667048], # tau_t <= 6
				[ 1.25772036, -0.35609818, -0.59024459,  3.84002667], # tau_t <= 7
				[ 1.24058056, -0.33290104, -0.5907199 ,  3.78264772], # tau_t <= 8
				[ 1.22113859, -0.30834976, -0.60007301,  3.82209053], # tau_t <= 9
				[ 1.20685848, -0.30610587, -0.59529843,  3.970225  ], # tau_t <= 10
				[ 1.19614703, -0.27885995, -0.59757553,  3.79007963], # tau_t <= 11
				[ 1.18550901, -0.26553828, -0.61629813,  3.95148033], # tau_t <= 12
				[ 1.17369848, -0.23615239, -0.62423467,  3.76668054], # tau_t <= 13
				[ 1.16712712, -0.24975551, -0.62754506,  4.1340121 ], # tau_t <= 14
				[ 1.16077959, -0.23867797, -0.63290292,  4.15252798], # tau_t <= 15
				[ 1.15602513, -0.2439741 , -0.6402916 ,  4.38802029], # tau_t <= 16
				[ 1.14394247, -0.16354218, -0.69043576,  3.84412933], # tau_t <= 17
				[ 1.13989116, -0.15202603, -0.70568381,  3.89785306], # tau_t <= 18
				[ 1.13516762, -0.12676313, -0.74381222,  3.9069203 ], # tau_t <= 19
				[ 1.12791251, -0.10562627, -0.79108487,  4.22406844], # tau_t <= 20
                    ]
		ps = np.array(ps)
		for depth in t:
			tau_t = depth * beta
			if tau_t <= 1:
				c.append(np.concatenate([ps[0], [0, 0, 0]]))
			elif tau_t <= 2:
				c.append(np.concatenate([ps[1], [0, 0, 0]]))
			elif tau_t <= 3:
				c.append(np.concatenate([ps[2], [0, 0, 0]]))
			elif tau_t <= 4:
				c.append(np.concatenate([ps[3], [0, 0, 0]]))
			elif tau_t <= 5:
				c.append(np.concatenate([ps[4], [0, 0, 0]]))
			elif tau_t <= 6:
				c.append(np.concatenate([ps[5], [0, 0, 0]]))
			elif tau_t <= 7:
				c.append(np.concatenate([ps[6], [0, 0, 0]]))
			elif tau_t <= 8:
				c.append(np.concatenate([ps[7], [0, 0, 0]]))
			elif tau_t <= 9:
				c.append(np.concatenate([ps[8], [0, 0, 0]]))
			elif tau_t <= 10:
				c.append(np.concatenate([ps[9], [0, 0, 0]]))
			elif tau_t <= 11:
				c.append(np.concatenate([ps[10], [0, 0, 0]]))
			elif tau_t <= 12:
				c.append(np.concatenate([ps[11], [0, 0, 0]]))
			elif tau_t <= 13:
				c.append(np.concatenate([ps[12], [0, 0, 0]]))
			elif tau_t <= 14:
				c.append(np.concatenate([ps[13], [0, 0, 0]]))
			elif tau_t <= 15:
				c.append(np.concatenate([ps[14], [0, 0, 0]]))
			elif tau_t <= 16:	
				c.append(np.concatenate([ps[15], [0, 0, 0]]))
			elif tau_t <= 17:
				c.append(np.concatenate([ps[16], [0, 0, 0]]))
			elif tau_t <= 18:
				c.append(np.concatenate([ps[17], [0, 0, 0]]))
			elif tau_t <= 19:
				c.append(np.concatenate([ps[18], [0, 0, 0]]))
			elif tau_t <= 20:
				c.append(np.concatenate([ps[19], [0, 0, 0]]))
			else: # tau_t > 20
				# print(f"Warning: tau_t={tau_t:0.1f} > 20 at depth of {depth*1000:0.3f}mm, using coefficients for tau [19:20]")
				c.append(np.concatenate([ps[19], [0, 0, 0]]))
			
               
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
	