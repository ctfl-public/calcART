import numpy as np
from calcART import *

class base:
	def __init__(self, beta:float, omega:float, SF:str, g1:float, D:float):
		self.beta = beta
		self.omega = omega
		self.SF = SF
		self.g1 = g1
		self.D = D

		self.sigma = omega * beta
		self.kappa = beta * (1 - omega)

		self._abs = None
		self._ref = None
		self._trans = None
		
	@property
	def abs(self):
		if self._abs is None:
			self._abs = calc_abs(thickness=self.D,
								ext=self.beta,
								omega=self.omega,
								SF=self.SF,
								g1=self.g1)
		return self._abs

	@property
	def ref(self):
		if self._ref is None:
			self._ref = calc_ref(thickness=self.D,
								ext=self.beta,
								omega=self.omega,
								SF=self.SF,
								g1=self.g1)
		return self._ref

	@property
	def trans(self):
		if self._trans is None:
			self._trans = calc_trans(thickness=self.D,
									ext=self.beta,
									omega=self.omega,
									SF=self.SF,
									g1=self.g1)
		return self._trans


class art(base):
	def __init__(self, beta:float, omega:float, SF:str, g1:float, D:float, 
			  T:float|str, size:int, nRays:int, in_rad:float):
		super().__init__(beta=beta, omega=omega, SF=SF, g1=g1, D=D)

		self.size = size
		self.dx = self.D / self.size
		self.x = np.linspace(self.dx/2, self.D-self.dx/2, self.size)
		self.xnorm = self.x / self.D
		self.t = self.D - self.x # distance from the radiating surface
		self.tau_t = self.t * self.beta

		if isinstance(T, str):
			self.MR_profile = T
			self.T = None
			self.xMR, self.T_xMR = readMRProfile(self.MR_profile)
			self.T_x = np.interp(self.x, self.xMR, self.T_xMR)
		else:
			self.MR_profile = None
			self.T = T
			self.xMR = self.x
			self.T_xMR = np.ones(self.size) * T
			self.T_x = np.ones(self.size) * T

		self.nRays = nRays
		self.in_rad = in_rad # applied at xhi

		# All are func of self.x (not xMR), corresponds to self.size
		self._dq_model = None
		self._dq_cooling_term = None
		self._dq_gas_term = None
		self._dq_medium_term = None
		self._dq_RMCRT = None 
		self._dq_cooling_RMCRT = None
		self._dq_gas_RMCRT = None
		self._dq_medium_RMCRT = None


	@property
	def dq_RMCRT(self):
		if self._dq_RMCRT is None:
			print("Calculating dq_RMCRT...")
			print(f"\text={self.beta}, omega={self.omega}, SF={self.SF}, g1={self.g1}, D={self.D}\n" + \
		 			f"\tT={self.T if self.T is not None else self.MR_profile}, in_rad={self.in_rad}\n" + \
					f"\tsize={self.size}, nRays={self.nRays}")	
			_, self._dq_RMCRT = calc_dq(kappa=self.kappa,
											sigma_sca=self.sigma,
											T=self.MR_profile if self.MR_profile else self.T,
											limits=[0,self.D],
											size=self.size,
											nRays=self.nRays,
											SF=self.SF,
											g1=self.g1,
											nonhomogeneous=False,
											in_rad=self.in_rad)
			# multiply by -1 to reverse negation enforced in read_dqrad() caled by calc_dq()
			self._dq_RMCRT = -np.array(self._dq_RMCRT)
		return self._dq_RMCRT

	@property
	def dq_cooling_term(self):
		if self._dq_cooling_term is None:
			eps = calc_abs(thickness=self.dx,
								ext=self.beta,
								omega=self.omega,
								SF=self.SF,
								g1=self.g1)
			self._dq_cooling_term = np.array(2*eps*SIGMA*np.power(self.T_x, 4)/self.dx)
		return self._dq_cooling_term
	

	@property
	def dq_cooling_RMCRT(self):
		if self._dq_cooling_RMCRT is None:
			dq_cooling = []
			for i in range(self.size):
				dq_cooling.append(
					calc_dq_cooling(kappa=self.kappa,
									sigma_sca=self.sigma,
									T=self.T_x[i],
									D=self.dx,
									SF=self.SF,
									g1=self.g1,
									nRays=self.nRays)
				)
			self._dq_cooling_RMCRT = np.array(dq_cooling)
		return self._dq_cooling_RMCRT
	

	@property
	def dq_medium_RMCRT(self):
		if self._dq_medium_RMCRT is None:
			print("Calculating dq_medium_RMCRT...")
			print(f"\text={self.beta}, omega={self.omega}, SF={self.SF}, g1={self.g1}, D={self.D}\n" + \
		 			f"\tT={self.T if self.T is not None else self.MR_profile}\n" + \
					f"\tsize={self.size}, nRays={self.nRays}")
			self._dq_medium_RMCRT = calc_dq_medium(
										kappa=self.kappa,
										sigma_sca=self.sigma,
										T=self.MR_profile if self.MR_profile else self.T,
										limits=[0,self.D],
										size=self.size,
										nRays=self.nRays,
										SF=self.SF,
										g1=self.g1,
										nonhomogeneous=False,
										)
			print("\tDone.")
		return self._dq_medium_RMCRT

	@property
	def dq_gas_RMCRT(self):
		if self._dq_gas_RMCRT is None:
			_, self._dq_gas_RMCRT = calc_dq(
										kappa=self.kappa,
										sigma_sca=self.sigma,
										T=0.0,
										limits=[0,self.D],
										size=self.size,
										nRays=self.nRays,
										SF=self.SF,
										g1=self.g1,
										in_rad=self.in_rad
										)
			# multiply by -1 to reverse negation enforced in read_dqrad() called by calc_dq()
			self._dq_gas_RMCRT = -np.array(self._dq_gas_RMCRT)
		return self._dq_gas_RMCRT

	@property
	def dq_gas_term(self):
		if self._dq_gas_term is None:
			self._dq_gas_term = calc_dq_ED(q=self.in_rad, 
										t=self.t, 
										beta=self.beta, 
										omega=self.omega, 
										SF=self.SF, 
										g1=self.g1, 
										rho=self.ref)
		return self._dq_gas_term
	
	@property
	def dq_medium_term(self):
		if self._dq_medium_term is None:
			dq_medium = []
			for i in range(self.size):
				q_pos = 0.0
				q_neg = 0.0
				dq_medium.append(0.0)
				if i != self.size-1: # omit last cell
					D = self.dx * (self.size-i-1) # slab thickness
					t = self.x[i+1:] - (self.x[i] + self.dx/2)
					# q_pos = self._calc_emission(t=t,
					# 							T=self.T_x[i+1:],
					# 							D=D,
					# 							tag=f"{i}pos")
					q_pos = calc_EWET(	t=t,
										T=self.T_x[i+1:],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
					rho_plus = calc_ref(thickness=D,
										ext=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1)
					rho_minus = calc_ref(thickness=self.D - D,
										ext=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1)
					q_pos /= (1 - rho_plus * rho_minus)
					# dq_medium[-1] += calc_first_cell_dq(
					# 						kappa=self.kappa,
					# 						sigma=self.sigma,
					# 						nRays=self.nRays,
					# 						SF=self.SF,
					# 						g1=self.g1,
					# 						q=q_pos,
					# 						D=(i+1)*self.dx,
					# 						size=i+1
					# 						)
					dq_medium[-1] += calc_dq_ED(q_pos,
											t=self.dx/2,
											beta=self.beta,
											omega=self.omega,
											SF=self.SF,
											g1=self.g1,
											rho=rho_minus
											)[0]
				if i != 0: # omit first cell
					D = self.dx * i
					t = (self.x[i] - self.dx/2) - self.x[0:i]
					# q_neg = self._calc_emission(t=t,
					# 							T=self.T_x[0:i],
					# 							D=D,
					# 							tag=f"{i}neg")
					q_neg = calc_EWET(	t=t,
										T=self.T_x[0:i],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
					rho_plus = calc_ref(thickness=D,
										ext=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1)
					rho_minus = calc_ref(thickness=self.D - D,
										ext=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1)
					q_neg /= (1 - rho_plus * rho_minus)
					# dq_medium[-1] += calc_first_cell_dq(
					# 						kappa=self.kappa,
					# 						sigma=self.sigma,
					# 						nRays=self.nRays,
					# 						SF=self.SF,
					# 						g1=self.g1,
					# 						q=q_neg,
					# 						D=(self.size-i)*self.dx,
					# 						size=self.size-i
					# 						)
					dq_medium[-1] += calc_dq_ED(q_neg,
											t=self.dx/2,
											beta=self.beta,
											omega=self.omega,
											SF=self.SF,
											g1=self.g1,
											rho=rho_minus
											)[0]
				# dq_medium.append(
				# 	calc_dq_ED(q_pos+q_neg,
				# 			t=self.dx/2,
				# 			beta=self.beta,
				# 			omega=self.omega,
				# 			SF=self.SF,
				# 			g1=self.g1)[0]
				# )
			self._dq_medium_term = np.array(dq_medium)
		return self._dq_medium_term
	
	def _calc_emission(self, t:np.ndarray, T:np.ndarray, D:float, tag:str):
		""" 
		Wrapper for calc_emission() array of distances t and temperatures T.
		It saves the temperature profile to a file if t is an array within a 'dump_profiles' directory.
		
		Args:
			t (array): distances from the radiating surface
			T (array): temperatures at distances t
			D (float): slab thickness

		Returns:
			emission (float): calculated emission
		"""
		# print(f"Calculating emission for i{tag} = ", end=' ')

		profileName = None
		profileDir = "dump_profiles"
		if not os.path.exists(profileDir):
			os.makedirs(profileDir)

		if not self.MR_profile:
			if T[0] != T[-1]:
				raise ValueError("T must be constant when MR_profile is not provided.")
			MR_profile = f"{T[0]:0.0f}const"
		else:
			MR_profile = self.MR_profile

		yc = D - t
		if yc.size != 1:
			# ensure yc is ascending
			if yc[0] > yc[-1]:
				yc = yc[::-1]
				T = T[::-1]
			profileName = os.path.join(profileDir, f"T{MR_profile}-i{tag}.T")
			with open(profileName, 'w') as f:
				f.write(f"#\n")
				f.write(f"#\n")
				f.write(f"#\n")
				for yi, Ti in zip(yc, T):
					f.write(f"{yi} {Ti}\n")

		outputName = f"T{MR_profile}-i{tag}-abs{self.kappa:0.0f}-sca{self.sigma:0.0f}-{self.SF}-g1{self.g1:0.3f}-D{D*1e6:0.3f}microns-size{len(t)}-nRays{self.nRays}.emi"
		emission = calc_emission(	ext=self.beta,
									omega=self.omega,
									T=profileName if profileName else T[0],
									limits=[0,D],
									size=len(t),
									nRays=self.nRays,
									SF=self.SF,
									g1=self.g1,
									outputName=outputName,
							)
		# print(emission)
		return emission

	@property
	def dq_model(self):
		if self._dq_model is None:
			# print("Calculating dq_model...")
			# print(f"\text={self.beta}, omega={self.omega}, SF={self.SF}, g1={self.g1}, D={self.D}\n" + \
		 	# 		f"\tT={self.T if self.T is not None else self.MR_profile}, in_rad={self.in_rad}\n" + \
			# 		f"\tsize={self.size}, nRays={self.nRays}")
			self._dq_model = self.dq_cooling_term - self.dq_gas_term - self.dq_medium_term
		return self._dq_model
