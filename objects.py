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
			  T:float|str, size:int, nRays:int, in_rad:float, size_RMCRT:int=None):
		super().__init__(beta=beta, omega=omega, SF=SF, g1=g1, D=D)

		self.size = size
		self.size_RMCRT = size_RMCRT if size_RMCRT else size
		self.dx = self.D / self.size
		self.x = np.linspace(self.dx/2, self.D-self.dx/2, self.size)
		self.dx_RMCRT = self.D / self.size_RMCRT
		self._x_RMCRT = None
		self.xnorm = self.x / self.D
		self._xnorm_RMCRT = None
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

		self._dq_RMCRT = None # corresponds to self.size_RMCRT
		self._dq_model_xRMCRT = None
		# All are func of self.x (not xMR), corresponds to self.size
		self._dq_model = None
		self._dq_cooling_term = None
		self._dq_gas_term = None
		self._dq_medium_term = None

	@property
	def x_RMCRT(self):
		if self._x_RMCRT is None:
			self._x_RMCRT = np.linspace(self.dx_RMCRT/2, self.D-self.dx_RMCRT/2, self.size_RMCRT)
		return self._x_RMCRT
	
	@property
	def xnorm_RMCRT(self):
		if self._xnorm_RMCRT is None:
			self._xnorm_RMCRT = self.x_RMCRT / self.D
		return self._xnorm_RMCRT

	@property
	def dq_RMCRT(self):
		if self._dq_RMCRT is None:
			print("Calculating dq_RMCRT...")
			print(f"\text={self.beta}, omega={self.omega}, SF={self.SF}, g1={self.g1}, D={self.D}\n" + \
		 			f"\tT={self.T if self.T is not None else self.MR_profile}, in_rad={self.in_rad}\n" + \
					f"\tsize={self.size_RMCRT}, nRays={self.nRays}")	
			_, self._dq_RMCRT = calc_dq(kappa=self.kappa,
											sigma_sca=self.sigma,
											T=self.MR_profile if self.MR_profile else self.T,
											limits=[0,self.D],
											size=self.size_RMCRT,
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
	def dq_gas_term(self):
		if self._dq_gas_term is None:
			self._dq_gas_term = calc_ED(q=self.in_rad, 
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
				if i != self.size-1: # omit last cell
					D = self.dx * (self.size-i-1) # slab thickness
					t = self.x[i+1:] - (self.x[i] + self.dx/2)
					q_pos = calc_EWET(	t=t,
										T=self.T_x[i+1:],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
				if i != 0: # omit first cell
					D = self.dx * i
					t = (self.x[i] - self.dx/2) - self.x[0:i]
					q_neg = calc_EWET(	t=t,
										T=self.T_x[0:i],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
				dq_medium.append(
					calc_ED(q_pos+q_neg,
							t=self.dx/2,
							beta=self.beta,
							omega=self.omega,
							SF=self.SF,
							g1=self.g1,
							D=self.dx)[0]
				)
			self._dq_medium_term = np.array(dq_medium)
		return self._dq_medium_term

	@property
	def dq_model(self):
		if self._dq_model is None:
			# print("Calculating dq_model...")
			# print(f"\text={self.beta}, omega={self.omega}, SF={self.SF}, g1={self.g1}, D={self.D}\n" + \
		 	# 		f"\tT={self.T if self.T is not None else self.MR_profile}, in_rad={self.in_rad}\n" + \
			# 		f"\tsize={self.size}, nRays={self.nRays}")
			self._dq_model = self.dq_cooling_term - self.dq_gas_term - self.dq_medium_term
		return self._dq_model
	
	@property
	def dq_model_xRMCRT(self):
		if self._dq_model_xRMCRT is None:
			if self.size_RMCRT > self.size:
				raise ValueError("size_RMCRT must be less than or equal to size")
			self._dq_model_xRMCRT = np.interp(self.x_RMCRT, self.x, self.dq_model)
		return self._dq_model_xRMCRT