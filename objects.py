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
		self._dq_RMCRT = None
		self._dq_model = None
		self._dq_cooling_term = None
		self._dq_gas_term = None
		self._dq_medium_term = None

	@property
	def dq_RMCRT(self):
		if self._dq_RMCRT is None:
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
			self._dq_cooling_term = 4 * self.kappa * SIGMA * np.power(self.T_x, 4)
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
			dt = self.D/(len(self.xMR)) # MR cell thickness
			dq_medium_xMR = []
			for i in range(len(self.xMR)):
				q_pos = 0.0
				q_neg = 0.0
				if i != len(self.xMR)-1: # omit last cell
					D = dt * (len(self.xMR)-i-1) # slab thickness
					t = self.xMR[i+1:] - (self.xMR[i] + dt/2)
					q_pos = calc_EWET(	t=t,
										T=self.T_xMR[i+1:],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
				if i != 0: # omit first cell
					D = dt * i
					t = (self.xMR[i] - dt/2) - self.xMR[0:i]
					q_neg = calc_EWET(	t=t,
										T=self.T_xMR[0:i],
										beta=self.beta,
										omega=self.omega,
										SF=self.SF,
										g1=self.g1,
										D=D)
				dq_medium_xMR.append(
					calc_ED(q_pos+q_neg,
							t=dt/2,
							beta=self.beta,
							omega=self.omega,
							SF=self.SF,
							g1=self.g1,
							D=dt)
				)
			dq_medium_xMR = np.array(dq_medium_xMR)
			dq_medium_x = np.interp(self.x, self.xMR, dq_medium_xMR) # this should not be linear interpolation
			self._dq_medium_term = dq_medium_x
		return self._dq_medium_term

	@property
	def dq_model(self):
		if self._dq_model is None:
			self._dq_model = self.dq_cooling_term - self.dq_gas_term - self.dq_medium_term
		return self._dq_model