import numpy as np
from astropy import units, constants
import pandas as pd

from matplotlib import pyplot as plt

class SEDClass(object):
	def __init__(self, *args,**kwargs):
		self.filterFilesRoot = '../input/filters/'
		#['u_band_Response.dat','g_band_Response.dat','r_band_Response.dat','i_band_Response.dat','z_band_Response.dat','y_band_Response.dat']
		self.filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']
		self.filterThroughput = dict()

		self.BCf = dict()

	def readFilters(self):
		for f in self.filters:
			# #https://github.com/lsst-pst/syseng_throughputs/tree/master/components/camera/filters
			# fname = self.filterFilesRoot + f + 'band_Response.dat'
			# df = pd.read_csv(fname, delim_whitespace=True, header=None, names=['w','t'])

			#https://github.com/lsst/throughputs/tree/master/baseline
			fname = self.filterFilesRoot + 'filter_'+f[0]+'.dat' 
			df = pd.read_csv(fname, delim_whitespace=True, header=None, names=['w','t'], skiprows = 6)
			self.filterThroughput[f] = {'w':df['w'].values*units.nm, 't':df['t'].values}

	def getBCf(self, T, dw = 0.1, wmin =10, wmax = 10000):

		#could improve this with a Kurucz model
		def blackbody(w, T):
			#erg/s/cm^2/AA/steradian
			Bl = 2.*constants.h*constants.c**2./w**5. / (np.exp( constants.h*constants.c / (w*constants.k_B*T)) -1.)
			return Bl.decompose()

		w = np.arange(wmin, wmax, dw)*units.nm
		f = blackbody(w,T).value
		#norm = np.sum(f)*dw
		norm = max(f)
		fn = f/norm
		
		ftot = np.sum(fn)

		for f in self.filters:
			ft = np.interp(w, self.filterThroughput[f]['w'], self.filterThroughput[f]['t'])
			# plt.plot(w,ft,'--')
			tot = np.sum(fn*ft)
			self.BCf[f] = tot/ftot

		print(self.BCf)

		# plt.plot(w,fn)
		# for f in self.filters:
		# 	plt.plot(self.filterThroughput[f]['w'], self.filterThroughput[f]['t'])
		# plt.xlim(100, 2000)
		# plt.show()

if __name__ == "__main__":


	sed = SEDClass()
	sed.readFilters()
	sed.getBCf(5000*units.K)
