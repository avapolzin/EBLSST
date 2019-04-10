

import pandas as pd
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.modeling import models, fitting

#for Quest
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

def fitRagfb():
	x = [0.05, 0.1, 1, 8, 15]  #estimates of midpoints in bins, and using this: https://sites.uni.edu/morgans/astro/course/Notes/section2/spectralmasses.html
	y = [0.20, 0.35, 0.50, 0.70, 0.75]
	init = models.PowerLaw1D(amplitude=0.5, x_0=1, alpha=-1.)
	fitter = fitting.LevMarLSQFitter()
	fit = fitter(init, x, y)

	return fit

def saveHist(histAll, histObs, histRec, bin_edges, xtitle, fname):
	c1 = '#0294A5'  #turqoise
	c2 = '#d95f02' #orange from color brewer
	c3 = '#00353E' #slate
	f,(ax1, ax2) = plt.subplots(2,1,figsize=(5, 8), sharex=True)

	#PDF
	ax1.step(bin_edges, histAll/np.sum(histAll), color=c1)
	ax1.step(bin_edges, histObs/np.sum(histObs), color=c2)
	ax1.step(bin_edges, histRec/np.sum(histRec), color=c3)
	ax1.set_ylabel('PDF')
	ax1.set_yscale('log')
	#CDF
	cdfAll = []
	cdfObs = []
	cdfRec = []
	for i in range(len(histAll)):
		cdfAll.append(np.sum(histAll[:i])/np.sum(histAll))
	for i in range(len(histObs)):
		cdfObs.append(np.sum(histObs[:i])/np.sum(histObs))
	for i in range(len(histRec)):
		cdfRec.append(np.sum(histRec[:i])/np.sum(histRec))
	ax2.step(bin_edges, cdfAll, color=c1)
	ax2.step(bin_edges, cdfObs, color=c2)
	ax2.step(bin_edges, cdfRec, color=c3)
	ax2.set_ylabel('CDF')

	ax2.set_xlabel(xtitle)
	f.subplots_adjust(hspace=0)
	f.savefig(fname+'.pdf',format='pdf', bbox_inches = 'tight')

	#write to a text file
	with open(fname+'.txt','w') as f:
		f.write('binEdges, histAll, histObs, histRec\n')
		for (b,a,o,r) in zip(bin_edges, histAll, histObs, histRec):
			f.write(str(b)+','+str(a)+','+str(o)+','+str(r)+'\n')

if __name__ == "__main__":

	#get the Raghavan binary fraction fit
	fbFit= fitRagfb()
	print(fbFit)
		
	#cutoff in percent error for "recovered"
	Pcut = 0.1

	#minimum number of lines to consider in file
	Nlim = 3

	#bins for all the histograms
	Nbins = 50
	mbins = np.linspace(0,5, Nbins)
	qbins = np.linspace(0,2, Nbins)
	ebins = np.linspace(0,1, Nbins)
	lpbins = np.linspace(-2, 6, Nbins)
	dbins = np.linspace(0, 20, Nbins)

	#blanks for the histograms
	#All
	m1hAll = np.zeros_like(mbins)[1:]
	qhAll = np.zeros_like(qbins)[1:]
	ehAll = np.zeros_like(ebins)[1:]
	lphAll = np.zeros_like(lpbins)[1:]
	dhAll = np.zeros_like(dbins)[1:]
	#Observable
	m1hObs = np.zeros_like(mbins)[1:]
	qhObs = np.zeros_like(qbins)[1:]
	ehObs = np.zeros_like(ebins)[1:]
	lphObs = np.zeros_like(lpbins)[1:]
	dhObs = np.zeros_like(dbins)[1:]
	#Recovered
	m1hRec = np.zeros_like(mbins)[1:]
	qhRec = np.zeros_like(qbins)[1:]
	ehRec = np.zeros_like(ebins)[1:]
	lphRec = np.zeros_like(lpbins)[1:]
	dhRec = np.zeros_like(dbins)[1:]

	RA = []
	Dec = []
	recFrac = []
	recN = []

	rawN = []
	obsN = []
	fileN = []
	fileObsN = []
	fileRecN = []

	#Read in all the data and make the histograms
	d = "../output_files/"
	files = os.listdir(d)
	IDs = []
	for i, f in enumerate(files):
		print(round(i/len(files),4), f)

		#read in the header
		header = pd.read_csv(d+f, nrows=1)
######################
#NEED TO ACCOUNT FOR THE BINARY FRACTION when combining histograms
#####################
# Also, this weights those near the galactic plane sooo highly (and these are usually poorly recovered), that the resulting histograms are VERY noisy (since we're basically just looking at a few fields new to galactic plane)
		Nmult = header['NstarsTRILEGAL'][0]
		#Nmult = 1.

		RA.append(header['OpSimRA'])
		Dec.append(header['OpSimDec'])

		#read in rest of the file
		data = pd.read_csv(d+f, header = 2).dropna()
		rF = 0.
		rN = 0.
		if (len(data.index) >= Nlim):
			#create histograms
			#All
			m1hAll0, m1b = np.histogram(data["m1"], bins=mbins, density=True)
			qhAll0, qb = np.histogram(data["m2"]/data["m1"], bins=qbins, density=True)
			ehAll0, eb = np.histogram(data["e"], bins=ebins, density=True)
			lphAll0, lpb = np.histogram(np.ma.log10(data["p"].values).filled(-999), bins=lpbins, density=True)
			dhAll0, db = np.histogram(data["d"], bins=dbins, density=True)

			#account for the binary fraction, as a function of mass
			dm1 = np.diff(m1b)
			m1val = m1b[:-1] + dm1/2.
			fb = np.sum(m1hAll0*dm1*fbFit(m1val))
			Nmult *= fb

						
			m1hAll += m1hAll0*Nmult
			qhAll += qhAll0*Nmult
			ehAll += ehAll0*Nmult
			lphAll += lphAll0*Nmult
			dhAll += dhAll0*Nmult

			#Obs
			obs = data.loc[data['LSM_PERIOD'] != -999]
			if (len(obs.index) >= Nlim):
				ofrac = len(obs.index)/len(data.index)
				m1hObs0, m1b = np.histogram(obs["m1"], bins=mbins, density=True)
				qhObs0, qb = np.histogram(obs["m2"]/obs["m1"], bins=qbins, density=True)
				ehObs0, eb = np.histogram(obs["e"], bins=ebins, density=True)
				lphObs0, lpb = np.histogram(np.ma.log10(obs["p"].values).filled(-999), bins=lpbins, density=True)
				dhObs0, db = np.histogram(obs["d"], bins=dbins, density=True)
				m1hObs += m1hObs0*Nmult*ofrac
				qhObs += qhObs0*Nmult*ofrac
				ehObs += ehObs0*Nmult*ofrac
				lphObs += lphObs0*Nmult*ofrac
				dhObs += dhObs0*Nmult*ofrac

				#Rec
				fullP = abs(data['LSM_PERIOD'] - data['p'])/data['LSM_PERIOD']
				halfP = abs(0.5*data['LSM_PERIOD'] - data['p'])/(0.5*data['LSM_PERIOD'])
				twiceP = abs(2.*data['LSM_PERIOD'] - data['p'])/(2.*data['LSM_PERIOD'])
				rec = data.loc[(data['LSM_PERIOD'] != -999) & ( (fullP < Pcut) | (halfP < Pcut) | (twiceP < Pcut))]
				if (len(rec.index) >= Nlim):
					rfrac = len(rec.index)/len(data.index)
					m1hRec0, m1b = np.histogram(rec["m1"], bins=mbins, density=True)
					qhRec0, qb = np.histogram(rec["m2"]/rec["m1"], bins=qbins, density=True)
					ehRec0, eb = np.histogram(rec["e"], bins=ebins, density=True)
					lphRec0, lpb = np.histogram(np.ma.log10(rec["p"].values).filled(-999), bins=lpbins, density=True)
					dhRec0, db = np.histogram(rec["d"], bins=dbins, density=True)
					m1hRec += m1hRec0*Nmult*rfrac
					qhRec += qhRec0*Nmult*rfrac
					ehRec += ehRec0*Nmult*rfrac
					lphRec += lphRec0*Nmult*rfrac
					dhRec += dhRec0*Nmult*rfrac

					#for the mollweide
					rF = len(rec.index)/len(data.index)
					rN = len(rec.index)/len(data.index)*Nmult
					raN = Nmult
					obN = len(obs.index)/len(data.index)*Nmult
					fiN = len(data.index)
					fioN = len(obs.index)
					firN = len(rec.index)

		recFrac.append(rF)
		recN.append(rN)
		rawN.append(raN)
		obsN.append(obN)
		fileN.append(fiN)
		fileObsN.append(fioN)
		fileRecN.append(firN)

	#plot and save the histograms
	saveHist(np.append(m1hAll,0), np.append(m1hObs,0), np.append(m1hRec,0), m1b, 'm1 (Msolar)', 'EBLSST_m1hist')
	saveHist(np.append(qhAll,0), np.append(qhObs,0), np.append(qhRec,0), qb, 'q (m2/m1)', 'EBLSST_qhist')
	saveHist(np.append(ehAll,0), np.append(ehObs,0), np.append(ehRec,0), eb, 'e', 'EBLSST_ehist')
	saveHist(np.append(lphAll,0), np.append(lphObs,0), np.append(lphRec,0), lpb, 'log(P [days])', 'EBLSST_lphist')
	saveHist(np.append(dhAll,0), np.append(dhObs,0), np.append(dhRec,0), db, 'd (kpc)', 'EBLSST_dhist')


	#make the mollweide
	coords = SkyCoord(RA, Dec, unit=(units.degree, units.degree),frame='icrs')	
	lGal = coords.galactic.l.wrap_at(180.*units.degree).degree
	bGal = coords.galactic.b.wrap_at(180.*units.degree).degree
	RAwrap = coords.ra.wrap_at(180.*units.degree).degree
	Decwrap = coords.dec.wrap_at(180.*units.degree).degree

	f, ax = plt.subplots(subplot_kw={'projection': "mollweide"}, figsize=(8,5))
	ax.grid(True)
	#ax.set_xlabel(r"$l$",fontsize=16)
	#ax.set_ylabel(r"$b$",fontsize=16)
	#mlw = ax.scatter(lGal.ravel()*np.pi/180., bGal.ravel()*np.pi/180., c=np.log10(np.array(recFrac)*100.), cmap='viridis_r', s = 4)
	ax.set_xlabel("RA",fontsize=16)
	ax.set_ylabel("Dec",fontsize=16)
	mlw = ax.scatter(np.array(RAwrap).ravel()*np.pi/180., np.array(Decwrap).ravel()*np.pi/180., c=np.array(recFrac)*100., cmap='viridis_r', s = 4)
	cbar = f.colorbar(mlw, shrink=0.7)
	cbar.set_label(r'% recovered')
	f.savefig('mollweide_pct.pdf',format='pdf', bbox_inches = 'tight')

	f, ax = plt.subplots(subplot_kw={'projection': "mollweide"}, figsize=(8,5))
	ax.grid(True)
	#ax.set_xlabel(r"$l$",fontsize=16)
	#ax.set_ylabel(r"$b$",fontsize=16)
	#mlw = ax.scatter(lGal.ravel()*np.pi/180., bGal.ravel()*np.pi/180., c=np.log10(np.array(recN)), cmap='viridis_r', s = 4)
	ax.set_xlabel("RA",fontsize=16)
	ax.set_ylabel("Dec",fontsize=16)
	mlw = ax.scatter(np.array(RAwrap).ravel()*np.pi/180., np.array(Decwrap).ravel()*np.pi/180., c=np.log10(np.array(recN)), cmap='viridis_r', s = 4)
	cbar = f.colorbar(mlw, shrink=0.7)
	cbar.set_label(r'log10(N) recovered')
	f.savefig('mollweide_N.pdf',format='pdf', bbox_inches = 'tight')

	print("###################")
	print("number of binaries in input files (raw, log):",np.sum(fileN), np.log10(np.sum(fileN)))
	print("number of binaries in tested with gatspy (raw, log):",np.sum(fileObsN), np.log10(np.sum(fileObsN)))
	print("number of binaries in recovered with gatspy (raw, log):",np.sum(fileRecN), np.log10(np.sum(fileRecN)))
	print("###################")
	print("total in sample (raw, log):",np.sum(rawN), np.log10(np.sum(rawN)))
	print("total observable (raw, log):",np.sum(obsN), np.log10(np.sum(obsN)))
	print("total recovered (raw, log):",np.sum(recN), np.log10(np.sum(recN)))
