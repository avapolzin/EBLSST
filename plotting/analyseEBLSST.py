######################
#NEED TO ACCOUNT FOR THE BINARY FRACTION when combining histograms
#####################

import pandas as pd
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units

#for Quest
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

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
	f.savefig(fname,format='pdf', bbox_inches = 'tight')

if __name__ == "__main__":

	#cutoff in percent error for "recovered"
	Pcut = 0.1

	#minimum number of lines to consider in file
	Nlim = 3

	#bins for all the histograms
	mbins = np.linspace(0,2, 100)
	qbins = np.linspace(0,2, 100)
	ebins = np.linspace(0,1, 100)
	lpbins = np.linspace(-2, 6, 100)
	dbins = np.linspace(0, 20, 100)

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

	#Read in all the data and make the histograms
	d = "output_files/"
	files = os.listdir(d)
	IDs = []
	for i, f in enumerate(files):
		print(round(i/len(files),4), f)

		#read in the header
		header = pd.read_csv(d+f, nrows=1)
		#Nmult = header['NstarsTRILEGAL'][0]
		Nmult = 1.

		RA.append(header['OpSimRA'])
		Dec.append(header['OpSimDec'])

		#read in rest of the file
		data = pd.read_csv(d+f, header = 2)
		rF = 0.
		rN = 0.
		if (len(data.index) >= Nlim):
			#create histograms
			#All
			m1hAll0, m1b = np.histogram(data["m1"], bins=mbins, density=True)
			qhAll0, qb = np.histogram(data["m2"]/data["m1"], bins=qbins, density=True)
			ehAll0, eb = np.histogram(data["e"], bins=ebins, density=True)
			lphAll0, lpb = np.histogram(np.log10(data["p"]), bins=lpbins, density=True)
			dhAll0, db = np.histogram(data["d"], bins=dbins, density=True)
			m1hAll += m1hAll0*Nmult
			qhAll += qhAll0*Nmult
			ehAll += ehAll0*Nmult
			lphAll += lphAll0*Nmult
			dhAll += dhAll0*Nmult

			#Obs
			obs = data.loc[data['LSM_PERIOD'] != -999]
			if (len(obs.index) >= Nlim):
				m1hObs0, m1b = np.histogram(obs["m1"], bins=mbins, density=True)
				qhObs0, qb = np.histogram(obs["m2"]/obs["m1"], bins=qbins, density=True)
				ehObs0, eb = np.histogram(obs["e"], bins=ebins, density=True)
				lphObs0, lpb = np.histogram(np.log10(obs["p"]), bins=lpbins, density=True)
				dhObs0, db = np.histogram(obs["d"], bins=dbins, density=True)
				m1hObs += m1hObs0*Nmult
				qhObs += qhObs0*Nmult
				ehObs += ehObs0*Nmult
				lphObs += lphObs0*Nmult
				dhObs += dhObs0*Nmult

				#Rec
				fullP = abs(data['LSM_PERIOD'] - data['p'])/data['LSM_PERIOD']
				halfP = abs(0.5*data['LSM_PERIOD'] - data['p'])/(0.5*data['LSM_PERIOD'])
				twiceP = abs(2.*data['LSM_PERIOD'] - data['p'])/(2.*data['LSM_PERIOD'])
				rec = data.loc[(data['LSM_PERIOD'] != -999) & ( (fullP < Pcut) | (halfP < Pcut) | (twiceP < Pcut))]
				if (len(rec.index) >= Nlim):
					m1hRec0, m1b = np.histogram(rec["m1"], bins=mbins, density=True)
					qhRec0, qb = np.histogram(rec["m2"]/rec["m1"], bins=qbins, density=True)
					ehRec0, eb = np.histogram(rec["e"], bins=ebins, density=True)
					lphRec0, lpb = np.histogram(np.log10(rec["p"]), bins=lpbins, density=True)
					dhRec0, db = np.histogram(rec["d"], bins=dbins, density=True)
					m1hRec += m1hRec0*Nmult
					qhRec += qhRec0*Nmult
					ehRec += ehRec0*Nmult
					lphRec += lphRec0*Nmult
					dhRec += dhRec0*Nmult

					#for the mollweide
					rF = len(rec.index)/len(data.index)
					rN = len(rec.index)/len(data.index)*Nmult


		recFrac.append(rF)
		recN.append(rN)

	#plot and save the histograms
	saveHist(np.append(m1hAll,0), np.append(m1hObs,0), np.append(m1hRec,0), m1b, 'm1 (Msolar)', 'EBLSST_m1hist.pdf')
	saveHist(np.append(qhAll,0), np.append(qhObs,0), np.append(qhRec,0), qb, 'q (m2/m1)', 'EBLSST_qhist.pdf')
	saveHist(np.append(ehAll,0), np.append(ehObs,0), np.append(ehRec,0), eb, 'e', 'EBLSST_ehist.pdf')
	saveHist(np.append(lphAll,0), np.append(lphObs,0), np.append(lphRec,0), lpb, 'log(P [days])', 'EBLSST_lphist.pdf')
	saveHist(np.append(dhAll,0), np.append(dhObs,0), np.append(dhRec,0), db, 'd (kpc)', 'EBLSST_dhist.pdf')


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
	mlw = ax.scatter(np.array(RAwrap).ravel()*np.pi/180., np.array(Decwrap).ravel()*np.pi/180., c=np.log10(np.array(recFrac)*100.), cmap='viridis_r', s = 4)
	cbar = f.colorbar(mlw, shrink=0.7)
	cbar.set_label(r'log10(% recovered)')
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
