import pandas as pd
import numpy as np

#for Quest
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import matplotlib.lines as mlines

def plotHist(histAll, histObs, histRec, bin_edges, histAllOD, histObsOD, histRecOD, bin_edgesOD, xlim, xtitle, fname):
	c1 = '#0294A5'  #turqoise
	c2 = '#d95f02' #orange from color brewer
	c3 = '#00353E' #slate
	c4 = '#508201' #olive

	f,(ax1, ax2, ax3) = plt.subplots(3,1,figsize=(5, 12), sharex=True)

	binHalf = (bin_edges[1] - bin_edges[0])/2.
	binHalfOD = (bin_edgesOD[1] - bin_edgesOD[0])/2.

	#CDF
	cdfAll = []
	cdfObs = []
	cdfRec = []
	cdfAllOD = []
	cdfObsOD = []
	cdfRecOD = []
	for i in range(len(histAll)):
		cdfAll.append(np.sum(histAll[:i])/np.sum(histAll))
	for i in range(len(histAllOD)):
		cdfAllOD.append(np.sum(histAllOD[:i])/np.sum(histAllOD))
	for i in range(len(histObs)):
		cdfObs.append(np.sum(histObs[:i])/np.sum(histObs))
	for i in range(len(histObsOD)):
		cdfObsOD.append(np.sum(histObsOD[:i])/np.sum(histObsOD))
	for i in range(len(histRec)):
		cdfRec.append(np.sum(histRec[:i])/np.sum(histRec))
	for i in range(len(histRecOD)):
		cdfRecOD.append(np.sum(histRecOD[:i])/np.sum(histRecOD))
	ax1.step(bin_edges, cdfAll, color=c1, label='All')
	ax1.step(bin_edges, cdfObs, color=c2, label='Observable')
	ax1.step(bin_edges, cdfRec, color=c3, label='Recoverable')
	ax1.step(bin_edgesOD, cdfAllOD, color=c1, linestyle=':')
	ax1.step(bin_edgesOD, cdfObsOD, color=c2, linestyle=':')
	ax1.step(bin_edgesOD, cdfRecOD, color=c3, linestyle=':')
	ax1.set_ylabel('CDF', fontsize=16)
	ax1.set_ylim(-0.01,1.01)
	ax1.set_xlim(xlim[0],xlim[1])

	#PDF
	ax2.step(bin_edges, histAll/np.sum(histAll), color=c1, label='All')
	ax2.step(bin_edges, histObs/np.sum(histObs), color=c2, label='Observable')
	ax2.step(bin_edges, histRec/np.sum(histRec), color=c3, label='Recoverable')
	ax2.step(bin_edgesOD, histAllOD/np.sum(histAllOD), color=c1, linestyle=':')
	ax2.step(bin_edgesOD, histObsOD/np.sum(histObsOD), color=c2, linestyle=':')
	ax2.step(bin_edgesOD, histRecOD/np.sum(histRecOD), color=c3, linestyle=':')
	ax2.set_ylabel('PDF', fontsize=16)
	ax2.set_yscale('log')
	ax2.set_ylim(0.5e-5, 0.9)
	ax2.set_xlim(xlim[0],xlim[1])

	#ratio (I don't need these use bits anymore, but will keep just in case.  If I set > 0 then I miss the first bin with the line, presumably because of the step plotting method)
	use = np.where(histAll > -1)[0]
	ax3.step(bin_edges[use], histObs[use]/histAll[use], color=c2, label='Observable/All')
	ax3.plot(bin_edges[use] - binHalf, histObs[use]/histAll[use], 'o',color=c1, markersize=5, markeredgecolor=c2)

	ax3.step(bin_edges[use], histRec[use]/histAll[use], color=c3, label='Recoverable/All')
	ax3.plot(bin_edges[use] - binHalf, histRec[use]/histAll[use], 'o',color=c1, markersize=5, markeredgecolor=c3)

	use = np.where(histObs > -1)[0]
	ax3.step(bin_edges[use], histRec[use]/histObs[use], color=c3, label='Recoverable/Observable')
	ax3.plot(bin_edges[use] - binHalf, histRec[use]/histObs[use], 'o',color=c2, markersize=5, markeredgecolor=c3)

	use = np.where(histAllOD > -1)[0]
	ax3.step(bin_edgesOD[use], histObsOD[use]/histAllOD[use], color=c2, linestyle=':')
	ax3.plot(bin_edgesOD[use] - binHalfOD, histObsOD[use]/histAllOD[use], 'o',color=c1, markersize=3.5, markeredgecolor=c2)

	ax3.step(bin_edgesOD[use], histRecOD[use]/histAllOD[use], color=c3, linestyle=':')
	ax3.plot(bin_edgesOD[use] - binHalfOD, histRecOD[use]/histAllOD[use], 'o',color=c1, markersize=3.5, markeredgecolor=c3)

	use = np.where(histObsOD > -1)[0]
	ax3.step(bin_edgesOD[use], histRecOD[use]/histObsOD[use], color=c3, linestyle=':')
	ax3.plot(bin_edgesOD[use] - binHalfOD, histRecOD[use]/histObsOD[use], 'o',color=c2, markersize=3.5, markeredgecolor=c3)


	ax3.set_ylabel('Ratio', fontsize=16)
	ax3.set_yscale('log')
	ax3.set_ylim(10**-4,1)
	ax3.set_xlabel(xtitle, fontsize=16)
	ax3.set_xlim(xlim[0],xlim[1])

	# ax1.legend()
	# ax2.legend()
	# ax3.legend()
	lAll = mlines.Line2D([], [], color=c1, label='All')
	lObs = mlines.Line2D([], [], color=c2, label='Obs.')
	lRec = mlines.Line2D([], [], color=c3, label='Rec.')
	lObsAll = mlines.Line2D([], [], color=c2, marker='o', markerfacecolor=c1, markersize=5, markeredgecolor=c2, label='Obs./All')
	lRecAll = mlines.Line2D([], [], color=c3, marker='o', markerfacecolor=c1, markersize=5, markeredgecolor=c3, label='Rec./All')
	lRecObs = mlines.Line2D([], [], color=c3, marker='o', markerfacecolor=c2, markersize=5, markeredgecolor=c3, label='Rec./Obs.')
	ax1.legend(handles=[lAll, lObs, lRec, lObsAll, lRecAll, lRecObs], loc='lower right')

	f.subplots_adjust(hspace=0)
	f.savefig(fname+'_ratio.pdf',format='pdf', bbox_inches = 'tight')



if __name__ == "__main__":

	data = pd.read_csv('EBLSST_lphist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_lphist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[-2, 6], 'log(Period [days])', 'EBLSST_lphist_final')

	data = pd.read_csv('EBLSST_ehist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_ehist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[0, 1], 'Eccentricity', 'EBLSST_ehist_final')

	data = pd.read_csv('EBLSST_dhist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_dhist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[0, 20], 'Distance (kpc)', 'EBLSST_dhist_final')

	data = pd.read_csv('EBLSST_m1hist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_m1hist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[0, 2], r'm$_1$ (M$_\odot$)', 'EBLSST_m1hist_final')

	data = pd.read_csv('EBLSST_qhist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_qhist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[0, 3], r'q (m$_2$/m$_1$)', 'EBLSST_qhist_final')

	data = pd.read_csv('EBLSST_maghist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_maghist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[11, 24], 'mag', 'EBLSST_maghist_final')

	data = pd.read_csv('EBLSST_rhist.csv')
	dataOD = pd.read_csv('obsDist/EBLSST_rhist.csv')
	plotHist(data['histAll'].values, data['histObs'].values, data['histRec'].values, data['binEdges'].values, \
			dataOD['histAll'].values, dataOD['histObs'].values, dataOD['histRec'].values, dataOD['binEdges'].values, \
			[0, 3], r'r$_2$/r$_1$', 'EBLSST_rhist_final')