#!/software/anaconda3.6/bin/python

from mpi4py import MPI
import os
import pickle

from OpSim import OpSim
from astropy.coordinates import SkyCoord
from astropy import units
import numpy as np


if __name__ == "__main__":

	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

	sendbuf = None
	root = 0


	#nfields = 2295# "primary" fields
	nfields = 3339# all observed fields
	#nfields = 5292 #total number of fields from OpSim
	nfieldsPerCore = int(np.floor(nfields/size))
	print(f"nfields={nfields}, nfieldsPerCore={nfieldsPerCore}")

	sendbuf = np.empty((size, nfieldsPerCore), dtype='float64')
	recvbuf = np.empty(nfieldsPerCore, dtype='float64')

	if (rank == root):

		if not os.path.exists('output_files'):
			os.makedirs('output_files')

		OpS = OpSim()
		OpS.dbFile = '/projects/p30137/ageller/EBLSST/input/db/minion_1016_sqlite.db' #for the OpSim database	
		OpS.getAllOpSimFields()

		#primary = np.where(OpS.Nobs > 800)
		#OpS.fieldID = OpS.fieldID[primary]
		observed = np.where(OpS.Nobs > 0)
		OpS.fieldID = OpS.fieldID[observed]

		nfields = len(OpS.fieldID)
		print(f"rank 0 nfields={nfields}")
		print(OpS.fieldID)

		#scatter the fieldID 
		#get as close as we can to having everything scattered
		maxIndex = min(nfieldsPerCore*size, nfields-1)
		output = OpS.fieldID[:maxIndex].T

		print("reshaping to send to other processes")
		sendbuf = np.reshape(output, (size, nfieldsPerCore))

	#scatter to the all of the processes
	comm.Scatter(sendbuf, recvbuf, root=root) 
	#now reshape again to get back to the right format
	fieldData = np.reshape(recvbuf, (nfieldsPerCore, 1))	

	#add on any extra fields to rank =0
	if (rank == 0):
		if (nfieldsPerCore*size < nfields):
			print("adding to rank 0")
			extra = OpS.fieldID[maxIndex:].T
			fieldData = np.append(fieldData, extra)

	#redefine the OpSim fieldID and the run through the rest of the code
	fields = fieldData.T
	OpS = OpSim()
	OpS.dbFile = '/projects/p30137/ageller/EBLSST/input/db/minion_1016_sqlite.db' #for the OpSim database	
	OpS.getCursors()
	OpS.fieldID = np.squeeze(fields)
	OpS.obsDates = np.full_like(OpS.fieldID, dict(), dtype=dict)
	OpS.NobsDates = np.full_like(OpS.fieldID, dict(), dtype=dict)
	OpS.m_5 = np.full_like(OpS.fieldID, dict(), dtype=dict)
	OpS.totalNobs = np.full_like(OpS.fieldID, 0)

	filters = ['u_', 'g_', 'r_', 'i_', 'z_', 'y_']
	Nbins = 100
	OpSimi = fieldData.astype("int")
	OpS.dt = np.array([{} for i in OpSimi])
	for i, ID in enumerate(OpSimi): 
		print(rank, i, ID)
		dt = {}
		OpS.setDates(i, filters)
		for f in filters:
			dt[f] = np.diff(OpS.obsDates[i][f])
		OpS.dt[i] = dt

		dist = {}
		for f in filters:
			pdf,bin_edges = np.histogram(np.log10(OpS.dt[i][f], where=(OpS.dt[i][f] > 0)), bins=Nbins)

			bins = bin_edges[:-1] + np.diff(bin_edges)/2.
			cdf = np.cumsum(pdf)
			dist[f] = {}
			dist[f]['bins']= bins
			dist[f]['cdf']= cdf
			dist[f]['pdf']= pdf

		oname = 'output_files/'+str(int(OpS.fieldID[i])).zfill(4) + "dist.pickle"
		pickle.dump( dist, open( oname, "wb" ) )
