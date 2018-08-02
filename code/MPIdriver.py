#!/software/anaconda3.6/bin/python

from EBLSST import LSSTEBworker, OpSim
import csv
import argparse
import numpy as np
from mpi4py import MPI


def define_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("-c", "--n_bin", 		type=int, help="Number of binaries per process [100000]")
	parser.add_argument("-o", "--output_file", 	type=str, help="output file name")
	parser.add_argument("-a", "--n_band", 		type=int, help="Nterms_band input for gatspy [2]")
	parser.add_argument("-b", "--n_base", 		type=int, help="Nterms_base input for gatspy [2]")
	parser.add_argument("-s", "--seed", 		type=int, help="random seed []")
	parser.add_argument("-v", "--verbose", 		action='store_true', help="Set to show verbose output")
	parser.add_argument("-l", "--opsim", 		action='store_false', help="set to run LSST OpSim, else run nobs =")

	#https://docs.python.org/2/howto/argparse.html
	args = parser.parse_args()
	#to print out the options that were selected (probably some way to use this to quickly assign args)
	opts = vars(args)
	options = { k : opts[k] for k in opts if opts[k] != None }
	print(options)

	return args

def apply_args(worker, args):


	if (args.n_bin is not None):
		worker.n_bin = args.n_bin
		
	if (args.output_file is not None):
		worker.ofile = args.output_file

	if (args.n_band is not None):
		worker.n_band = args.n_band
	if (args.n_base is not None):
		worker.n_base = args.n_base

	worker.verbose = args.verbose
	worker.doOpSim = args.opsim

	#set the random seed
	if (args.seed is not None):
		worker.seed = args.seed


if __name__ == "__main__":

	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

	sendbuf = None
	root = 0

	args = define_args()
	if (args.n_bin == None):
		args.n_bin = size

	nfields = 5292 #total number of fields from OpSim
	nfieldsPerCore = int(np.floor(nfields/size))

	sendbuf = np.empty((size, 3*nfieldsPerCore), dtype='float64')
	recvbuf = np.empty(3*nfieldsPerCore, dtype='float64')

	if (rank == 0):
		OpS = OpSim()
		OpS.dbFile = '/projects/p30137/ageller/EBLSST/input/db/minion_1016_sqlite.db' #for the OpSim database	
		OpS.getAllOpSimFields()

		#scatter the fieldID, RA, Dec 
		nfields = len(OpS.fieldID)
		#get as close as we can to having everything scattered
		maxIndex = nfieldsPerCore*size
		output = np.vstack((OpS.fieldID[:maxIndex], OpS.RA[:maxIndex], OpS.Dec[:maxIndex])).T

		print("reshaping to send to other processes")
		sendbuf = np.reshape(output, (size, 3*nfieldsPerCore))


	#scatter to the all of the processes
	comm.Scatter(sendbuf, recvbuf, root=root) 
	#now reshape again to get back to the right format
	fieldData = np.reshape(recvbuf, (nfieldsPerCore, 3))	

	print("rank", rank, fieldData)

	#add on any extra fields to rank =0
	if (rank == 0):
		if (nfieldsPerCore*size < nfields):
			print("adding to rank 0")
			extra = np.vstack((OpS.fieldID[maxIndex:], OpS.RA[maxIndex:], OpS.Dec[maxIndex:])).T
			fieldData = np.append(fieldData, extra)

	#define the worker
	worker = LSSTEBworker()
	worker.filterFilesRoot = '/projects/p30137/ageller/EBLSST/input/filters/'
	#check for command-line arguments
	apply_args(worker, args)	
	if (worker.seed == None):
		worker.seed = 1234
	worker.seed += rank

	#redefine the OpSim fieldID, RA, Dec and the run through the rest of the code
	fields = fieldData.T
	OpS = OpSim()
	OpS.dbFile = '/projects/p30137/ageller/EBLSST/input/db/minion_1016_sqlite.db' #for the OpSim database	
	OpS.getCursors()
	OpS.fieldID = fields[0]
	OpS.RA = fields[1]
	OpS.Dec = fields[2]
	worker.OpSim = OpS

	#initialize
	worker.initialize() #sets the random seed and reads in the filter files
	#worker.doLSM = False

	print(worker.BreivikGal)
	
	raise

	#set up the output file
	worker.ofile = 'output_files/'+str(rank).zfill(3) + worker.ofile
	csvfile = open(worker.ofile, 'wt')	
	worker.csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	#write header
	worker.writeOutputLine(None, header=True)

	j = 0 #this is used with multiprocessing to keep track of the return_dict, but not needed here
	for i, line in enumerate(gxDat):
		line = gxDat[i]

		#define the binary parameters
		EB = worker.getEB(line, i)
		EB.lineNum = i
		print(rank, i, EB.period)

		if (EB.observable):
			worker.return_dict[j] = EB
			worker.run_ellc_gatspy(j)
			EB = worker.return_dict[j]

		worker.writeOutputLine(EB)
		csvfile.flush()




	csvfile.close()


