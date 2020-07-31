import dendropy
import sys
import random
from dendropy.simulate import treesim
import math
import numpy as np
from test_lib import *


try:
	nspecies=int(sys.argv[1])
	genespersp=int(sys.argv[2])
	mean_len=float(sys.argv[3])
except:
	sys.stderr.write("Usage: "+sys.argv[0]+" n_species n_genes_per_species mean_branch_len_in_ST")
	exit(0)

stree=generateSpeciesTree(nspecies,mean_len)
gtree=generateMonoConcordantTree(stree,genespersp)

try:
	st_filename = 'ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)
	with open(st_filename,'w') as FST:
		FST.write(str(stree)+';\n')

	gt_filename = 'GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)
	with open(gt_filename,'w') as FGT:
		FGT.write(str(gtree)+';\n')

except:
	sys.stderr.write("Couldn't create output files")
	exit(0)

