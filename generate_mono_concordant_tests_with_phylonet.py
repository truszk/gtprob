import dendropy
from test_lib import *
import subprocess
import os

stells_path='/home/jakub/stells/STELLS2-master/stells-v2-1-0-linux64'
DP_path='python ./calcProbConcordant.py'
phylonet_path = "java -jar /home/jakub/phylonet/PhyloNet_3.8.2.jar"
test_dir = 'tests'


random.seed(1)

for nspecies in [4,8,12,16]:
	for genespersp in [1,2,3]:
		for mean_len in [0.1,0.2,0.5]:
			stree=generate_species_tree(nspecies,mean_len)
			gtree=generate_mono_concordant_tree(stree,genespersp)

			try:
				stname=os.path.join(test_dir,'ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len))
				with open(stname,'w') as FST:
					FST.write(str(stree)+';\n')
				
				
				gtname=os.path.join(test_dir,'GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len))
				with open(gtname,'w') as FGT:
					FGT.write(str(gtree)+';\n')
				

				gtname_stells=os.path.join(test_dir,'GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'.nonums')
				with open(gtname_stells,'w') as FGTS:
					ss_gtree_stells=''.join([ i for i in str(gtree) if not i.isdigit()])
					FGTS.write(ss_gtree_stells+';\n')


			except:
				sys.stderr.write("Couldn't create output files")
				exit(0)

			#now run stells and ours
			stells_out=subprocess.check_output(stells_path+" -g "+gtname_stells+" -s "+stname+" -B -O |grep Log-prob",shell=True)
			stells_prob=float(stells_out.split()[-1])

			ours_out=subprocess.check_output(DP_path+" "+gtname+" "+stname,shell=True)
			ours_prob=float(ours_out)
			
			phylonet_cmd_path = os.path.join("phylonet",gtname+".cmd")
			make_phylonet_file(phylonet_cmd_path,gtree,stree)

			phylonet_out = subprocess.check_output(phylonet_path+" "+phylonet_cmd_path+" |grep probability",shell=True)
			phylonet_prob=float(phylonet_out.split()[-1])
			prob_diff=ours_prob-stells_prob
			prob_diff_phylonet = ours_prob - phylonet_prob
			print(gtname+" DP prob="+str(ours_prob)+" STELLS exact prob="+str(stells_prob)+" diff="+str(prob_diff))
			print(gtname+" DP prob="+str(ours_prob)+" Phylonet prob="+str(phylonet_prob)+" diff="+str(prob_diff_phylonet))





