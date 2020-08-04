import dendropy
from test_lib import *
import subprocess
import os

stells_path='/home/jakub/stells/STELLS2-master/stells-v2-1-0-linux64'
DP_path='./calc_prob_concordant.py'
phylonet_path = "java -jar /home/jakub/phylonet/PhyloNet_3.8.2.jar"
test_dir = 'appendix'
phylonet_cmd_prefix = ''

sp_examples = ['(((a:1.0,b:1.0):0.5,c:1.5):0.3,(d:0.9,e:0.9):0.9);',
		'(((a:0.1,b:0.1):0.1,(c:0.15,d:0.15):0.05):0.05,((e:0.02,f:0.02):0.13,(g:0.04,h:0.04):0.11):0.1);',
		'(((((((a:0.03,b:0.03):0.03,c:0.06):0.03,d:0.09):0.03,e:0.12):0.03,f:0.15):0.03,g:0.18):0.03,h:0.21);',]

dataset_names = ['small','balanced','caterpillar']

generate_gt = False

def get_output_and_time(arglist):
	arglist_with_time = ['time']+arglist
	p = subprocess.Popen(arglist_with_time,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out = p.stdout.readlines()
	err = p.stderr.readlines()
	timestr = err[0].split(b'user')[0]
	time = float(timestr)
	return out,time

random.seed(1)

for name in dataset_names:
	for genespersp in [1,4]:
		stname = os.path.join(test_dir,'ST.'+name+str(genespersp))
		stree = dendropy.Tree.get(path = stname, schema = 'newick')
		
		if generate_gt:
			gtree=generate_mono_concordant_tree(stree,genespersp)
		else:			
			gtname=os.path.join(test_dir,'GT.'+name+str(genespersp))
			gtree=dendropy.Tree.get(path=gtname,schema='newick')

		try:

			stname=os.path.join(test_dir,'ST.'+name+str(genespersp))
			gtname=os.path.join(test_dir,'GT.'+name+str(genespersp))
			gtname_stells=os.path.join(test_dir,'GT.'+name+str(genespersp)+'.nonums')
			if generate_gt:

				with open(stname,'w') as FST:
					FST.write(str(stree)+';\n')
				
				with open(gtname,'w') as FGT:
					FGT.write(str(gtree)+';\n')

			with open(gtname_stells,'w') as FGTS:
				ss_gtree_stells=''.join([ i for i in str(gtree) if not i.isdigit()])
				FGTS.write(ss_gtree_stells+';\n')


		except:
			sys.stderr.write("Couldn't create output files")
			exit(0)

		#now run stells and ours
		stells_out, stells_time = get_output_and_time([stells_path,"-g",gtname_stells,"-s",stname,"-B","-O"])
		stells_out_line = stells_out[6]
		stells_prob=float(stells_out_line.split()[-1])



		ours_out, ours_time=get_output_and_time(['python',DP_path,gtname,stname])
		ours_prob=float(ours_out[0])
		
		phylonet_cmd_path = os.path.join(phylonet_cmd_prefix,gtname+".cmd")
		make_phylonet_file(phylonet_cmd_path,gtree,stree)

		#phylonet_out = subprocess.check_output(phylonet_path+" "+phylonet_cmd_path+" |grep probability",shell=True)
		phylonet_out_whole, phylonet_time = get_output_and_time(phylonet_path.split()+[phylonet_cmd_path])
		phylonet_out = phylonet_out_whole[-1]
		phylonet_prob=float(phylonet_out.split()[-1])

		prob_diff=ours_prob-stells_prob
		prob_diff_phylonet = ours_prob - phylonet_prob
		print(gtname+" DP prob="+str(ours_prob)+" STELLS exact prob="+str(stells_prob)+" diff="+str(prob_diff))
		print(gtname+" DP prob="+str(ours_prob)+" Phylonet prob="+str(phylonet_prob)+" diff="+str(prob_diff_phylonet))
		print(gtname+" DP time = "+str(ours_time)+" STELLS time = "+str(stells_time)+" Phylonet time = "+str(phylonet_time))

		





