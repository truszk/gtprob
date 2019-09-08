import dendropy
import sys
import random
from dendropy.simulate import treesim
import math
import numpy as np
import subprocess

stells_path='/home/jakub/stells/STELLS2-master/stells-v2-1-0-linux64'
DP_path='python /home/jakub/speciestrees/calcProbConcordant.py'

justcommandsch=True

letters=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def generateNames(nnames):
	length=int(math.ceil(math.log(nnames)/math.log(len(letters))))
	names=[]
	for i in range(nnames):
		x=i
		digits=[]
		for j in range(length):
			digits.append(x%len(letters))
			x=x/len(letters)
		assert len(digits)==length
		lettered=map(lambda d:letters[d],digits)
		name=''.join(lettered)
		names.append(name)
	return names




#assigns genes to species. A gene is assigned to a species if the gene name contains the species name as a substring. 
#If one gene has many species names as substrings or if some genes are left unassigned, the function returns None.
def assignGenesToSpecies(G,S):
	spToGenes={}
	splabels=[leaf.taxon.label for leaf in S.leaf_nodes()]
	
	for species in splabels:
		spToGenes[species]=[]

	genelabels=[leaf.taxon.label for leaf in G.leaf_nodes()]

	for gene in genelabels:
		assigned=False
		for species in splabels:
			if species in gene:
				if assigned==True:
					return None
				else:
					spToGenes[species].append(gene)
					assigned=True

		if assigned==False: #still no species for this gene
			return None
	return spToGenes




def isConcordant(G,S,spToGenes):
	G2=dendropy.Tree(G)
	for s in S.leaf_nodes():
		genes=spToGenes[s.taxon.label]
		mrca=G2.mrca(taxon_labels=genes)
		leaf_taxa=map(lambda x: x.taxon.label,mrca.leaf_nodes())
		if set(leaf_taxa)!=set(genes):
			return False
		children=get_children(mrca)
		for child in children:
			mrca.remove_child(child)

		mrca.taxon=G.taxon_namespace.get_taxon(label=s.taxon.label)
		

	# maybe update bipartitions before?
	
	S.encode_bipartitions()
	G2.encode_bipartitions()
	diff=treecompare.symmetric_difference(G2,S)
	if diff==0:
		return True
	else:
		return False


def splitLeaf(node):
	assert node.is_leaf()
	c1=dendropy.Node()
	c2=dendropy.Node()
	node.add_child(c1)
	node.add_child(c2)
	return c1,c2


def generateMonoConcordantTree(sptree,genesPerSp):

	gtree=dendropy.Tree(sptree)
	for node in gtree.postorder_node_iter():
		node.edge_length=None
	for leaf in gtree.leaf_nodes():
		name=leaf.taxon.label
		leaves=[leaf]
		for i in range(genesPerSp-1):
			toSplit=random.choice(range(len(leaves)))
			c1,c2=splitLeaf(leaves[toSplit])
			del leaves[toSplit]
			leaves.append(c1)
			leaves.append(c2)
		for i in range(len(leaves)):
			tx=dendropy.Taxon(label=name+str(i))
			leaves[i].taxon=tx

	gtree.unassign_taxa(exclude_leaves=True)
	return gtree

def sample_edge_length(mean_len):
	l=np.random.exponential(mean_len)
	return l

def generateSpeciesTree(nspecies,mean_edge_length):
	node=dendropy.Node()
	leaves=[node]
	for i in range(nspecies-1):
		toSplit=random.choice(range(len(leaves)))
		c1,c2=splitLeaf(leaves[toSplit])
		c1.edge_length=sample_edge_length(mean_edge_length)
		c2.edge_length=sample_edge_length(mean_edge_length)
		del leaves[toSplit]
		leaves.append(c1)
		leaves.append(c2)

	names=generateNames(nspecies)
	for i in range(len(leaves)):
		tx=dendropy.Taxon(label=names[i])
		leaves[i].taxon=tx
	tree=dendropy.Tree(seed_node=node)
	return tree


random.seed(1)


#few samples, many species

for nspecies in [4,8,16,24,32,40]:
	for genespersp in [1,3,5]:
		for mean_len in [0.2]:
			for ind in range(10):
				if justcommandsch:
					if nspecies>16:
						continue
					if nspecies>4 and genespersp>1:
						continue
					gtname_stells='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)+'.nonums'
					stname='ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					gtname='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					
					print 'qsub -cwd -b y -l "hostname=fischer0'+str(random.randint(1,4))+'" ./runStellsCH.sh '+gtname_stells+" "+stname
					print "sleep 3"
					continue
				stree=generateSpeciesTree(nspecies,mean_len)
				gtree=generateMonoConcordantTree(stree,genespersp)
				try:
					stname='ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					FST=file(stname,'w')
					print >>FST,str(stree)+';'
					FST.close()
					
					gtname='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					FGT=file(gtname,'w')
					print >>FGT,str(gtree)+';'
					FGT.close()

					gtname_stells='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)+'.nonums'
					FGTS=file(gtname_stells,'w')
					ss_gtree_stells=''.join([ i for i in str(gtree) if not i.isdigit()])
					print >>FGTS,ss_gtree_stells+';'
					FGTS.close()

				except:
					print >>sys.stderr,"Couldn't create output files"
					exit(0)

#many samples few species

for nspecies in [2,3,4]:
	for genespersp in [5,10,15,20,50]:
		for mean_len in [0.2]:
			for ind in range(10):
				if justcommandsch:
					gtname_stells='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)+'.nonums'
					stname='ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					gtname='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					
					print 'qsub -cwd -b y -l "hostname=fischer0'+str(random.randint(1,4))+'" ./runStellsCH.sh '+gtname_stells+" "+stname
					print "sleep 3"
					continue
				stree=generateSpeciesTree(nspecies,mean_len)
				gtree=generateMonoConcordantTree(stree,genespersp)

				try:
					stname='ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					FST=file(stname,'w')
					print >>FST,str(stree)+';'
					FST.close()
					
					gtname='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)
					FGT=file(gtname,'w')
					print >>FGT,str(gtree)+';'
					FGT.close()

					gtname_stells='GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len)+'ind'+str(ind)+'.nonums'
					FGTS=file(gtname_stells,'w')
					ss_gtree_stells=''.join([ i for i in str(gtree) if not i.isdigit()])
					print >>FGTS,ss_gtree_stells+';'
					FGTS.close()

				except:
					print >>sys.stderr,"Couldn't create output files"
					exit(0)

			#now run stells and ours
			#stells_out=subprocess.check_output(stells_path+" -g "+gtname_stells+" -s "+stname+" -B -O |grep Log-prob",shell=True)
			#stells_prob=float(stells_out.split()[-1])
			#ours_out=subprocess.check_output(DP_path+" "+gtname+" "+stname,shell=True)
			#ours_prob=float(ours_out)
			#prob_diff=ours_prob-stells_prob
			#print gtname+" DP prob="+str(ours_prob)+" STELLS exact prob="+str(stells_prob)+" diff="+str(prob_diff)





