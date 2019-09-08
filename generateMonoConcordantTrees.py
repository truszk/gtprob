import dendropy
import sys
import random
from dendropy.simulate import treesim
import math
import numpy as np

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


try:
	nspecies=int(sys.argv[1])
	genespersp=int(sys.argv[2])
	mean_len=float(sys.argv[3])
except:
	print >>sys.stderr,"Usage: "+sys.argv[0]+" n_species n_genes_per_species mean_branch_len_in_ST"
	exit(0)

stree=generateSpeciesTree(nspecies,mean_len)
gtree=generateMonoConcordantTree(stree,genespersp)

try:
	FST=file('ST.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len),'w')
	print >>FST,str(stree)+';'
	FST.close()
	
	FGT=file('GT.sp'+str(nspecies)+'.gps'+str(genespersp)+'.len'+str(mean_len),'w')
	print >>FGT,str(gtree)+';'
	FGT.close()

except:
	print >>sys.stderr,"Couldn't create output files"
	exit(0)

