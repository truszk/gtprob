import dendropy
from dendropy.calculate import treecompare
import numpy as np
import scipy.linalg
import math
from math import log,exp
import sys

err_not_concordant=1000000
err_threshold_probcoals=0.00001
infty_time=1000000.0

delta=0.001

loggedsum_range=500.0 # if a prob is exp(500) times smaller than the largest element, we drop it

probhash={}
pathhash={}
lpathhash={}
fac_table=[]
choosehash={}
lchoosehash={}
probspechash={}
lprobspechash={}

def loggedsum(logprobs_o):
	global loggedsum_range
	if len(logprobs_o)==0:
		return float('-inf')
	threshold=max(logprobs_o)-loggedsum_range
	if threshold==float('-inf'):
		return float('-inf')
	logprobs=filter(lambda x: x>=threshold,logprobs_o)
	lprobsscaled=map(lambda x:x-threshold,logprobs)
	probsscaled=map(exp,lprobsscaled)
	# now we can sort them to help avoid numerical underflows...should we?
	probssorted=sorted(probsscaled)
	esum=sum(probssorted)
	lsum=log(esum)+threshold
	return lsum

# sums adjacent pairs in list
# not log scale
def sum_adjacent_pairs(lst):
	pairsums=[]
	for i in range(len(lst)):
		if i%2==0:
			pairsums.append(lst[i])
		else:
			pairsums[-1]+=lst[i]
	return pairsums

# sums alternating positive-negative numbers by first summing pairs of adjacent numbers to improve precision
def sum_list_by_pairs(lst):
	pairsums=sum_adjacent_pairs(lst)
	total=sum(pairsums)
	return total



def precalc_log_factorials(n):
	global fac_table
	acc=log(1.0)#log of 0 factorial
	for i in range(n):
		fac_table.append(acc)
		acc+=log(i+1)

def logfact(n):
	return fac_table[n]

def log_mult_consecutive(smaller,larger):
	assert smaller-1<=larger
	return logfact(larger)-logfact(smaller-1)

# n choose k
def choose(n,k):
	global choosehash
	if (n,k) in choosehash:
		return choosehash[(n,k)]
	chooz=math.factorial(n)/math.factorial(k)/math.factorial(n-k)
	choosehash[(n,k)]=chooz
	return chooz

def log_choose(n,k):
	global lchoosehash
	if (n,k) in lchoosehash:
		return lchoosehash[(n,k)]
	lchooz=logfact(n)-logfact(k)-logfact(n-k)
	lchoosehash[(n,k)]=lchooz
	return lchooz

#computes smaller*(smaller+1)*...*larger=larger!/(smaller-1)!
def mult_consecutive(smaller,larger):
	assert smaller-1<=larger
	return math.factorial(larger)/math.factorial(smaller-1)

#returns a list of child nodes, not sure if needed, probably depends on the version of dendropy
def get_children(node):
	return [child for child in node.child_nodes()]


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

#for monoconcordance, checking if the species names agree in the two files
def isSameSpeciesTree(S,spToGeneNums):
	for s in S.leaf_nodes():
		if not s.taxon.label in spToGeneNums:
			return False
	if len(spToGeneNums)!=len(S.leaf_nodes()):
		return False
	return True





# calculates the probability of the gene tree G given species tree S. G has to be monophyletically concordant with S.
# G - gene tree
# S - species tree
# spToGenes - a map with a list of samples for each species 
def calcProbConcordant(G,S,spToGenes):
	if not isConcordant(G,S,spToGenes):
		return err_not_concordant
	for s in S.leaf_node_iter():
		genes=spToGenes[s.taxon.label]
		mrca=G.mrca(taxon_labels=genes)
		Gmrca=dendropy.Tree(seed_node=mrca.extract_subtree())
		prob=geneTreeProbOneSpecies(Gmrca,infty_time)# this also annotates Gmrca with n_int_nodes and n_leaf_nodes
		ncoals=Gmrca.seed_node.n_int_nodes
		s.P=[0.0 for i in range(ncoals+1)]
		s.P[ncoals]=prob
		#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)

	nodecnt=0

	for s in S.postorder_node_iter():
		if s.is_leaf():
			continue
		else:
			child=get_children(s)
			s.P=[0.0 for i in range(len(child[0].P)+len(child[1].P))]
			for k in range(len(s.P)):
				for k1p in range(0,k):
					
					k2p=k-1-k1p

					for k1 in range(k1p,len(child[0].P)):
						for k2 in range(k2p,len(child[1].P)):
							#print "nodecnt="+str(nodecnt)+"k,k1,k2="+str(k)+','+str(k1)+','+str(k2)
							partialmain=child[0].P[k1]*child[1].P[k2]*choose(k1p+k2p,k1p)
							partialc1=probSpecificCoalPath(k1+1,k1p+1,child[0].edge_length)/probSpecificCoalPath(k1+1,1,infty_time)#remove some coals from c1
							partialc2=probSpecificCoalPath(k2+1,k2p+1,child[1].edge_length)/probSpecificCoalPath(k2+1,1,infty_time)#remove coals from c2
							partialtop=probSpecificCoalPath(k+1,1,infty_time)#add coals to the top species
							partial=partialmain*partialc1*partialc2*partialtop
							s.P[k]+=partial
			#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)
		nodecnt+=1
	#sum everything to get the overall prob

	totalprob=sum(S.seed_node.P)
	return totalprob



# calculates the log-probability of the gene tree G given species tree S. G has to be monophyletically concordant with S.
# G - gene tree
# S - species tree
# spToGenes - a map with a list of samples for each species 
def logProbConcordant(G,S,spToGenes):
	if not isConcordant(G,S,spToGenes):
		return err_not_concordant
	for s in S.leaf_node_iter():
		genes=spToGenes[s.taxon.label]
		mrca=G.mrca(taxon_labels=genes)
		Gmrca=dendropy.Tree(seed_node=mrca.extract_subtree())
		prob=geneTreeProbOneSpecies(Gmrca,infty_time)# this also annotates Gmrca with n_int_nodes and n_leaf_nodes
		ncoals=Gmrca.seed_node.n_int_nodes
		s.P=[float('-inf') for i in range(ncoals+1)]
		s.P[ncoals]=log(prob)
		#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)

	nodecnt=0

	for s in S.postorder_node_iter():
		if s.is_leaf():
			continue
		else:
			child=get_children(s)
			s.P=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k in range(len(s.P)):
				partials=[]
				for k1p in range(0,k):
					
					k2p=k-1-k1p

					for k1 in range(k1p,len(child[0].P)):
						for k2 in range(k2p,len(child[1].P)):
							#print "nodecnt="+str(nodecnt)+"k,k1,k2="+str(k)+','+str(k1)+','+str(k2)
							partialmain=child[0].P[k1]+child[1].P[k2]+log_choose(k1p+k2p,k1p)
							partialc1=logProbSpecificCoalPath(k1+1,k1p+1,child[0].edge_length)-logProbSpecificCoalPath(k1+1,1,infty_time)#remove some coals from c1
							partialc2=logProbSpecificCoalPath(k2+1,k2p+1,child[1].edge_length)-logProbSpecificCoalPath(k2+1,1,infty_time)#remove coals from c2
							partialtop=logProbSpecificCoalPath(k+1,1,infty_time)#add coals to the top species
							partial=partialmain+partialc1+partialc2+partialtop
							assert not math.isnan(partial)
							partials.append(partial)
				s.P[k]=loggedsum(partials)
				assert not math.isnan(s.P[k])
			#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)
		nodecnt+=1
	#sum everything to get the overall prob

	totallogprob=loggedsum(S.seed_node.P)
	return totallogprob


# calculates the log-probability of the gene tree G given species tree S. G has to be monophyletically concordant with S.
# G - gene tree
# S - species tree
# spToGenes - a map with a list of samples for each species 
def logProbConcordantFast(G,S,spToGenes):
	if not isConcordant(G,S,spToGenes):
		return err_not_concordant
	for s in S.leaf_node_iter():
		genes=spToGenes[s.taxon.label]
		mrca=G.mrca(taxon_labels=genes)
		Gmrca=dendropy.Tree(seed_node=mrca.extract_subtree())
		lprob=geneTreeLogProbOneSpecies(Gmrca,infty_time)# this also annotates Gmrca with n_int_nodes and n_leaf_nodes
		ncoals=Gmrca.seed_node.n_int_nodes
		s.P=[float('-inf') for i in range(ncoals+1)]
		s.P[ncoals]=lprob
		#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)

	nodecnt=0

	partialhash={}

	for s in S.postorder_node_iter():
		if s.is_leaf():
			continue
		else:
			child=get_children(s)
			
			
			partialsums1=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k1p in range(0,len(child[0].P)):
				
				partials1=[]

				for k1 in range(k1p,len(child[0].P)):
					partialmain1=child[0].P[k1]
					partialc1=logProbSpecificCoalPath(k1+1,k1p+1,child[0].edge_length)-logProbSpecificCoalPath(k1+1,1,infty_time)
					partial1=partialmain1+partialc1
					partials1.append(partial1)

				partialsum1=loggedsum(partials1)
				partialsums1[k1p]=partialsum1

			partialsums2=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k2p in range(0,len(child[1].P)):

				partials2=[]
				
				for k2 in range(k2p,len(child[1].P)):
					partialmain2=child[1].P[k2]
					partialc2=logProbSpecificCoalPath(k2+1,k2p+1,child[1].edge_length)-logProbSpecificCoalPath(k2+1,1,infty_time)
					partial2=partialmain2+partialc2
					partials2.append(partial2)

				partialsum2=loggedsum(partials2)
				partialsums2[k2p]=partialsum2







			s.P=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k in range(len(s.P)):
				partials=[]
				for k1p in range(0,k):
					
					k2p=k-1-k1p
					partial=partialsums1[k1p]+partialsums2[k2p]+log_choose(k1p+k2p,k1p)+logProbSpecificCoalPath(k+1,1,infty_time)
					partials.append(partial)
				
					
				s.P[k]=loggedsum(partials)
				assert not math.isnan(s.P[k])
			#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)
		nodecnt+=1
	#sum everything to get the overall prob

	totallogprob=loggedsum(S.seed_node.P)
	return totallogprob




# calculates the log-probability of the gene tree G given species tree S. G has to be monophyletically concordant with S.
# G - gene tree
# S - species tree
# spToGenes - a map with a list of samples for each species 
def logProbMonoConcordance(S,spToGeneNums):
	if not isSameSpeciesTree(S,spToGeneNums):
		return err_not_concordant
	for s in S.leaf_node_iter():
		ngenes=spToGeneNums[s.taxon.label]
		prob=1.0# this is the prob. of all possible topos under infty_time
		ncoals=ngenes-1
		s.P=[float('-inf') for i in range(ncoals+1)]
		s.P[ncoals]=log(prob)
		#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)

	nodecnt=0

	partialhash={}

	for s in S.postorder_node_iter():
		if s.is_leaf():
			continue
		else:
			child=get_children(s)
			
			
			partialsums1=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k1p in range(0,len(child[0].P)):
				
				partials1=[]

				for k1 in range(k1p,len(child[0].P)):
					partialmain1=child[0].P[k1]
					partialc1=logProbSpecificCoalPath(k1+1,k1p+1,child[0].edge_length)-logProbSpecificCoalPath(k1+1,1,infty_time)
					partial1=partialmain1+partialc1
					partials1.append(partial1)

				partialsum1=loggedsum(partials1)
				partialsums1[k1p]=partialsum1

			partialsums2=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k2p in range(0,len(child[1].P)):

				partials2=[]
				
				for k2 in range(k2p,len(child[1].P)):
					partialmain2=child[1].P[k2]
					partialc2=logProbSpecificCoalPath(k2+1,k2p+1,child[1].edge_length)-logProbSpecificCoalPath(k2+1,1,infty_time)
					partial2=partialmain2+partialc2
					partials2.append(partial2)

				partialsum2=loggedsum(partials2)
				partialsums2[k2p]=partialsum2



			s.P=[float('-inf') for i in range(len(child[0].P)+len(child[1].P))]
			
			for k in range(len(s.P)):
				partials=[]
				for k1p in range(0,k):
					
					k2p=k-1-k1p
					partial=partialsums1[k1p]+partialsums2[k2p]+log_choose(k1p+k2p,k1p)+logProbSpecificCoalPath(k+1,1,infty_time)
					partials.append(partial)
				
					
				s.P[k]=loggedsum(partials)
				assert not math.isnan(s.P[k])
			#print str(map(lambda x:x.taxon.label,s.leaf_nodes()))+' s.P=   '+str(s.P)
		nodecnt+=1
	#sum everything to get the overall prob

	totallogprob=loggedsum(S.seed_node.P)
	return totallogprob








def geneTreeProbOneSpecies(G,length):
	nranked=getNumRankedTopos(G)
	pranked=probSpecificCoalPath(len(G.leaf_nodes()),1,length)
	prob=pranked*nranked
	return prob

def geneTreeLogProbOneSpecies(G,length):
	nranked=getNumRankedTopos(G)
	lpranked=logProbSpecificCoalPath(len(G.leaf_nodes()),1,length)
	lprob=lpranked+log(nranked)
	return lprob

# we assume a constant delta
def getCoalProbFromIntegral(density,nup,length):
	global delta

	lastepoch=nup
	lambdaa=float(lastepoch*(lastepoch-1)/2)

	probtable=[]
	for i in range(len(density)):
		xi=delta*i
		complement=length-xi
		probnocoal=exp(-lambdaa*complement)
		probtable.append(density[i]*delta*probnocoal)

	prob=sum(probtable)
	return prob


	
def probCoalsByIntegration(ndown,nup,length):
	global denshash
	global delta

	assert ndown>=nup

	npoints=int(length/delta)
	density=[[0.0 for i in range(npoints)] for j in range(ndown-nup)]

	#initialize exponential
	lambdaa=float(ndown*(ndown-1)/2)
	for i in range(npoints):
		xi=delta*i
		density[0][i]=lambdaa*exp(-lambdaa*xi)

	probnocoal=exp(-lambdaa*length)
	#probhash[(ndown,ndown,length)]=probnocoal

	dnum=0
	for epoch in range(ndown-1,nup,-1):
		probcoal=getCoalProbFromIntegral(density[dnum],epoch,length)
		#probhash[(ndown,epoch,length)]=probcoal
		dnum+=1
		lambdaa=float(epoch*(epoch-1)/2)
		for i in range(npoints):
			xi=delta*i
			for j in range(i):
				xj=delta*j
				complement=xi-xj
				densexp=lambdaa*exp(-lambdaa*complement)
				density[dnum][i]+=density[dnum-1][j]*densexp*delta
	probcoal=getCoalProbFromIntegral(density[dnum],nup,length)
	#probhash[(ndown,nup,length)]=probcoal
	return probcoal


# builds rate matrix for probCoals - treats the number of lineages as a continuous-time Markov chain
def buildRateMatrix(ndown,length):
	ratemat=np.eye(ndown)
	for i in range(ndown):
		nlins=i+1
		rate=length*nlins*(nlins-1)*0.5
		if nlins==1:
			ratemat[nlins-1,nlins-1]=0.0 # absorbing state
		else:
			ratemat[nlins-1,nlins-2]=rate
			ratemat[nlins-1,nlins-1]=-rate
	return ratemat

def getProbCoalMatrix(ndown,length):
	ratemat=buildRateMatrix(ndown,length)
	pcoalmat=scipy.linalg.expm(ratemat)

	for i in range(ndown):# this is to cope with numerical issues which still seem to happen occasionally
		for j in range(i+1):
			if pcoalmat[i,j]<0.0:
				pcoalmat[i,j]=0.0
	return pcoalmat


def probCoals(ndown,nup,length):
	global probhash
	global err_threshold_probcoals
	
	assert ndown>=nup

	if (ndown,nup,length) in probhash:
		return probhash[(ndown,nup,length)]

	if length==None:
		length=infty_time #for the top node that continues forever, usually

	if length==infty_time:
		if nup==1:
			prob=1.0
		else:
			prob=0.0
		probhash[(ndown,nup,length)]=prob
		return prob

	pmat=getProbCoalMatrix(ndown,length)
	prob=pmat[ndown-1,nup-1]

	for i in range(ndown):
		for j in range(i+1):
			probhash[(i+1,j+1,length)]=pmat[i,j]
	assert prob>=0.0
	return prob


def probCoalsOld(ndown,nup,length):
	global probhash
	global err_threshold_probcoals
	
	assert ndown>=nup

	if (ndown,nup,length) in probhash:
		return probhash[(ndown,nup,length)]

	if length==None:
		length=infty_time #for the top node that continues forever, usually
	
	prob=0.0
	probtable=[]
	for k in range(nup,ndown+1):
		sign=((-1)**(k-nup))
		#partial=math.exp(-length*k*(k-1)/2)*(2*k-1)*((-1)**(k-nup))
		#partial/=(math.factorial(nup)*math.factorial(k-nup)*(nup+k-1))
		lpartial=-length*k*(k-1)/2.0+log(2*k-1)
		lpartial=lpartial-logfact(nup)-logfact(k-nup)-log(nup+k-1)
		#prodpartial=(float(mult_consecutive(ndown-k+1,ndown))/mult_consecutive(ndown,ndown+k-1))*mult_consecutive(nup,nup+k-1)
		lprodpartial=log_mult_consecutive(ndown-k+1,ndown)-log_mult_consecutive(ndown,ndown+k-1)+log_mult_consecutive(nup,nup+k-1)
		lpartial+=lprodpartial
		partial=sign*exp(lpartial)
		probtable.append(partial)
		#add all that to prob
		#print "ndown="+str(ndown)+"nup="+str(nup)+"length="+str(length)+"k="+str(k)+"-->  partial="+str(partial)
		#prob+=partial

	prob=sum_list_by_pairs(probtable)

	if prob<0.0 and prob>-err_threshold_probcoals:
		prob=0.0
	probhash[(ndown,nup,length)]=prob
	#print probtable
	#summ=0.0
	#for i in range(len(probtable)-1,-1,-1):
	#	summ+=probtable[i]
	#print summ

	return prob

# not used, perhaps will come back to that
def probCoalsWithLoggedSum(ndown,nup,length):
	global probhash
	global err_threshold_probcoals
	
	assert ndown>=nup

	if (ndown,nup,length) in probhash:
		return probhash[(ndown,nup,length)]

	if length==None:
		length=infty_time #for the top node that continues forever, usually
	
	prob=0.0
	poslpartials=[]
	neglpartials=[]
	for k in range(nup,ndown+1):
		sign=((-1)**(k-nup))
		#partial=math.exp(-length*k*(k-1)/2)*(2*k-1)*((-1)**(k-nup))
		#partial/=(math.factorial(nup)*math.factorial(k-nup)*(nup+k-1))
		lpartial=-length*k*(k-1)/2.0+log(2*k-1)
		lpartial=lpartial-logfact(nup)-logfact(k-nup)-log(nup+k-1)
		#prodpartial=(float(mult_consecutive(ndown-k+1,ndown))/mult_consecutive(ndown,ndown+k-1))*mult_consecutive(nup,nup+k-1)
		lprodpartial=log_mult_consecutive(ndown-k+1,ndown)-log_mult_consecutive(ndown,ndown+k-1)+log_mult_consecutive(nup,nup+k-1)
		lpartial+=lprodpartial
		if sign>0:
			poslpartials.append(lpartial)
		else:
			neglpartials.append(lpartial)
		#partial=sign*exp(lpartial)
		#probtable.append(partial)
		#add all that to prob
		#print "ndown="+str(ndown)+"nup="+str(nup)+"length="+str(length)+"k="+str(k)+"-->  partial="+str(partial)
		#prob+=partial
	poslprob=loggedsum(poslpartials)
	neglprob=loggedsum(neglpartials)
	prob=exp(poslprob)-exp(neglprob)

	if prob<0.0 and prob>-err_threshold_probcoals:
		prob=0.0
	probhash[(ndown,nup,length)]=prob
	#print probtable
	#summ=0.0
	#for i in range(len(probtable)-1,-1,-1):
	#	summ+=probtable[i]
	#print summ
	#assert prob>=0.0
	return prob



def numCoalPaths(ndown,nup):
	global pathhash
	assert ndown>=nup
	if (ndown,nup) in pathhash:
		return pathhash[(ndown,nup)]
	#npaths=math.factorial(ndown)*math.factorial(ndown-1)*2**(nup-ndown)/(math.factorial(nup)*math.factorial(nup-1))
	npaths=(mult_consecutive(nup,ndown-1))**2*ndown*2**(nup-ndown)/nup

	pathhash[(ndown,nup)]=npaths
	return npaths

def logNumCoalPaths(ndown,nup):
	global lpathhash
	assert ndown>=nup
	if (ndown,nup) in lpathhash:
		return lpathhash[(ndown,nup)]
	lnpaths=log_mult_consecutive(nup,ndown-1)*2.0+log(ndown)-log(nup)+log(2.0)*(nup-ndown)
	lpathhash[(ndown,nup)]=lnpaths
	return lnpaths

def probSpecificCoalPath(ndown,nup,length):
	global probspechash
	if (ndown,nup,length) in probspechash:
		return probspechash[(ndown,nup,length)]
	p=probCoals(ndown,nup,length)
	if p==0.0:
		return 0.0
	lpspec=log(p)-logNumCoalPaths(ndown,nup)
	pspec=exp(lpspec)
	probspechash[(ndown,nup,length)]=pspec
	return pspec

def logProbSpecificCoalPath(ndown,nup,length):
	global lprobspechash
	if (ndown,nup,length) in lprobspechash:
		return lprobspechash[(ndown,nup,length)]
	p=probCoals(ndown,nup,length)
	lpspec=0.0
	if p==0.0:
		lpspec=float('-inf')-logNumCoalPaths(ndown,nup)
	else:
		lpspec=log(p)-logNumCoalPaths(ndown,nup)
	lprobspechash[(ndown,nup,length)]=lpspec
	return lpspec




# the number of ranked topologies corresponding to the unranked topology G
def getNumRankedTopos(G):
	annotateNumCoals(G)
	for node in G.postorder_node_iter():
		if node.is_leaf():
			node.n_ranked=1
		else:
			children=get_children(node)
			assert len(children)==2
			node.n_ranked=choose(node.n_int_nodes-1,children[0].n_int_nodes)*children[0].n_ranked*children[1].n_ranked

	return G.seed_node.n_ranked




#for each node, records # internal nodes down that node (including the node itself)
def annotateNumCoals(G):
	
	for node in G.postorder_node_iter():
		if node.is_leaf():
			node.n_int_nodes=0
			node.n_leaf_nodes=1
		else:
			children=get_children(node)
			assert len(children)==2
			node.n_int_nodes=children[0].n_int_nodes+children[1].n_int_nodes+1
			node.n_leaf_nodes=children[0].n_leaf_nodes+children[1].n_leaf_nodes


#for probability of concordance
def parseSamplesFile(lines):
	spToGeneNums={}
	for line in lines:
		tokens=line.split()
		nsamples=int(tokens[1])
		spToGeneNums[tokens[0]]=nsamples
	return spToGeneNums



#-----tests----

def testProbCoals(length=0.006,ndown=70):
	probs=[]
	for i in range(1,ndown+1):
		prob=probCoals(ndown,i,length)
		probint=probCoalsByIntegration(ndown,i,length)
		print("prob Coals"+str(prob)+" "+str(probint))
		probs.append(prob)

	summ=sum(probs)

	print(summ)
	return


	

#def testAllTopos(ngenes,length):

#-------MAIN--------------


tree_option=True

try:
	stree=dendropy.Tree.get(path=sys.argv[2],schema='newick',rooting='force-rooted')
	with open(sys.argv[1],'r') as gfile:
		lines=gfile.readlines()
		tokens=lines[0].split()

	if len(tokens)>1 and tokens[1].isdigit():
		tree_option=False
	else:
		gtree=dendropy.Tree.get(path=sys.argv[1],schema='newick',taxon_namespace=stree.taxon_namespace,rooting='force-rooted')

except:
	print('Usage: '+str(sys.argv[0])+' genetreefile speciestreefile')
	print('or')
	print(str(sys.argv[0])+' genesamplesfile speciestreefile')
	exit(0)

if tree_option:

	precalc_log_factorials(len(gtree.leaf_nodes())*2+100)

	spToGenes=assignGenesToSpecies(gtree,stree)
	if spToGenes==None:
		sys.stderr.write("Error: Cannot figure out the assignment of gene tree leaves to species tree leaves - each gene tree leaf name must have exactly one species name as a substring")
		exit(0)



	lprob=logProbConcordantFast(gtree,stree,spToGenes)
	#prob=calcProbConcordant(gtree,stree,spToGenes)
	if lprob==err_not_concordant:
		sys.stderr.write("Error: The gene tree is not monophyletically concordant with the species tree.")
	else:
		#print prob
		print(lprob)
else: # calc monophyletic concordance
	
	spToGeneNums=parseSamplesFile(lines)
	nsamples=sum([spToGeneNums[s] for s in spToGeneNums])

	precalc_log_factorials(nsamples*2+100)

	lprob=logProbMonoConcordance(stree,spToGeneNums)
	if lprob==err_not_concordant:
		sys.stderr.write("Error: The sample file does not match the species tree.")
	else:
		#print prob
		print(lprob)
	






