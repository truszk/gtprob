import sys
import math
from math import log,exp
import os
import subprocess
import numpy as np

def get_res_params(filename):
	tokens=filename.split('.')
	nsp=int(tokens[1][2:])
	ngps=int(tokens[2][3:])
	bl=float(str(tokens[3][-1])+'.'+tokens[4].split('ind')[0])
	ind=int(tokens[4][-1])
	stells=('nonums' in filename)*1
	if ('outch' in filename) or ('timech' in filename):
		stells=2
	return (nsp,ngps,bl,stells,ind)


def get_time(timefile):
	F=file(timefile,'r')
	lines=F.readlines()
	rtime=float(lines[0].split('user')[0])
	return rtime

def print_result(HDL,filename,lprob,time):
	params=get_res_params(filename)
	outstr=''
	for param in params:
		outstr+=(str(param)+' ')
	outstr+=(str(lprob)+' '+str(time))
	print >>HDL,outstr

def process_result(ldic,tdic,filename,lprob,time):
	params=get_res_params(filename)
	if not params[:4] in ldic:
		ldic[params[:4]]=[]
		tdic[params[:4]]=[]
	ldic[params[:4]].append(lprob) #without the index to do means
	tdic[params[:4]].append(time)


def stformat(num):
	if math.isnan(num) or num==0.0:
		return "0.000"
	ndigits=max(int(math.log10(abs(num))),0)
	nafter=max([3-ndigits,0])
	return ("%."+str(nafter)+"f") % num

def calc_stats_for_dic(dic):
	meandic={}
	stddic={}
	for key in dic:
		#entry in dic is a list of numbers
		meandic[key]=stformat(np.mean(dic[key]))
		stddic[key]=stformat(np.std(dic[key]))
	return meandic,stddic

def calc_fun_for_dic(dic,func):
	resdic={}
	for key in dic:
		resdic[key]=func(dic[key])
	return resdic

def calc_diffs(dic):
	

#just puts '-' where there's no value;)
def enlarge_dic(dic):
	keys=dic.keys()
	for (nsp,ngps,bl,i) in keys:
		if not (nsp,ngps,bl,(i+1)%3) in keys:
			dic[(nsp,ngps,bl,(i+1)%3)]='-'
		if not (nsp,ngps,bl,(i+2)%3) in keys:
			dic[(nsp,ngps,bl,(i+2)%3)]='-'


def write_row_species(meantdic,stdtdic,nspecies):
	ss=str(nspecies)+"	&	"+meantdic[(nspecies,1,0.2,1)]+" $\pm$ "+stdtdic[(nspecies,1,0.2,1)]+"	&	"
	ss+=meantdic[(nspecies,1,0.2,2)]+" $\pm$ "+stdtdic[(nspecies,1,0.2,2)]+"	&	"
	ss+=meantdic[(nspecies,1,0.2,0)]+" $\pm$ "+stdtdic[(nspecies,1,0.2,0)]+"	&	"
	# 3genes/sp
	ss+=meantdic[(nspecies,3,0.2,1)]+" $\pm$ "+stdtdic[(nspecies,3,0.2,1)]+"	&	"
	ss+=meantdic[(nspecies,3,0.2,2)]+" $\pm$ "+stdtdic[(nspecies,3,0.2,2)]+"	&	"
	ss+=meantdic[(nspecies,3,0.2,0)]+" $\pm$ "+stdtdic[(nspecies,3,0.2,0)]+"	&	"
	# 5genes/sp
	ss+=meantdic[(nspecies,5,0.2,1)]+" $\pm$ "+stdtdic[(nspecies,5,0.2,1)]+"	&	"
	ss+=meantdic[(nspecies,5,0.2,2)]+" $\pm$ "+stdtdic[(nspecies,5,0.2,2)]+"	&	"
	ss+=meantdic[(nspecies,5,0.2,0)]+" $\pm$ "+stdtdic[(nspecies,5,0.2,0)]+"	\\\\"
	ss2=ss.replace("- $\pm$ -","NA")

	return ss2

def write_row_species2(meantdic,stdtdic,nspecies,ngenes):
	ss=str(nspecies)+" species, "+str(ngenes)+ " gene/sp.	&	"+meantdic[(nspecies,ngenes,0.2,1)]+" $\pm$ "+stdtdic[(nspecies,ngenes,0.2,1)]+"	&	"
	ss+=meantdic[(nspecies,ngenes,0.2,2)]+" $\pm$ "+stdtdic[(nspecies,ngenes,0.2,2)]+"	&	"
	ss+=meantdic[(nspecies,ngenes,0.2,0)]+" $\pm$ "+stdtdic[(nspecies,ngenes,0.2,0)]+"	\\\\"

	ss2=ss.replace("- $\pm$ -","NA")

	return ss2

def write_row_genes(meantdic,stdtdic,ngenes):
	ss=str(ngenes)+"	&	"+meantdic[(2,ngenes,0.2,1)]+" $\pm$ "+stdtdic[(2,ngenes,0.2,1)]+"	&	"
	ss+=meantdic[(2,ngenes,0.2,2)]+" $\pm$ "+stdtdic[(2,ngenes,0.2,2)]+"	&	"
	ss+=meantdic[(2,ngenes,0.2,0)]+" $\pm$ "+stdtdic[(2,ngenes,0.2,0)]+"	&	"

	ss+=meantdic[(3,ngenes,0.2,1)]+" $\pm$ "+stdtdic[(3,ngenes,0.2,1)]+"	&	"
	ss+=meantdic[(3,ngenes,0.2,2)]+" $\pm$ "+stdtdic[(3,ngenes,0.2,2)]+"	&	"
	ss+=meantdic[(3,ngenes,0.2,0)]+" $\pm$ "+stdtdic[(3,ngenes,0.2,0)]+"	&	"
	# 5genes/sp
	ss+=meantdic[(4,ngenes,0.2,1)]+" $\pm$ "+stdtdic[(4,ngenes,0.2,1)]+"	&	"
	ss+=meantdic[(4,ngenes,0.2,2)]+" $\pm$ "+stdtdic[(4,ngenes,0.2,2)]+"	&	"
	ss+=meantdic[(4,ngenes,0.2,0)]+" $\pm$ "+stdtdic[(4,ngenes,0.2,0)]+"	\\\\"
	ss2=ss.replace("- $\pm$ -","NA")
	return ss2

def display_latex_table_species(meantdic,stdtdic):


	enlarge_dic(meantdic)
	enlarge_dic(stdtdic)

	print "\\begin{sidewaystable}"
	print "\\centering"
	print "\\begin{tabular}{|c||c|c|c||c|c|c||c|c|c|}"
	print "\# species & \\multicolumn{3}{||c||}{1 gene/species} & \\multicolumn{3}{||c||}{3 genes/species} & \\multicolumn{3}{||c||}{5 genes/species} \\\\"
	print "\\hline"
	print " 	& STELLS & CompactCH & GTProb & STELLS & CompactCH & GTProb & STELLS & CompactCH & GTProb \\\\"
	print "\\hline"
	print write_row_species(meantdic,stdtdic,4)
	print write_row_species(meantdic,stdtdic,8)
	print write_row_species(meantdic,stdtdic,16)
	print write_row_species(meantdic,stdtdic,24)
	print write_row_species(meantdic,stdtdic,32)
	print write_row_species(meantdic,stdtdic,40)


	print "\\end{tabular}"
	print "\\caption{The running times of the three exact algorithms for varying numbers of species.}"
	print "\\end{sidewaystable}"


def display_latex_table_genes(meantdic,stdtdic):


	enlarge_dic(meantdic)
	enlarge_dic(stdtdic)

	print "\\begin{sidewaystable}"
	print "\\centering"
	print "\\begin{tabular}{|c||c|c|c||c|c|c||c|c|c|}"
	print "\# samples & \\multicolumn{3}{||c||}{2 species} & \\multicolumn{3}{||c||}{3 species} & \\multicolumn{3}{||c||}{4 species} \\\\"
	print "\\hline"
	print " 	& STELLS & CompactCH & GTProb & STELLS & CompactCH & GTProb & STELLS & CompactCH & GTProb \\\\"
	print "\\hline"
	print write_row_genes(meantdic,stdtdic,5)
	print write_row_genes(meantdic,stdtdic,10)
	print write_row_genes(meantdic,stdtdic,15)
	print write_row_genes(meantdic,stdtdic,20)
	print write_row_genes(meantdic,stdtdic,50)


	print "\\end{tabular}"
	print "\\caption{The running times of the three exact algorithms for varying numbers of samples per species.}"
	print "\\end{sidewaystable}"

def display_latex_table_species2(meantdic,stdtdic):


	enlarge_dic(meantdic)
	enlarge_dic(stdtdic)

	print "\\begin{table}"
	print "\\centering"
	print "\\begin{tabular}{|c||c|c|c||}"
	#print "Number of species & "
	#print "\\hline"
	print write_row_species2(meantdic,stdtdic,4,1)
	print write_row_species2(meantdic,stdtdic,4,3)
	print write_row_species2(meantdic,stdtdic,4,5)

	print write_row_species2(meantdic,stdtdic,8,1)
	print write_row_species2(meantdic,stdtdic,8,3)
	print write_row_species2(meantdic,stdtdic,8,5)


	print write_row_species2(meantdic,stdtdic,16,1)
	print write_row_species2(meantdic,stdtdic,16,3)
	print write_row_species2(meantdic,stdtdic,16,5)


	print write_row_species2(meantdic,stdtdic,24,1)
	print write_row_species2(meantdic,stdtdic,24,3)
	print write_row_species2(meantdic,stdtdic,24,5)


	print write_row_species2(meantdic,stdtdic,32,1)
	print write_row_species2(meantdic,stdtdic,32,3)
	print write_row_species2(meantdic,stdtdic,32,5)

	print write_row_species2(meantdic,stdtdic,40,1)
	print write_row_species2(meantdic,stdtdic,40,3)
	print write_row_species2(meantdic,stdtdic,40,5)



	print "\\end{tabular}"
	print "\\caption{Running times}"
	print "\\end{table}"









ldic={}
tdic={}

Results=file('results','w')
timefiles=subprocess.check_output('ls|grep time',shell=True).split('\n')[:-1]
for f in timefiles:
	if 'timech' in f:
		outfile=f[:-6]+'outch'
	else:
		outfile=f[:-4]+'out'
	OF=file(outfile,'r')
	lines=OF.readlines()
	OF.close()
	if 'nonums' in f:
		finished=False
		lprob=100000000.0
		for line in lines:
			if 'Log-prob' in line:
				lpstr=line.split()[-1]
				lprob=float(lpstr)
				finished=True
		if finished:
			t=get_time(f)
			process_result(ldic,tdic,f,lprob,t)
	else:
		if len(lines)==0 or len(lines[0])<=2:
			continue
		finished=True
		lprob=float(lines[0])
		t=get_time(f)
		process_result(ldic,tdic,f,lprob,t)
					
meanldic,stdldic=calc_stats_for_dic(ldic)
meantdic,stdtdic=calc_stats_for_dic(tdic)

#display_latex_table_species(meantdic,stdtdic)

display_latex_table_genes(meantdic,stdtdic)

#print meantdic

