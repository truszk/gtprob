import sys
import math
from math import log,exp
import os
import subprocess

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
			print_result(Results,f,lprob,t)
	else:
		if len(lines)==0 or len(lines[0])<=2:
			continue
		finished=True
		lprob=float(lines[0])
		t=get_time(f)
		print_result(Results,f,lprob,t)
					
Results.close()


