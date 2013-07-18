#!/usr/local/env python
"""
MODULE

Collection of wrapper functions for running DCON, IPEC, and PENT.

""" 
"""
	@package pypec
	@author NC Logan
	@email nlogan@pppl.gov
"""

#basics
import os                               # operating system commands
import sys
import subprocess
import time
import pickle
import numpy as np                      # math
from scipy import optimize
from string import join                 # string manipulation

#in this package
import data                             # read/write .out files
import _defaults_
import namelist                         # read/write fortran namelist .in files
from _bashjob_ import bashjob            # bash script skeleton

if os.getenv('HOST').startswith('sunfire'):
	print 'WARNING: USING PORTAL - "use -w 24:00:00 -l mem=8gb" recomended for heavy computation'

# make sure package is in path
packagedir = join(os.path.abspath(__file__).split('/')[:-2],'/')+'/' # one directory up
sys.path.insert(0,packagedir+'pypec')
sys.path.insert(0,packagedir+'rundir/'+_defaults_.rundir)


##################################################### DEFAULTS

class InputSet:
	"""
	Class designed to hold the standard inputs for a run.
	attr  - rundir : str. location of executables, .dat files, and extra .in files.
	        indir  : str. location of .in files read in and stored here
	        infiles: dict. tree like storage of namelist objects.
	"""
	def __init__(self,rundir=packagedir+'rundir/'+_defaults_.rundir,indir=packagedir+'input'):
		"""
		Initializes the intance.
		kwargs - rundir : str. location of IPEC package rundirectory.
		         indir  : str. location of .in namelist files to read.
		"""
		self.rundir = rundir
		self.indir = indir.rstrip('/')+'/' # want ending /

		self.infiles ={}
		for f in os.listdir(self.indir):
			if f.endswith('.in'):
				try:
					self.infiles[f[:-3]]=namelist.read(self.indir+f)
				except:
					print('Failed to read '+f)
					pass
			

default = InputSet()


##################################################### WRAPPER

def _newloc(loc):
	"""
	Utility function for making and moving to a new location.
	Note: First 2 directories counted pre-/ must exist already.
	
	args - loc : str. Directory to make and/or move to.
	"""
	dirs = os.path.abspath(loc).split('/')
	for i in range(2,len(dirs)+1):
		tmp = join(dirs[0:i],'/')
		if not os.path.exists(tmp):
			try:
				os.system('mkdir '+tmp)
			except:
				print('Could not make directory '+tmp)
				raise
	os.chdir(loc)
	return
		


def run(loc='.',rundir=default.rundir,qsub=True,return_on_complete=False,rundcon=True,runipec=True,runpent=True,mailon=_defaults_.mailon,email=_defaults_.email,**kwargs):
	"""
	Python wrapper for running ipec package.
	
	kwargs  - loc      : str. Directory location for run.
		  rundir   : str. IPEC package directory with executables and *.dat's.
	          qsub     : bool. Submit job to cluster.
		  return_on_complete : bool. Return only after job is finished on cluster (irrelevant if qsub=False).
		  rundcon  : bool. Run dcon.
		  runipec  : bool. Run ipec.
		  runpent  : bool. Run pent.
		  mailon   : str. Any combination of 'a','b','e' for on interruption
		                  execution, and termination respectively.
		  email    : str. Email address.
		  **kwargs : dict. "namelist" instance(s) written to <kwarg>.in file(s).
		            -NOTE: Namelists taken from (in order of priority): 
			           1) **kwargs 
				   2) *.in files from loc
				   3) *.in files from rundir

	returns -          : bool. True
	"""
	# housekeeping
	loc = os.path.abspath(loc)
	rundir = os.path.abspath(rundir)
       
	print('Running IPEC in '+loc)
	_newloc(loc)
	if rundcon:
		os.system('rm *.out')
		os.system('rm *.dat')
		os.system('rm *.bin')
	elif runipec:
		os.system('rm ipec_*')
		os.system('rm pent_*')
		os.system('rm ipdiag_*')
	elif runpent:
		os.system('rm pent_*')

	# set up run
	for key in kwargs:
		print('Writting '+key+'.in from keyword argument.')
		namelist.write(kwargs[key],key+'.in') #over-write if exists
	
	print('Using base package '+rundir)
	if rundir != loc:
		os.system('cp '+packagedir+'coil/*.dat .')
		os.system('cp '+packagedir+'pent/*.dat .')
		localins=[f for f in os.listdir('.') if f[-3:]=='.in' and f[:4]!='draw']
		rundirins=[f for f in os.listdir(rundir) if f[-3:]=='.in' and f[:4]!='draw']
		for infile in [f for f in rundirins if f not in localins]:
			if rundcon==False and runipec==False and f!='pent.in':
				continue # allows clean torque scan loops
			print('Copying '+infile+' from rundir.')
			os.system('cp '+rundir+'/'+infile+' .')

	# actual run
	if qsub:
		os.system('rm *.sh')
		jobname=os.path.abspath(loc).split('/')[-1]
		os.system('rm ../'+jobname+'.o*')
		exelist=''
		if rundcon: exelist+=rundir+'/dcon \n'
		if runipec: exelist+=rundir+'/ipec \n'
		if runpent: exelist+=rundir+'/pent \n'
		jobstr = bashjob.replace('jobnamehere',jobname)
		if mailon: jobstr = jobstr.replace('# --- emailoptionhere','#PBS -m '+mailon)
		if email:  jobstr = jobstr.replace('# --- emailhere','#PBS -M '+email)
		jobstr = jobstr.replace('runlochere',loc)
		jobstr = jobstr.replace('exelisthere',exelist)
		with open(jobname+'.sh','w') as f:
			f.write(jobstr)
		os.system('qsub '+jobname+'.sh')
		#wait for job completion to return
		if return_on_complete:
			while not os.path.exists('../'+jobname+'.oe'):
				time.sleep(1)	
			# clean up
			os.system('rm *.dat')
			os.system('rm draw*')
	else:
		if rundcon: os.system(rundir+'/dcon')
		if runipec: os.system(rundir+'/ipec')
		if runpent: os.system(rundir+'/pent')    
		# clean up
		os.system('rm *.dat')
		os.system('rm draw*')

	return True





##################################################### COMMON LOOP FUNCTIONS


def optntv(mlist,maxiter=50,loc='.',rundir=default.rundir,**kwargs):
	"""
	Python wrapper for running ipec package.
	args    - mlist    : list. Integer poloidal modes included on plasma surface.
	kwargs  - loc      : str. Directory location for run.
	         rundir   : str. IPEC package directory with executables and *.dat's.
		  **kwargs : dict. "namelist" instance(s) written to <kwarg>.in file(s).
		            -NOTE: Namelists taken from (in order of priority): 
			           1) **kwargs 
				   2) *.in files from loc
				   3) *.in files from rundir

	returns -          : bool. True
	"""
	mlist = map(int,mlist)
	mlist.sort()
	loc=os.path.abspath(loc)
		
	# Clean ipec.in of any error fields
	if 'ipec' not in kwargs:
		print('WARNING: No ipec.in specidied. Using copy from run directory')
		kwargs['ipec']=namelist.read(rundir+'/ipec.in',combine_arrays=0)
	kwargs['ipec']['IPEC_INPUT']['coil_flag']=False
	kwargs['ipec']['IPEC_INPUT']['data_flag']=False
	kwargs['ipec']['IPEC_INPUT']['harmonic_flag']=True
	kwargs['ipec']['IPEC_OUTPUT']['ntv_flag']=True
	for key in kwargs['ipec']['IPEC_INPUT'].keys():
		if key.startswith('cos') or key.startswith('sin'):
			del kwargs['ipec']['IPEC_INPUT'][key]

	# Include each m errorfield namelist variable in ipec input
	for m in mlist:
		for cs in ['cosmn','sinmn']:
			kwargs['ipec']['IPEC_INPUT'][cs+'('+str(m)+')'] = 1e-4

	# Set up base for run
	run(loc=loc,rundir=rundir,qsub=False,
		rundcon=True,runipec=False,runpent=False,**kwargs)
	inputs = InputSet(indir='.').infiles
	nn = inputs['dcon']['DCON_CONTROL']['nn']

	def totntv(cs,mlist=mlist):
		"""
		Returns total ntv torque given a list of the cosine coefficients 
		and sine coefficients each ordered largest m to lowest.
		"""
		cs = np.array(cs)
		cs = cs.reshape(2,-1)
		for m,c,s in zip(mlist,cs[0,:],cs[1,:]):
			kwargs['ipec']['IPEC_INPUT']['cosmn('+str(m)+')'] = c
			kwargs['ipec']['IPEC_INPUT']['sinmn('+str(m)+')'] = s
		run(loc=loc,rundcon=False,qsub=False,**kwargs)
		ntv = data.read('pent_n'+str(nn)+'.out')
		tot = ntv[0].params['total(T_phi)']
		return tot

	norm = len(mlist)*(1e-4)**2.0
	def sumsqr(x,norm=norm):
		return sum(np.array(x)**2.)-norm
	x0 = np.array([1e-4]*(len(mlist)*2))
	opt = optimize.fmin_slsqp(totntv,x0,iter=maxiter,eqcons=[sumsqr],full_output=1)
	with open('optimization.pkl','wb') as f:
		pickle.dump(opt,f)
	return opt
	


def wescan(base='.',scale=np.linspace(-2,2,20),**kwargs):
	"""
	Run pent for series of scaled omegaE values. 
	The function then waits  for jobs to complete and cleans the
	copied .bin files from the subdirectories of the scan.

	kwargs - base : str. Top level directory containing dcon and ipec runs.
	                -- Must contain euler.bin and ipec_order1 files.
	         scale: ndarray. The scale factors iterated over.
	returns -     : True.
	"""
	base = os.path.abspath(base)+'/'
	pent = namelist.read(base+'pent.in')
	for s in scale:
		name = 'we{0:.5}'.format(s)
		_newloc(base+name)
		os.system('cp '+base+'euler.bin .')
		os.system('cp '+base+'vacuum.bin .')
		os.system('cp '+base+'ipec_*.bin .')
		pent['PENT_INPUT']['ipd'] = s
		run(rundcon=False,runipec=False,mailon='',pent=pent)

	#while not os.path.isfile('../'+name+'.oe'):
	#	time.sleep(5)
	#for s in scale:
	#	name = 'we{0:.5}'.format(s)
	#	os.system('rm '+base+name+'/*.bin')
		
	return True



