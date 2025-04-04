#!/usr/local/env python
"""
:mod:`pypec.gpec` -- Python Wrappers for FORTRAN Codes
======================================================

This is a collection of wrapper functions for running DCON, GPEC, and PENT.

If run from the command line, this module will call an associated GUI.

This module is for more advanced operators who want to control/run 
the fortran package for DCON, GPEC, and PENT. To get started, lets
explore the package

>>> gpec.
gpec.InputSet       gpec.join           gpec.optimize       gpec.run
gpec.bashjob        gpec.namelist       gpec.optntv         gpec.subprocess
gpec.data           gpec.np             gpec.os             gpec.sys
gpec.default        gpec.ntv_dominantm  gpec.packagedir     gpec.time

many of these are simply standard python modules imported by gpec
lets submit a full run in just 4 lines:
steal the settings from a previous run (leaving run director as deafult)

>>> new_inputs = gpec.InputSet(indir='/p/gpec/users/nlogan/data/d3d/ipecout/145117/ico3_n3/')

change what we want to

>>> new_inputs['dcon']['DCON_CONTROL']['nn']=2
>>> new_dir = new_inputs.indir[:-2]+'2'

*Double check that all inputs are what you expect* (not shown here)
and submit the job to the cluster.

>>> gpec.run(loc=new_dir,**new_inputs.infiles)

Running the package using this module will 

1. Create the new location (if not already existing), and cd to it
2. rm all .out, .dat, and .bin files already in the directory (not dcon ones if rundcon=False, not pent ones for runpent=False, etc).
3. cp all .dat and .in files from the run directoy to the new location
4. Write all namelists supplied in kwargs to file (overwriting)
5. Write a unique bash script for the run, and submit it OR run each fortran routine from the commandline.
6. Remove all .dat and xdraw files from the directory

You will not that there is significant oportunity for overwriting files,
and the user must be responsible for double checking that all inputs 
are correct and will not result in life crushing disaster for themselves
or a friend when that super special file gets over-written.

"""
"""
    @package pypec
    @author NC Logan
    @email nlogan@pppl.gov
"""

# basics
import os                               # operating system commands
import sys,copy
import time,pickle
import subprocess
import numpy as np                      # math
from scipy import optimize

#in this package
from . import data                             # read/write .out files
from . import namelist                         # read/write fortran namelist .in files
from ._bashjob_ import bashjob            # bash script skeleton
try:
    import gui
except ImportError: # not supported at GA
    print('Failed to import gui - must have enthought.traits')
    
if os.getenv('HOST',default='').startswith('sunfire'):
    print('WARNING: USING PORTAL - "use -w 24:00:00 -l mem=8gb" recomended for heavy computation')

# make sure package is in path
packagedir = '/'.join(os.path.abspath(__file__).split('/')[:-2])+'/' # one directory up
sys.path.insert(0,packagedir+'pypec')
sys.path.insert(0,packagedir+'bin')


np.norm = np.linalg.norm
iteration =1

##################################################### DEFAULTS

class InputSet:
    """
    Class designed to hold the standard inputs for a run.

    :Attributes:

        - rundir: str. Location of executables, .dat files, and extra .in files.
        - indir: str. Location of .in files read in and stored here
        - infiles: dict. Tree like storage of namelist objects.

    """
    def __init__(self,rundir=packagedir+'bin',indir=packagedir+'input'):
        """
        Initializes the instance.

        :param rundir: str. Location of GPEC package rundirectory.
        :param indir : str. Location of .in namelist files to read.

        """
        self.rundir = rundir.rstrip('/')+'/' # want ending /
        self.indir = indir.rstrip('/')+'/' # want ending /

        self.infiles ={}
        for f in os.listdir(self.indir):
            if f.endswith('.in') and not f.startswith('draw'):
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
    
    :param loc: str. Directory to make and/or move to.

    """
    dirs = os.path.abspath(loc).split('/')
    for i in range(2,len(dirs)+1):
        tmp = '/'.join(dirs[0:i])
        if not os.path.exists(tmp):
            try:
                os.system('mkdir '+tmp)
            except:
                print('Could not make directory '+tmp)
                raise
    os.chdir(loc)
    return


def run(loc='.', rundir=default.rundir, submit=True, return_on_complete=False, rerun=False,
        rundcon=True, rungpec=True, runpentrc=True, runrdcon=False, runstride=False, cleandcon=False, fill_inputs=False, version='1.3',
        mailon='NONE', email='', mem=3e3, hours=36, partition='general', runipec=False, qsub=None, **kwargs):
    """
    Python wrapper for running gpec package.
    
    :param loc: str. Directory location for run.
    :param rundir: str. GPEC package directory with executables and .dat files.
    :param submit: bool. Submit job to cluster.
    :param return_on_complete: bool. Return only after job is finished on cluster (irrelevant if qsub=False).
    :param rerun: bool. Does not delete existing .out files.
    :param rundcon: bool. Run dcon.
    :param runipec: bool. Run ipec.
    :param rungpec: bool. Run gpec.
    :param runpentrc: bool. Run pentrc.
    :param runrdcon: bool. Run rdcon and rmatch.
    :param runstride: bool. Run stride.
    :param cleandcon: bool. Remove euler.bin file after run is complete.
    :param fill_inputs: bool. Use inputs from rundir (see kwargs).
    :param version: str. Version of GPEC loaded (may set compilers if 1.2 or later)
    :param mailon: str. Choose from NONE, BEGIN, END, FAIL, REQUEUE, or ALL.
    :param email: str. Email address (default is submitting user).
    :param mem: floatMemory request of q-submission in megabytes (converted to integer).
    :param hours: int. Number of hours requested from job manager.
    :param partition: str. Specify a specific computing queue.
    :param kwargs: dict. namelist instance(s) written to <kwarg>.in file(s).

    .. note::
      Namelists taken from (in order of priority): 
      1) Key word arguments 
      2) .in files from loc
      3) .in files from rundir

    :returns: bool. True.

    """
    # housekeeping
    loc = os.path.abspath(loc)
    rundir = os.path.abspath(rundir)       
    print('Running Purturbed Equilibrium Package in '+loc)
    _newloc(loc)
    locfiles = os.listdir('.')
    if runipec:
        print("WARNING: runipec is deprecated in GPEC 1.0. Use rungpec.")
    if qsub:
        print("WARNING: qsub is being deprecated in GPEC 1.0. Use submit in the future.")
        submit = qsub

    if not rerun:
        print(' > Cleaning old run out, dat, and bin files')
        if rundcon:
            if np.any([f.endswith('.out') for f in locfiles]):
                os.system('rm *.out')
            if np.any([f.endswith('.dat') for f in locfiles]):
                os.system('rm *.dat')
            if np.any([f.endswith('.bin') for f in locfiles]):
                os.system('rm *.bin')
        elif rungpec or runipec:
            if np.any([f.startswith('gpec_') for f in locfiles]):
                os.system('rm gpec_*')
            if np.any([f.startswith('ipec_') for f in locfiles]):
                os.system('rm ipec_*')
            if np.any([f.startswith('gpec_diagnostics_') for f in locfiles]):
                os.system('rm gpdiag_*')
            if np.any([f.startswith('pent_') for f in locfiles]):
                os.system('rm pent_*')
        elif runpentrc:
            if np.any([f.startswith('pentrc_') for f in locfiles]):
                os.system('rm pentrc_*')

    # set up run
    for key in kwargs:
        print(' > Writting '+key+'.in from keyword argument.')
        namelist.write(kwargs[key],key+'.in') #over-write if exists
    
    print(' > Using package from '+rundir)
    if rundir != loc:
        # fill missing input files 
        if fill_inputs:
            localins=[f for f in os.listdir('.') if f[-3:]=='.in' and f[:4]!='draw']
            rundirins=[f for f in os.listdir(rundir) if f[-3:]=='.in' and f[:4]!='draw']
            missedins=[f for f in rundirins if f not in localins]
            for infile in missedins:
                print(' > Copying '+infile+' from rundir.')
                os.system('cp '+rundir+'/'+infile+' .')
        # path to package data files
        if os.path.isfile('coil.in'):
            coil = namelist.read('coil.in')
            if 'data_dir' not in coil['COIL_CONTROL']:
                coil['COIL_CONTROL']['data_dir'] = packagedir+'coil'
        if os.path.isfile('pentrc.in'):
            pentrc = namelist.read('pentrc.in')
            if 'data_dir' not in pentrc['PENT_INPUT']:
                pentrc['PENT_INPUT']['data_dir'] = packagedir+'pentrc'
        # request enough threads for openmpi
        ntasks = 1
        if rundcon and os.path.isfile('dcon.in'):
            dcon = namelist.read('dcon.in')
            if 'dcon_kin_threads' in dcon.get('DCON_CONTROL', {}):
                ntasks = dcon['DCON_CONTROL']['dcon_kin_threads']
                print(' > Requesting {:} threads'.format(ntasks))

    # actual run
    if submit:
        # name the job based on the directory it is run in
        dloc = os.path.abspath(loc).split('/')[-1]
        jobname=(data.getshot(loc)+'_'+dloc).lstrip('_')
        # clean up old job submission if it exists
        if os.path.exists('../{:}.oe'.format(dloc)):
            os.system('rm ../{:}.oe'.format(dloc))
        if os.path.exists('{:}.sh'.format(jobname)):
            os.system('rm {:}.sh'.format(jobname))
        # fill in and write shell script
        exelist='module purge; module load intel openmpi netcdf acml/5.3.1/ifort64 git \n'
        if rundcon: exelist+=rundir+'/dcon \n'
        if runrdcon: exelist+=rundir+'/rdcon \n' + rundir+'/rmatch \n'
        if runipec: exelist+=rundir+'/ipec \n'
        if rungpec: exelist+=rundir+'/gpec \n'
        if runpentrc: exelist+=rundir+'/pentrc \n'
        if runstride: exelist+=rundir+'/stride \n'
        if cleandcon: exelist+='rm euler.bin \n'

        jobstr = bashjob.format(name=jobname, ntasks=ntasks, mem=str(int(mem)), days=0, hours=int(hours),
                                partition=partition, location= loc, exes=exelist,
                                mailtype=mailon.upper(), mailuser=email, version=version)
        jobstr = jobstr.replace('#SBATCH --mail-type=NONE', '') # robust to sbatch versions
        with open(jobname+'.sh','w') as f:
            f.write(jobstr)
        # submit job to cluster
        os.system('sbatch '+jobname+'.sh')
        # wait for job completion to return
        if return_on_complete:
            print(' > Waiting for job to complete...')
            while not os.path.exists('../'+jobname+'.oe'):
                time.sleep(60)    
            # clean up
            os.system('rm *.dat')
    else:
        if rundcon: os.system(rundir+'/dcon')
        if runipec: os.system(rundir+'/ipec')
        if rungpec: os.system(rundir+'/gpec')
        if runpentrc: os.system(rundir+'/pentrc')
        # clean up
        if cleandcon: os.system('rm euler.bin')
        os.system('rm *.dat')

    return True


##################################################### COMMON LOOP FUNCTIONS


def optntv(mlist,maxiter=50,loc='.',rundir=default.rundir,**kwargs):
    """
    Python wrapper for running gpec package multiple times in
    a search for the normalized spectrum that maximizes NTV.
    
    :param mlist: list. Integer poloidal modes included on plasma surface.
    :param maxiter: int.Maximum number of iterations.
    :param loc: str. Directory location for run.
    :param rundir: str. GPEC package directory with executables and .dat files.
    :param kwargs: dict. namelist instance(s) written to <kwarg>.in file(s).
        
     .. note::
       Namelists taken from (in order of priority):
       1) Key word arguments 
       2) .in files from loc
       3) .in files from rundir

    :returns: bool. True.

    """
    mlist = map(int,mlist)
    mlist.sort()
    loc=os.path.abspath(loc)
        
    # Clean gpec.in of any error fields
    if 'gpec' not in kwargs:
        print('WARNING: No gpec.in specidied. Using copy from run directory')
        kwargs['gpec']=namelist.read(rundir+'/gpec.in',combine_arrays=0)
    kwargs['gpec']['GPEC_INPUT']['coil_flag']=False
    kwargs['gpec']['GPEC_INPUT']['data_flag']=False
    kwargs['gpec']['GPEC_INPUT']['harmonic_flag']=True
    kwargs['gpec']['GPEC_OUTPUT']['ntv_flag']=True
    for key in kwargs['gpec']['GPEC_INPUT'].keys():
        if key.startswith('cos') or key.startswith('sin'):
            del kwargs['gpec']['GPEC_INPUT'][key]

    # Include each m errorfield namelist variable in gpec input
    for m in mlist:
        for cs in ['cosmn','sinmn']:
            kwargs['gpec']['GPEC_INPUT'][cs+'('+str(m)+')'] = 1e-4

    # Set up base for run
    run(loc=loc,rundir=rundir,qsub=False,
        rundcon=True,rungpec=False,runpent=False,**kwargs)
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
            kwargs['gpec']['GPEC_INPUT']['cosmn('+str(m)+')'] = c
            kwargs['gpec']['GPEC_INPUT']['sinmn('+str(m)+')'] = s
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

def optpentrc(ms=range(-5,25),ttype='tgar',tfac=-1,perp1=False,norm=1e-3,qsub=True,
              method='Powell',maxiter=1000,loc='.',rundir=default.rundir,**kwargs):
    """
    Python wrapper for running gpec packages multiple times in
    a search for the normalized spectrum that maximizes NTV.
    
    Approach is as follows:
    1) Run DCON and GPEC in base directory with singcoup_flag True
    2) Run GPEC twice for each m, creating new subdirectories cosmn* and sinmn*
    3) Optimize a functional wrapper for PENTRC using optimize.fmin_slsqp
    - Spectra constrained to have 1 Gm^2 norm
    - wrapper uses data package to linearly combine gpec_xclebsch outputs
    - Initial guess is the first singcoup_svd mode from GPEC
    
    :param ms: list. Integer poloidal modes included on plasma surface.
    :param ttype: str.PENTRC_OUTPUT flag (fgar,tgar,rlar,etc).
    :param tfac: int.Optimization minimizes tfac*|T_phi|. Negative values maximize applied NTV.
    :param perp1: bool.Optimizes in space perpendicular to 1st GPEC SVD mode.
    :param norm: float.Amplitude of applied area normalized spectrum Phi_x (Tesla meter^2).
    :param qsub: bool.Submits individual jobs to cluster. Use gpec.run for qsub overall optimizer.
    :param method: str.Method used in scipy.optimize.minimize.
    :param maxiter: int.Maximum number of iterations.
    :param loc: str. Directory location for run.
    :param rundir: str. GPEC package directory with executables and .dat files.
    :param kwargs: dict. namelist instance(s) written to <kwarg>.in file(s).
        
     .. note::
       Namelists taken from (in order of priority):
       1) Key word arguments 
       2) .in files from loc
       3) .in files from rundir

    :returns: bool. True.

    """
    global iteration
    iteration=1
    # setup
    ms = map(int,ms)
    ms.sort()
    loc=os.path.abspath(loc)+'/'
    rundir=os.path.abspath(rundir)+'/'
    
    # get gpec.in
    if 'gpec' not in kwargs:
        if os.path.isfile(loc+'gpec.in'):
            kwargs['gpec']=namelist.read(loc+'gpec.in')
        else:
            print('WARNING: No gpec.in specified. Using copy from run directory')
            kwargs['gpec']=namelist.read(rundir+'gpec.in')
    if 'coil_flag' in kwargs['gpec']['GPEC_INPUT']:
        if kwargs['gpec']['GPEC_INPUT']['coil_flag']:
            print('WARNING: Includes additional field from coils as background.')
    if 'data_flag' in kwargs['gpec']['GPEC_INPUT']:
        if kwargs['gpec']['GPEC_INPUT']['data_flag']:
            print('WARNING: Includes additional field from data file as background.')
    # assert necessary flags
    kwargs['gpec']['GPEC_INPUT']['jsurf_in']=1
    kwargs['gpec']['GPEC_INPUT']['harmonic_flag']=True
    kwargs['gpec']['GPEC_OUTPUT']['xclebsch_flag']=True
    kwargs['gpec']['GPEC_OUTPUT']['singcoup_flag']=True
    kwargs['gpec']['GPEC_OUTPUT']['resp_flag']=True
    # remove any previous spectra
    for key in kwargs['gpec']['GPEC_INPUT'].keys():
        if key.startswith('cos') or key.startswith('sin'):
            del kwargs['gpec']['GPEC_INPUT'][key]
    # include flat cosine and sine component for each m
    for m in ms:
        for cs in ['cosmn','sinmn']:
            kwargs['gpec']['GPEC_INPUT'][cs+'('+str(m)+')'] = norm/np.sqrt(2.0*len(ms))
    
    # base run (runs the only DCON run and runs GPEC to get singcoup svd modes)
    print('-'*40+' Running base case')
    run(loc=loc,rundir=rundir,qsub=False,
        rundcon=True,rungpec=True,runpent=False,runpentrc=False,**kwargs)
    
    # retrieve base variables
    kwargs['dcon'] = namelist.read(loc+'dcon.in')
    nn = kwargs['dcon']['DCON_CONTROL']['nn']
    mem = (10/max(1,5-nn))*1e3
    # retrieve gpec svd modes
    print('Reading GPEC svd dominant modes...')
    sc1,sc2 = data.read(loc+'gpec_singcoup_svd_n{}.out'.format(nn),quiet=True)[:2]
    msc = sc1.x[0]
    Phi1,Phi2 = [],[]
    for m in ms:
        if m in msc:
            Phi1.append(sc1.y['b'][msc==m][0])
            Phi2.append(sc2.y['b'][msc==m][0])
        else:
            Phi1.append(0)
            Phi2.append(0)
    Phi1,Phi2 = np.array([Phi1,Phi2])
    amp1,amp2 = np.sqrt(np.sum(np.abs(Phi1)**2)),np.sqrt(np.sum(np.abs(Phi2)**2))
    print('Renormalizing |Phi1| = {:}, to {:} Tesla meter^2'.format(amp1,norm))
    Phi1 =(Phi1/amp1)*norm#/np.sqrt(2.0*len(ms))
    print('Renormalizing |Phi2| = {:}, to {:} Tesla meter^2'.format(amp2,norm))
    Phi2 =(Phi2/amp2)*norm#/np.sqrt(2.0*len(ms))
    
    # assert sub-run variables
    kwargs['gpec']['GPEC_INPUT']['idconfile']=loc+'euler.bin'
    kwargs['gpec']['GPEC_INPUT']['ieqfile']=loc+'psi_in.bin'
    for k in kwargs['gpec']['GPEC_OUTPUT']: # maximum speed
        if 'flag' in k: kwargs['gpec']['GPEC_OUTPUT'][k]=False
    kwargs['gpec']['GPEC_OUTPUT']['xclebsch_flag']=True
    kwargs['pentrc'] = namelist.read(loc+'pentrc.in')
    kwargs['pentrc']['PENT_INPUT']['peq_file'] = 'gpec_xclebsch_n{}.out'.format(nn)
    kwargs['pentrc']['PENT_INPUT']['idconfile']= loc+'euler.bin'
    # run gpec for each spectral component
    print('-'*40+' Submitting jobs for each m')
    
    
    for m in ms:
        # null all other components
        for m2 in ms:
            kwargs['gpec']['GPEC_INPUT']['cosmn('+str(m2)+')'] = 0
            kwargs['gpec']['GPEC_INPUT']['sinmn('+str(m2)+')'] = 0
        # apply pure component
        kwargs['gpec']['GPEC_INPUT']['cosmn('+str(m)+')'] = norm
        if not os.path.isfile(loc+'cosmn{}.oe'.format(m)):
            run(loc=loc+'cosmn{}'.format(m),rundir=rundir,qsub=True,mem=mem,
                rundcon=False,rungpec=True,runpent=False,runpentrc=True,**kwargs)
        kwargs['gpec']['GPEC_INPUT']['cosmn('+str(m)+')'] = 0
        kwargs['gpec']['GPEC_INPUT']['sinmn('+str(m)+')'] = norm
        if not os.path.isfile(loc+'sinmn{}.oe'.format(m)):
            run(loc=loc+'sinmn{}'.format(m),rundir=rundir,qsub=True,mem=mem,
                rundcon=False,rungpec=True,runpent=False,runpentrc=True,**kwargs)
    
    
    # wait for runs to finish
    print('Waiting for m runs to complete...')
    gpecdone = False
    while not gpecdone:
        test = [os.path.isfile(loc+'cosmn{}.oe'.format(m)) for m in ms]
        test+= [os.path.isfile(loc+'sinmn{}.oe'.format(m)) for m in ms]
        if np.all(test):
            gpecdone=True
        else:
            time.sleep(60*3)
    
    # read all the displacements (this takes time)
    print('Reading xclebsch files...')
    data.default_quiet=True
    xms = data.readall(loc,filetype='gpec_xclebsch_n{:}.out'.format(nn),quiet=True)
    
    # assert pentrc iteration variables
    kwargs['pentrc'] = namelist.read(loc+'pentrc.in')
    kwargs['pentrc']['PENT_INPUT']['peq_file'] = loc+'gpec_xclebsch_mopt_n{}.out'.format(nn)
    kwargs['pentrc']['PENT_INPUT']['idconfile']= loc+'euler.bin'
    for k in kwargs['pentrc']['PENT_OUTPUT'].keys():
        if 'flag' in k: kwargs['pentrc']['PENT_OUTPUT'][k]=False
    kwargs['pentrc']['PENT_OUTPUT'][ttype+'_flag']=True
    
    # functional wrapper of PENTRC for optimization
    def wrapper(cs,ms=ms):
        """
        Returns total ntv torque given a list of the cosine coefficients 
        and sine coefficients each ordered largest m to lowest.

        """
        ccs = condition(cs,direction=-1)
        # form linear superpostion of gpec results
        cosmn,sinmn = np.array(ccs).reshape(2,-1)#*np.array([[1,-1]]).T
        xclebsch = 0
        for m,c,s in zip(ms,cosmn,sinmn):
            xclebsch = xclebsch+(c*xms['cosmn'+str(m)][0]+s*xms['sinmn'+str(m)][0])
        # write the summed displacement
        #  - note read cannot distinguish between psi' and psi and numbers the second
        #  - note extra precision needed for psi values
        data.write(xclebsch,fname=loc+'gpec_xclebsch_mopt_n{}.out'.format(nn),
                   ynames=['xi^psi','xi^psi_1','xi^alpha'],fmt='%24.16E') 
        
        # run pentrc in que, but wait for result
        run(loc=loc,rundcon=False,rungpec=False,runpent=False,runpentrc=True,fill_inputs=False,
            mem=mem,qsub=qsub,return_on_complete=True,pentrc=kwargs['pentrc'])
        ntv = data.read('pentrc_{:}_n{:}.out'.format(ttype,nn),quiet=True)
        metric = tfac*np.abs(ntv[0].y['intT_phi'][-1]) # gets minimized
        #print('metric = {:.9E}'.format(metric))
        callback(cs)
        return metric
    
    # log progress
    with open('pentrc_optimization.log','w') as f:
        header = '{:>16} {:>16} {:>16} {:>16}'.format('i','T_phi','eqcons','overlap')
        for m in ms:
            header+=' {:>16} {:>16}'.format('real(m{})'.format(m),'imag(m{})'.format(m))
        f.write(header+'\n')
    def callback(x):
        global iteration
        #print('callback, iteration = {}'.format(iteration))
        ntv = data.read('pentrc_{:}_n{:}.out'.format(ttype,nn),quiet=True)
        tphi= ntv[0].y['intT_phi'][-1]
        check = checknorm(x)
        xin = condition(x,direction=-1)
        ovlp = overlap(xin)
        with open('pentrc_optimization.log','a') as f:
            f.write(('{:16.8E} {:16.8E} {:16.8E} {:16.8E}'+' {:16.8E}'*len(x)
                     +'\n').format(iteration,tphi,check,ovlp,*xin))
        if iteration%10==0:
            keys = sorted(ntv[0].y.keys()) # make order consistent
            ntv[0].params['iteration'] = iteration
            data.write(ntv[0],ynames=keys,
                       fname='pentrc_{:}_n{:}_i{:}.log'.format(ttype,nn,iteration))
        iteration+=1
    
    # initial guesses
    if not perp1:
        x0 = np.array([Phi1.real,Phi1.imag]).ravel()
    else:
        x0 = np.array([Phi2.real,Phi2.imag]).ravel()
    
    # Conditioning and constraints
    def condition(vec,x0=x0,direction=1):
        """Scale m components to have equal weighting (best guess)"""
        if direction>=0:
            convec = vec/x0
        else:
            convec = vec*x0
        convec*= norm/np.norm(convec)
        return convec
    def checknorm(vec,norm=norm,debug=False):
        """Equality or inequality constraint on the norm of the vector""" 
        vnorm = np.norm(vec)#np.sqrt(np.sum(np.abs(vec).ravel()**2))
        if debug: print("  vnorm,norm = {:}, {:}".format(vnorm,norm))
        diff = (vnorm-norm)/norm
        return diff
    def checknorm_prime(vec,norm=norm):
        """Jacobian of checknorm"""
        vnorm = np.norm(vec)
        return (np.abs(vm)/vnorm)/norm
    def overlap(vec,Phi1=Phi1,debug=False):
        """Equality constrain that solution be perpendicular to GPEC dominant mode"""
        compvec = np.sum(np.array(vec).reshape([2,-1])*np.array([[1,1j]]).T,axis=0)
        dot = compvec.dot(Phi1.conj()).real#np.sum(compvec.real*Phi1.real+compvec.imag*Phi1.imag)
        ovlp= dot/(np.norm(compvec)*np.norm(Phi1))
        return ovlp
    
    # Linear combination estimate
    print('Reading torques for linear estimate...')
    tm = data.readall(loc,filetype='pentrc_{:}_n{:}.out'.format(ttype,nn),quiet=True)
    tvec = np.array([tm['cosmn'+str(m)][0].y['intT_phi'][-1]
                     +1j*tm['sinmn'+str(m)][0].y['intT_phi'][-1] for m in ms])
    Phiq = tvec**0.5 * norm/np.norm(tvec**0.5)
    Phil = tvec * norm/np.norm(tvec)
    darray = np.array([ms,Phi1.real,Phi1.imag,Phi2.real,Phi2.imag,Phil.real,Phil.imag,Phiq.real,Phiq.imag]).T
    header ='GPEC: Area normalized poloidal GPEC spectrum that optimizes the PENTRC torque\n'
    header+='Phi_L estimated by linear approximation of PENTRC torque\n'
    header+='Phi_Q estimated by quadratic approximation of PENTRC torque\n'
    header+='Phi_T estimated by nonlinear {:} optimization of PENTRC torque\n\n'.format(method)
    header+='{:>16} {:>16} {:>16} {:>16} {:>16} {:>16} {:>16} {:>16} {:>16}'.format('m',
            'real(Phi_1)','imag(Phi_1)','real(Phi_2)','imag(Phi_2)',
            'real(Phi_L)','imag(Phi_L)','real(Phi_Q)','imag(Phi_Q)')
    fname = loc+'pentrc_optimized_spectrum_n{:}.log'.format(nn)    
    print('Writing estimated spectrum:\n  '+fname)
    np.savetxt(fname,darray,header=header,fmt='%16.8E',delimiter=' ',comments='')
    
    # ACTUAL OPTIMIZATION
    print('-'*40+' Optimizing PENTRC')
    #opt = optimize.fmin_slsqp(wrapper,x0,iter=maxiter,full_output=1,disp=2,
    #                          eqcons=[checknorm]+perp1*[checkperp1])
    cons = [{'type':'eq',
             'fun':checknorm,
             'jac':checknorm_prime}]#lambda x: np.array([np.sqrt(np.sum(np.abs(x).ravel()**2))-norm])}]
    if perp1: cons+=[{'type':'ineq','fun':checkperp1}]
    #opt = optimize.minimize(wrapper,condition(x0),method='SLSQP',options={'disp':True},constraints=cons)
    opt = optimize.minimize(wrapper,condition(x0),method=method,options={'disp':True})#,callback=callback)
    data.default_quiet=False
    
    # save output
    print('Pickling fmin_slsqp output -> pentrc_optimization.pkl')
    with open(loc+'pentrc_optimization.pkl','wb') as f:
        pickle.dump(opt,f)
    rspec,ispec = opt.x.reshape(2,-1)#opt[0].reshape(2,-1)
    darray = np.array([ms,Phi1.real,Phi1.imag,Phi2.real,Phi2.imag,Phil.real,Phil.imag,rspec,ispec]).T
    header+='{:>16} {:>16}'.format('real(Phi_T)','imag(Phi_T)')
    print('Writing nonlinear optimized spectrum -> '+fname)
    np.savetxt(fname,darray,header=header,fmt='%16.8E',delimiter=' ',comments='')
        
    return opt
    


def omegascan(omega='wp',base='.',scale=(-2,-1,-0.5,0.5,1,2),pentrcfile='pentrc.in',**kwargs):
    """
    Run pent for series of scaled nu/omega_phi/omega_E/omega_D/gamma_damp values. 
    User must specify full path to all files named in pent.in input file
    located in the base directory.
    Note that the rotation profile is scaled by manipulating omega_E
    on each surface to obtain the scaled rotation.

    :param omega: str. Choose from nu, wp, we, wd, or ga
    :param base: str. Top level directory containing dcon and gpec runs. Must contain euler.bin file.
    :param scale: ndarray. The scale factors iterated over.
    :param pentfile: str. Original input file.
    :param rundcon: bool. Set true if using hybrid kinetic MHD DCON. Will also attempt GPEC + PENT.

    :returns: bool. True.

    """
    # setup
    base = os.path.abspath(base)+'/'
    pentrc = namelist.read(pentrcfile)
    scale = map(float,scale)
    for k in ['dcon','gpec','pent','pentrc']:
        if 'run'+k not in kwargs:
            kwargs['run'+k] = k=='pentrc'
        elif kwargs['run'+k]:
            kwargs[k] = namelist.read(base+k+'.in')
            if k=='dcon':
                kwargs['equil'] = namelist.read(base+'equil'+'.in')
                kwargs['vac'] = namelist.read(base+'vac'+'.in')
            if k=='gpec':
                kwargs['coil'] = namelist.read(base+'coil'+'.in')
    
    # use some initial gpec run if already available
    if not kwargs['rungpec'] or kwargs['rundcon']:
        for k in pentrc['PENT_INPUT']:
            if 'file' in k:
                pentrc['PENT_INPUT'][k] = os.path.abspath(pentrc['PENT_INPUT'][k])
    
    # submit the scaled runs
    for s in scale:
        name = omega+'{0:.5}'.format(s)
        _newloc(base+name)
        pentrc['PENT_CONTROL'][omega+'fac'] = s
        run(pentrc=pentrc,**kwargs)
        
    return True

def nustarscan(base='.',scale=(0.1,1.0,10.0),pentfile='pent.in',
           scalen=True,scalet=True,**kwargs):
    """
    Run pent for series of *approximately* scaled collisionality values.
    Default scaling holds P=NT constant (NTV~P), but this is broken using
    scalen or scalet (scalet=False, for example, scales only density).
    
    User may want to include a non-zero ximag if going to low collisionality.
    
    User must specify full path to all files named in the pent input.
    
    Kinetic file must have header with 'n' and 'i' in the ion density label
    but no other label, both 't' and 'e' in only the electron temperature label,
    'omega' omega_EXB label, etc.
    
    .. note:: Approximate scaling ignores log-Lambda dependence, and uses
       nu_star ~ nu/v_th ~ (NT^-3/2)/(T^1/2) ~ N/T^2.

    :param base: str. Top level directory containing dcon and gpec runs. Must contain euler.bin file.
    :param scale: ndarray. The scale factors iterated over.
    :param pentfile: str. Original input file.
    :param scalen: bool.Use density scaling to change nu_star.
    :param scalet: bool.Use temperature scaling to change nu_star.

    :returns: bool. True.

    """
    # setup
    base = os.path.abspath(base) +'/'
    pent = namelist.read(pentfile)
    kfile = pent['PENT_INPUT']['kinetic_file']
    kin, = data.read(pent['PENT_INPUT']['kinetic_file'])
    markers = [['n','i'],['n','e'],['t','i'],['t','e'],['omega']]
    ynames = []
    for ms in markers:
        for key in kin.y.keys():
            unitless = key.replace('eV','').replace('m3','').replace('rads','')
            if np.all([m in unitless for m in ms]): ynames.append(key)
    # make sure names arent poluted
    assert(len(ynames)==5)
    #ynames = ['nim3', 'nem3',  'tieV', 'teeV', 'wexbrads']
    #for name in ynames:
    #    assert(name in kin.y)
    if not scalen and not scalet:
        raise ValueError('Must scale density, temperature, or both.')

    #loop through the scan
    for s in map(float,scale):
        # distribute nu scale to N and T (keeping NT constant if allowed)
        nscale = scalen*s**(scalet*(1.0/3)+1.0*(not scalet)) + (not scalen)
        tscale = scalet*s**(scalen*(-1./3)-0.5*(not scalen)) + (not scalet)
        print('nuscale = {:.3e} \n tscale = {:.3e}, nscale = {:.3e}, 1/tscale = {:.3e}'.format(s,nscale,tscale,1.0/tscale))
        if(nscale and tscale): assert(np.isclose(nscale,1.0/tscale)) # double check NT constant scaling
        newkfile = ('.'.join(kfile.split('.')[:-1])
                +'_nustar{:.2e}.dat'.format(s))
        kcop = copy.deepcopy(kin)
        kcop.y[ynames[0]]*=nscale
        kcop.y[ynames[1]]*=nscale
        kcop.y[ynames[2]]*=tscale
        kcop.y[ynames[3]]*=tscale
        data.write(kcop,fname=newkfile,ynames=ynames)

        # run pent with new kinetic profile
        pent['PENT_INPUT']['kinetic_file']=os.path.abspath(newkfile)
        
        run(loc=base+'nustar{:.2e}'.format(s),rundcon=False,
            rungpec=False,pent=pent,**kwargs)
    
    return True

if __name__ == '__main__':
    gui.main()
