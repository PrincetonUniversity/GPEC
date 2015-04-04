#!/usr/local/env python
"""
:mod:`pypec.optio` -- Interfacing Tools for IPECOPT
===================================================

The optio module is a specialized tool for manipulating data 
from IPEC standard ascii output files into IPECOPT standard input files.

Here, we show some examples using common workflows.


Examples
---------

First, the user is expected to run DCON and IPEC for the equilibrium of interest.
The applied surface perturbation is of no consequence in this run, as it will
be varied according to the IPECOPT optimizer. The important thing is to have
set the singcoup_flag true. This produces the ipec_singcoup_svd_n#.out output file.

Assuming nn=1 in the dcon.in namelist:

>>> svd = read('examples/example_singcoup_svd_n1.out')
Casting table 1 into Data object.
Casting table 2 into Data object.
Casting table 3 into Data object.
Casting table 4 into Data object.
Casting table 5 into Data object.
Casting table 6 into Data object.
Casting table 7 into Data object.
Casting table 8 into Data object.
Casting table 9 into Data object.
Casting table 10 into Data object.
Casting table 11 into Data object.
Casting table 12 into Data object.


.. note:: 
  These examples can be tested by developers using ipython 
  in this directory as follows:

    In [1]: import data,doctest
   
    In [2]: doctest.testmod(data,verbose=True)

"""
"""
	@package pypec
	@author NC Logan
    @email nlogan@pppl.gov
"""

import numpy as np                               # math

#in this package
import data

######################################################## Global Variables

default_quiet = False

######################################################## IO FOR DATA OBJECTs

def svd_transform(ipec_file,opt_file,phi=False,nperp=np.inf,
                  mlow=-np.inf,mhigh=np.inf,quiet=default_quiet):
    """
    Transform IPEC singcoup output file to form usable as IPECOPT
    input.

    **Arguments:**
        ipec_file : str.
            IPEC output file.
        opt_file : str.
            File writen for use in IPECOPT.

    **Key Word Arguments:**
        phi : bool.
            Write area normalized Phi spectrum (default is b normal).
        quiet : bool.
            Controls terminal display.

    **Returns:**
        obj.
            Data object containing all svd eigenvectors.

    """
    # read data
    allsvd = data.read(ipec_file,quiet=quiet)
    svd = allsvd[0]
    nperp = int(min(nperp,svd.params['msing']))
    mpert = int(svd.params['mpert'])
    
    # consolidate into one data object
    keys,perpkeys = [],[]
    v = np.vstack((np.identity(mpert),
                   np.identity(mpert)*1j))
    if phi:
        prefix = 'Phi'
    else:
        prefix = 'b'
    for i in range(int(svd.params['msing'])):#enumerate(allsvd):
        key = prefix+str(i+1)
        svd.y[key] = allsvd[i].y[prefix]
        keys.append(key)
        if i<nperp: # form basis perp to svd vector
            x = svd.y[key]
            xnorm = x.dot(x.conj()).real
            v = v-np.abs([v.dot(x.conj())]).T.dot([x])/xnorm
    # collect M-R new basis functions
    for i,vperp in enumerate(v[:-nperp]):
        key = 'nr'+prefix+str(i+1) # non-resonant
        #print(np.abs(vperp.dot(vperp.conj())),np.abs(vperp.dot(v[-1,:].conj())))
        svd.y[key] = vperp
        perpkeys.append(key)

    # Write in consolidated file
    header = '{:>16} {:>16}\n{:16} {:16}'.format('num_vec','num_m',int(svd.params['msing']),mpert)
    tmp = svd.params
    svd.params = {}
    data.write(svd,opt_file,ynames=keys,header=header)
    if nperp>0:
        perp_file = opt_file.replace('.','_perp{}.'.format(nperp))
        header = '{:>16} {:>16}\n{:16} {:16}'.format('num_vec','num_m',2*mpert-nperp,mpert)
        data.write(svd,perp_file,ynames=perpkeys,header=header)
    svd.params = tmp
    
    return svd

def eigenvectors(opt_file,ms=range(-30,30),ortho=None,phi=False,**kwargs):
    """
    Write simple IPECOPT eigenvector file corresponding to 
    unit m vectors. If ortho vector supplied, the eigenvectors
    are perpendicular components using 
    v_prp = v_m-|v_m dot ortho_m| ortho_m/|ortho_m|.

    """
    # setup
    ms = np.array(ms)
    if phi:
        prefix = 'Phi'
    else:
        prefix = 'b'
    if 'fmt' not in kwargs: kwargs['fmt'] = '%16.8E'
    if 'delimiter' not in kwargs: kwargs['delimiter']=' '
    if 'comments' not in kwargs: kwargs['comments']=''

    # intelegent ortho treatment - get array from data obj
    if type(ortho)==data.DataBase:
        tmp = []
        mtmp = ortho.x[0].tolist()
        for m in ms:
            if m in ortho.x[0]:
                tmp.append(np.real(ortho.y[prefix][mtmp.index(m)]))
                tmp.append(np.imag(ortho.y[prefix][mtmp.index(m)]))
            else:
                tmp.append(0)
                tmp.append(0)
        ortho = np.array(tmp)/np.linalg.norm(tmp)
    # normalize given ortho vector
    elif ortho!=None:
        ortho = np.array(ortho)/np.linalg.norm(ortho)

    # compute eigenvectors orthogonal to ortho
    v = np.identity(2*len(ms))
    if ortho!=None:
        v-= v.dot(ortho)*ortho.reshape(-1,1)
        #v = v[:-2,:] # reduced dimensionality

    cv = []
    for i,eig in enumerate(v):
        cv.append(eig[i%2::2])
    cv = np.array(cv)
    
    table = np.vstack([ms.repeat(2),v])
    table = np.vstack([ms,cv])
    header = '{:>16} {:>16}\n{:16} {:16}\n\n'.format('num_vec','num_m',cv.shape[0]/2,cv.shape[1])
    header+= '{:>16}'.format('m')
    for i in range(1,v.shape[0]/2+1): 
        header+= ' {:>16} {:>16}'.format('real('+prefix+str(i)+')','imag('+prefix+str(i)+')')
    kwargs['header'] = header
    np.savetxt(opt_file,table.T,**kwargs)

    return table




