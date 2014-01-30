#!/usr/local/env python
"""
MODULE

Collection base class types for DCON, IPEC,
and PENT outputs.
Base classes define operators for the output
appropriate to the independent and dependent
data contained within.

"""
"""
	@package pypec
	@author NC Logan
    @email nlogan@pppl.gov
"""

import os,copy,time

from string import join,capitalize,count         # faster than '+' loop
from StringIO import StringIO
from types import MethodType

import numpy as np                               # math
from scipy.interpolate import interp1d,interp2d,LinearNDInterpolator,griddata
from scipy.special import ellipk,ellipe

#in this package
import modplot as plt
import namelist 

try:
	import enthought.mayavi.mlab as mlab
except ImportError:
	mlab = False
	print('WARNING: Mayavi not in python path')




######################################################## IO FOR DATA OBJECTs

def read(fname,squeeze=False,forcex=[],rdim=2,quiet=False):
	"""
	Get the data from any ipec output as a list of python 
	class-type objects using numpy.genfromtxt.

	args
	---
	- fname	  : str. Path of ipec output file.
	
	kwargs
	------
	- squeeze : bool. Sets all attributes to single data object.
	- forcex  : list. Set independent data labels.
	- rdim    : int. Dimensionality of data assumed when reducing 
	                 extrodinarily (>1e6) large tables of data.
	- quiet   : bool. Prevent non-warning messages printed to terminal.
	
	returns
 	-------
	          : list. Class objects for each block of data in the file.
	"""
	#debug start_time = time.time()

	collection = []
	dcollection = []
	pcollection = []
	
	with open(fname) as f:
		#Get preamble
		preamble = ''
		lastline = ''
		for lstart,line in enumerate(f):
			try:
				test = float(line.split()[0])
				break
			except:
				preamble+= lastline
				lastline = line
		f.seek(0)
		try: #Simple file : single table
			data = np.genfromtxt(f,skip_header=lstart-1,names=True)
			pcollection.append(preamble)
			dcollection.append(data)
		#look through file manually
		except ValueError as e:
			data = []
			f.seek(0)
			# clean up messy dcon output formats
			if 'dcon.out' in fname:
				lines = f.read()
				lines=lines.replace('mu0 p','mu0p')
				lines=lines.replace('re w','realw')
				lines=lines.replace('im w','imagw')
				lines=lines.replace('abs w','absw')
				lines=lines.replace('*','')
				lines = StringIO(lines).readlines()
			else:
				lines = f.readlines()
			top,bot = 0,0
			length= len(lines)
			while bot<length and top<(length-2):
				preamble=''
				lastline=''
				for i,line in enumerate(lines[bot:]):
					try:
						test = float(line.split()[0])
						top = bot+i
						if lines[top+1]=='\n': 
							raise ValueError #single lines are bad
						break
					except:
						preamble+= lastline
						lastline = line
				if line == lines[-1]: break #end of file without another table
				try:
					bot = top+lines[top:].index(' \n')
				except ValueError:
					try:
						bot = top+lines[top:].index('\n')
					except:
						bot = length
				# include headers
				top-=1
				if not lines[top].translate(None,' \n\t'): #empty line
					top-=1
				skipfoot = length-bot
				f.seek(0)
				table = lines[top:bot]
				if '\n' in table: #empty space
					table.remove('\n')
				data = np.genfromtxt(StringIO(join(table)),names=True,dtype=None)
				pcollection.append(preamble)
				dcollection.append(data)
		
		#debug print("Finished parsing file in "
		#debug 	  +"{} seconds".format(time.time()-start_time))
		#turn arrays into classes
		for i,(data,preamble) in enumerate(zip(dcollection,pcollection)): 
			if not quiet: 
				print("Casting table "+str(i+1)+" into Data object.")
			collection.append(DataBase(data,preamble,forcex=forcex))
			
	# force all attributes to single object
	if squeeze:
		dc1 = collection[0]
		if len(collection)>1:
			for dc in collection[1:]:
				dc1+=dc
		return dc1

	return collection


def readall(base='.',filetype='pent_n1.out',**kwargs):
	"""
	Recursively searches base directory for files with name
	filetype, reads them into python objects, and stores 
	objects in a dictionary with directory names as keys.
	
	kwargs
	------
	- base     : str. Top level directory in which to begin search.
	- filetype : str. Files to read.
	**kwargs          :dict. Passed to appropriate reader function.
    
	returns
	-------
	-          : dict.
	"""

	if('.in' in filetype):
		reader = namelist.read
	else:
		reader = read
	base = os.path.abspath(base)+'/'
	#if(base[-1]!='/'): base+='/'
		
	subs = [name for name in os.listdir(base) if os.path.isdir(base+name)]
	if subs:
		results = {}
		for subd in subs:
			results[subd] = readall(base=base+subd,filetype=filetype,**kwargs)
			if os.path.isfile(base+filetype):
				results[filetype] = reader(base+filetype,**kwargs)
		_cleandict(results)
	else:
		try:
			results = reader(base+filetype,**kwargs)
		except:
			results = None

	return results

def _cleandict(d):
	""" 
	Helper function to clear empty keys from dictionary.
	"""
	ks = d.keys()
	for k in ks:
		if type(d[k])==dict:
			_cleandict(d[k])
		if not d[k]:
			del d[k]
	#return d


def write(dataobj,fname='',ynames=[],**kwargs):
	"""
	Write data object to file. The aim is to replicate the original file
	structure... the header text will have been lost in the read however,
	as will some column label symbols in the use of genfromtxt.
	
	args
	----
	- dataobj : obj. Data_Base type object.
	- fname   : str. File to write to.
	
	kwargs
	------
	**kwargs  : dict. Passed to numpy.savetxt

	returns
	-------
	-         : bool. True.
	"""
	if not fname: return False
	
	if 'fmt' not in kwargs: kwargs['fmt'] = '%16.8e'
	if 'delimiter' not in kwargs: kwargs['delimiter']=' '
	if 'comments' not in kwargs: kwargs['comments']=''
	if 'header' in kwargs:
		kwargs['header']+= '\n\n'
	else:
		kwargs['header']=''

	for k,v in dataobj.params.iteritems():
		kwargs['header']+= k+" = "+str(v)+"\n"

	nums = []
	if not ynames: ynames = dataobj.y.keys()
	for i,name in enumerate(dataobj.xnames):
		nums.append(dataobj.pts[:,i])
		kwargs['header']+=kwargs['delimiter']*(i!=0)+'{:>16}'.format(name)
	for i,name in enumerate(ynames):
		nums.append(dataobj.y[name])
		kwargs['header']+=kwargs['delimiter']+'{:>16}'.format(name)
	np.savetxt(fname,np.array(nums).T,**kwargs)

	return True




def plotall(results,xfun,yfun,label='',axes=None,**kwargs):
	"""
	Line plot of data gethered recursively from a dictionary 
	containing multiple namelists or data objects.

	args
	----
	- results : dict. A readall result.
	- xfun    : func. A function that takes args (key,val)
	                  for each dictionary item and whose result
		              is appended to the xaxis.
	- yfun    : func. A function that takes args (key,val)
	                  for each dictionary item and whose result
					  is appended to the yaxis.

	kwargs
	------
	- axes    : obj. matplotlib axes object.
	**kwargs         : dict. Passed to matplotlib plot function
	
	returns
	-------
	-         : fig.
	"""
	x = []
	y = []
	if not axes:
		f,ax = plt.subplots()
	else:
		ax = axes

	def newpoint(key,val):
		if type(val) == dict:
			for k,v in val.iteritems():
				newpoint(k,v)
		else:
			try:
				x.append(xfun(key,val))
				y.append(yfun(key,val))
			except:
				pass
	for key,val in results.iteritems():
		newpoint(key,val)

	z = zip(x,y)
	z.sort()
	x,y = zip(*z)
	ax.plot(x,y,label=label,**kwargs)
	ax.legend()
	f = ax.get_figure()
	f.show()
	return f




def set_machine(machine):
	"""
	Reads built-in geometric instances for the machine.
	These consist of surface objects, representing vessel
	walls, divertors, magnetic sensors, or any other 2D
	surface in the machine.
	
	args
	----
	- machine : str. Choose from "d3d", ...
	
	returns
	-------
	-         : obj. Instance for visualizing surface objects.
	"""
	machinefile = '_'+machine+'_.dat'
	class newmachine(SurfaceBase):
		surfs = {}

	# read basics from file
	with open(machinefile,'r') as f:
		lines = f.readlines()
		length= len(lines)
		top,bot = 0,0
		#skip preamble
		while lines[bot][0]!='#':
			bot+=1
		while bot<length:
			try:
				top = lines.index('"""\n',bot+1)
			except:
				break #stop reading if no more labels
			name = lines[top+1].translate(None,'\t\n').rstrip(' ')
			desc = lines[top+2].translate(None,'\t\n').rstrip(' ')
			top = lines.index('"""\n',top+1)
			try:
				bot = lines.index('\n',top+1)
			except ValueError:
				bot=length #in case no empty line at end
			f.seek(0)
			table = np.genfromtxt(StringIO(join(lines[top+1:bot])),names=True,dtype=None)
			newmachine.surfs[name]={'description':desc}
			for key in table.dtype.names:
				newmachine.surfs[name][key]=table[key]

	#fill out 2D surface in 3D space (and convert to rads)
	dtor = np.pi/180
	for key,val in newmachine.surfs.iteritems():
		#conversions
		for subkey in ['angle','width','phi']:
			val[subkey] *= dtor
		#3D fill
		ncm = int(val['length']*100)
		val['rzphi'] = []
		for r,z,a,p,w,l in zip(val['r'],val['z'],val['angle'],val['phi'],val['width'],val['length']):
			nphi = max(1,int(w/dtor)) #point per degree
			rs=np.linspace(r-l/2*np.cos(a),r+l/2*np.cos(a),ncm).repeat(nphi).reshape(-1,nphi)
			zs=np.linspace(z-l/2*np.sin(a),z+l/2*np.sin(a),ncm).repeat(nphi).reshape(-1,nphi)
			ps=np.linspace(a-w/2,a+w/2,nphi).repeat(ncm,axis=0)
			val['rzphi'].append([rs,zs,ps])
			
	return newmachine

######################################################## THE BASE IPEC DATA OBJECTS

class DataBase:
	"""
	An emtpy class object with the following operations overloaded:
	+,-,*,/,**
	such that the operation is only performed on the attributes 
	.real_* and .imag_* where * can be r,z, or phi. If the two
	factors in the operation have different .r or .z attributes
	data is interpolated to the wider spaced axis.
	"""
	
	def __init__(self,fromtxt,preamble,forcex=[]):
		"""
		Takes structured array from numpy.genfromtxt, and creates
		dictionaries of dependent and independent data. Also breaks
		down preamble to find any global parameters.
		
		args
		----
		- fromtxt : structured array. 
		- preamble: str.
		- forcex  : list. x-axis labels.
		"""
		#debug start_time = time.time()
		names = list(fromtxt.dtype.names)
		l = len(fromtxt)
		potxnames = names[:3]
		# Set independent/dependent variables (max 3D,min 3pt grid)
		if forcex:
			nd = len(forcex)
			potxnames = forcex
		else:
			nd = 1
			# hack for psi,q pairs
			if 'q' in potxnames:
				potxnames = names[:4]
				potxnames.remove('q')
			# hack for rzphi
			if 'l' in potxnames:
				potxnames = names[:4]
				potxnames.remove('l')
			if len(set(fromtxt[potxnames[1]]))<l/3: 
				nd+=1
				if len(set(fromtxt[potxnames[2]]))<l/3:
					nd+=1
			# hack for pent
			if nd==1 and potxnames[1]=='Lambda':
				nd+=1
			if nd==2 and potxnames[2]=='x':
				nd+=1
		#debug print("Determined dimensionality in {} seconds".format(time.time()-start_time))
		#debug start_time = time.time()
		self.xnames= potxnames[:nd]	
		x = fromtxt[potxnames[:nd]]
		self.x = [np.sort(list(set(x[name]))) for name in self.xnames]
		self.shape = [len(ax) for ax in self.x]
		#debug print("Set x in {} seconds".format(time.time()-start_time))
		#debug print("shape = "+str(self.shape))
		#debug start_time = time.time()		
		# reduce exesively large files, assuming 1st column is x
		if len(fromtxt)>1e6:
			name = self.xnames[0]
			step = l/int(1e6)+1
			if np.product(self.shape) != l:
				fromtxt = fromtxt[::step] # no axes to muck up
			else:
				mask = np.ones(self.shape,dtype=bool)
				mask[::step,...] = False
				#debug print('set mask in {} seconds'.format(time.time()-start_time))
				fromtxt = fromtxt[mask.ravel()]
			#debug print('masked fromtxt in {} seconds'.format(time.time()-start_time))
			print("WARNING: Reducing length of "
				  +"table from {:} to {:}".format(l,len(fromtxt))
				  +" by reducing x[0] by a factor of {}".format(step))
			l = len(fromtxt)
			x = fromtxt[potxnames[:nd]]
			self.x = [np.sort(list(set(x[name]))) for name in self.xnames]
			self.shape = [len(ax) for ax in self.x]
			#debug print("Set x in {} seconds".format(time.time()-start_time))
			#debug start_time = time.time()
		self.pts  = x.view(np.float).reshape(l,-1)
		#debug print("Set pts in {} seconds".format(time.time()-start_time))
		#debug start_time = time.time()
		self.nd    = np.shape(self.pts)[-1]
		if np.product(self.shape) != l: #probably psi to closely packed
			# maybe first axis is packed too tight, try trusting others
			if self.nd ==1:
				self.shape = [l]
				self.x[0] = self.pts[:,0]
			if self.nd ==2:
				self.shape[0] = l/self.shape[1] 
				self.x[0]  = self.pts[:,0][::self.shape[1]]
			if self.nd ==3:
				self.shape[0] = l/np.product(self.shape[1:])
				self.x[0] = self.pts[:,0][::np.product(self.shape[1:])]
			# test result
			if np.any(list(set(self.pts[:,0]).difference(self.x[0]))):
				print("WARNING: Irregular dependent data: can not form axes for interpolation.")
				self.x=None
			else:
				print("Warning: Axes contain repeated values. Forced axes my be incorrect.")
		#debug print("Checked axes in {} seconds".format(time.time()-start_time))
		#debug start_time = time.time()

		self.y={}
		ynames = [name for name in names if name not in self.xnames]
		for name in ynames:
			if 'real' in name:
				newname = name.replace('real','')
				self.y[newname] = fromtxt[name]+1j*fromtxt[name.replace('real','imag')]
			elif not 'imag' in name:
				self.y[name]=fromtxt[name]
		#debug print("Set y in {} seconds".format(time.time()-start_time))
		#debug start_time = time.time()

		#scipy interpolation function for ND (stored so user can override)
		self._interpfunc=[interp1d,interp2d,LinearNDInterpolator][self.nd-1]
		self._interpdict={}

		# Set any variables from preamble
		preamble=preamble.split()
		params=[]
		while preamble.count('='):	#paramter specifications
			idx = preamble.index('=')
			param = preamble[idx+1]
			name = preamble[idx-1]
			name.translate(None,'()[]{}\/*-+^%.,:!@#&')
			try: 
				params.append((name,float(param))) #if its a number
			except: 
				params.append((name,param))        #if its not a number
			preamble.remove(preamble[idx-1])
			preamble.remove(str(param))
			preamble.remove('=')
		self.params=dict(params)
		self.__doc__ = join(preamble,'')	#whatever is left

		#debug print("Set preamble in {} seconds".format(time.time()-start_time))
		#debug start_time = time.time()


	# Operations overrides
	def __add__(self,other): return _data_op(self,other,np.add,'+')
	def __sub__(self,other): return _data_op(self,other,np.subtract,'-')
	def __mul__(self,other): return _data_op(self,other,np.multiply,'*')
	def __div__(self,other): return _data_op(self,other,np.divide,'/')
	def __pow__(self,other): return _data_op(self,other,pow,'^')
	def __radd__(self,other): return self.__add__(other)
	def __rmul__(self,other): return self.__mul__(other)
	def __rsub__(self,other): return self.__sub__(other)
	def __rdiv__(self,other): return self.__div__(other)
	

	# Built in interpolation
	def interp(self,x,ynames=None,**kwargs):
		"""
		Interpolate data to point(s).
		
		args
		----
		- x  : tuple. Point(s) on dependent axis(axes).
		
		returns
		-------
		-    : dict. Each name contains array of interpolated data.  
		"""
		#housekeeping
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=(ynames,)
		if not type(x) in (list,tuple): x=(x,)
		
		#function to form interpolator if needed
		def new_interp(name):
			args = []
			if self.nd==1 or len(self.pts)<1e3:
				for n in range(self.nd):
					args.append(self.pts[:,n])
				args.append(self.y[name])
			else:
				myslices=[]
				for n in range(self.nd):
					myslices.append(np.s_[::max(1,self.shape[n]/(150/self.nd))])
				args = [self.pts[:,n].reshape(self.shape)[myslices] for n in range(self.nd)
						]+[self.y[name].reshape(self.shape)[myslices]]
				newlen = np.product(np.shape(args[0]))
				print("Warning: "+str(newlen)+" of "+str(len(self.pts))+" points being passed to interpolator")
			if np.any(np.iscomplex(self.y[name])):
				print("Forming real interpolator for "+name+".")
				freal = self._interpfunc(*np.real(args),**kwargs)
				print("Forming imaginary interpolator (assumed independent).")
				args[-1] = np.imag(args[-1])
				fimag = self._interpfunc(*args,**kwargs)
				return lambda *xy: freal(*np.real(xy))+1j*fimag(*np.imag(xy))
			else:
				print("Forming interpolator for "+name+".")
				return self._interpfunc(*args,**kwargs)
				
		# for each name check if interpolator is up to date and get values
		atts = ['y','z','values']
		att = atts[self.nd-1]
		values={}
		for name in ynames:
			if name in self._interpdict.keys():
				if not all(getattr(self._interpdict[name],att) == self.y[name]):
					self._interpdict[name] = new_interp(name)
			else:
				self._interpdict[name] = new_interp(name)
			print("Interpolating values for "+name+".")
			values[name]=self._interpdict[name](*x)
		return values
		
    
	# Built-in visualizations
	def plot1d(self,ynames=None,xname=None,x1rng=None,x2rng=None,x3rng=None,**kwargs):
		"""
		Line plots of from 1D or 2D data.

		kwargs
		------
		- ynames : list. Strings specifying dependent data displayed.
		- xname  : string. The independent axis if 2D data.
		- x1rng  : . Valid formats are:
		          -- x            : Plot single nearest point x-axis.
				  -- (xmin,xmax)  : Plots all data in range (min,max). 
				  -- (xmin,xmax,n): Plots n evenly spaced points in range.
		- x2rng  : . Same as x1rng. Applicable to 2D or 3D data.
		- x3rng  : . Same as x2rng for 3D data only. Order of x is 
		             xnames with xname placed in front if given.
		**kwargs         : Valid matplotlib pyplot.plot keyword arguments.
		
		returns
		-------
		-        : figure. poloidal cross sections of plasma response
		"""
		# housekeeping
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=(ynames,)
		if not xname: xname=self.xnames[0]

		indx = self.xnames.index(xname)
		if self.x:
			x=self.x[indx]
		else:
			x = self.pts[:,indx]

		# helper to sort through axes of ND data
		def maskmaker(ax,rng):
			if type(rng) in [float,int]: 
				mask = ax==ax[abs(ax-rng).argmin()]
			elif type(rng) in [tuple,list] and len(rng)==2:
				mask = (ax>=rng[0]) * (ax<=rng[1])
			elif type(rng) in [tuple,list] and len(rng)>2:
				print("Warning: mapping all axes to regular grid.")
				if self.nd==2:
					x,ax = np.mgrid[x.min():x.max():x.size,rng[0]:rng[1]:rng[2]]
					x = x[0]
					mask = range(rng[2])
				if self.nd==3:
					raise ValueError("Griding for 3D data is under developement.")
			else:
				mask = np.ones_like(ax,dtype=bool)
			return ax,mask
		
        # display data
		if 'figure' in kwargs:
			f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
		elif 'axes' in kwargs:
			if type(kwargs['axes'])==np.ndarray:
				f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
			else:
				f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
		else:
			nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
			f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=(6+2*ncol,4+2*nrow))


		if 'label' in kwargs:
			ulbl = kwargs['label']+', '
			del kwargs['label']
		else:
			ulbl = ''
		
		for a,name in zip(ax.ravel(),ynames):
			lbl = ulbl + _mathtext(name)
			if self.nd==1:
				a.plot(x,self.y[name],label=lbl,**kwargs)
			elif self.nd==2:
				if self.x:
					x2,rng2 = maskmaker(self.x[indx-1],x2rng)
					ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
					xlbl = _mathtext(self.xnames[indx-1])
					for y,x2 in zip(ys[rng2],x2[rng2]):
						x,rng1 = maskmaker(x,x1rng)
						a.plot(x[rng1],y[rng1],label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
				else:
					x2s,rng2 = maskmaker(self.pts[:,indx-1],x2rng)
					xlbl = _mathtext(self.xnames[indx-1])
					x2s = np.sort(list(set(x2s[rng2])))
					for x2 in x2s:
						x,rng2 = maskmaker(self.pts[:,indx-1],float(x2))
						x = self.pts[:,indx][rng2]
						y = self.y[name][rng2]
						x,y = np.array(zip(*sorted(zip(x,y))))
						x,rng1 = maskmaker(x,x1rng)
						a.plot(x[rng1],y[rng1],label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
			elif self.nd==3:
				indx23 = range(3)
				indx23.remove(indx)
				if self.x:
					x2,rng2 = maskmaker(self.x[indx23[0]],x2rng)
					x3,rng3 = maskmaker(self.x[indx23[1]],x3rng)
					ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
					xlbls = [_mathtext(name) for name in self.xnames if name!=xname]
					for y2,x2 in zip(ys[rng2],x2[rng2]):
						x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
						for y,x3 in zip(y2[rng3],x3[rng3]):
							x,rng1 = maskmaker(x,x1rng)
							a.plot(x[rng1],y[rng1],label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
				else:
					x2s,rng2 = maskmaker(self.pts[:,indx23[0]],x2rng)
					x2s = np.sort(list(set(x2s[rng2])))
					xs = self.pts[:,indx][rng2]
					ys = self.y[name][rng2]
					xlbls = [_mathtext(xlbl) for xlbl in self.xnames if xlbl!=xname]
					for x2 in x2s:
						x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
						x3s,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],x3rng)
						x3s = np.sort(list(set(x3s[rng3])))
						for x3 in x3s:
							x3d,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],float(x3))
							x = xs[rng3]
							y = ys[rng3]
							np.array(zip(*sorted(zip(x,y))))
							x,rng1 = maskmaker(x,x1rng)
							a.plot(x[rng1],y[rng1],label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
			a.set_xlabel(_mathtext(self.xnames[indx]))
			a.legend(ncol=max(1,len(a.lines)/6))
		f.show()
		return f
							  
	def plot2d(self,ynames=None,aspect='auto',plot_type='mesh',cbar=True,
			   grid_size=(124,124),swap=False,**kwargs):
		"""
		Matplotlib pcolormesh or contour plots of 2D data.
		
		kwargs
		------
		- ynames : list. Strings specifying dependent data displayed.
		- aspect : str. Choose from 'auto' or 'equal' (see matplotlib axes.set_aspect).
				   float. Stretches a circle to hight=aspect*width. 
		- plot_type: str. Choose one of
		                  - "mesh" : Use matplotlib pcolormesh to plot.
						  - "line" : Use matplotlib contour to plot
						  - "fill" : Use matplotlib contourf to plot
		- cbar   : bool. Show colorbar. 
		- grid_size : tuple. Size of grid used if no regular axes x 
		                     (order is order of xnames)
		- swap   : bool. Swap x/y axes.
		**kwargs : dict. Valid pcolormesh/contour/contourf key arguments.
		
		returns
		-------
		-        : figure. 
		"""
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=[ynames]
		if self.nd != 2: raise IndexError("Data not 2D.")
			
		if self.x:
			x1,x2 = map(np.reshape,self.pts.T,[self.shape]*self.nd)
		else:
			x1,x2 = np.mgrid[self.pts[:,0].min():self.pts[:,0].max():124j,
							 self.pts[:,1].min():self.pts[:,1].max():124j]
		if swap:
			x1,x2 = x2,x1

		for name in ynames:
			try:
				if np.any(np.iscomplex(self.y[name])):
					idx = ynames.index(name)
					ynames[idx] = 'imag '+name
					ynames.insert(idx,'real '+name)
			except: # looks at new names -> KeyError
				pass

        # display data
		if 'figure' in kwargs:
			f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
		elif 'axes' in kwargs:
			if type(kwargs['axes'])==np.ndarray:
				f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
			else:
				f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
		else:
			nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
			f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=(6+2*ncol,4+2*nrow))
		
		for a,name in zip(ax.ravel(),ynames):
			print("Plotting "+name+".")
			if(plot_type=='line'):
				plotter = a.contour
			elif(plot_type=='fill'):
				plotter=a.contourf
			else:
				plotter = a.pcolormesh
			reducedname = name.replace('real ','').replace('imag ','')
			if self.x:
				ry = np.real(self.y[reducedname]).reshape(self.shape)
				iy = np.imag(self.y[reducedname]).reshape(self.shape)
			else:
				ry = griddata(self.pts,np.real(self.y[reducedname]),
							  (x1,x2),method='nearest')
				iy = griddata(self.pts,np.imag(self.y[reducedname]),
							  (x1,x2),method='nearest')
			if 'real' in name:
				pcm = plotter(x1,x2,ry,**kwargs)
			elif 'imag' in name:
				pcm = plotter(x1,x2,iy,**kwargs)
			else:
				pcm = plotter(x1,x2,ry,**kwargs)

   			a.set_xlabel(_mathtext(self.xnames[0+swap]))
			a.set_ylabel(_mathtext(self.xnames[1-swap]))
			a.set_xlim(x1.min(),x1.max())
			a.set_title(_mathtext(name))
			a.set_aspect(aspect)
			if cbar: f.colorbar(pcm,ax=a)
		f.show()
		return f

	def plot3d(self,ynames=None,npts=100,override='',**kwargs):
		"""
		Three dimensional plots. Data with xnames r,z will be plotted
		with third dimension phi where Y = Re[y]cos(n*phi)+Im[y]sin(n*phi).
		    1D data > lines in a plane in 3-space.
		    2D data > z axis becomes y value
		    3D data > chosen display type
		Must have mayavi.mlab module available in pythonpath.
		
		kwargs
		------
		- ynames   : list. Strings specifying dependent data displayed.
		- npts     : int. Plot uses npts per axis.
		- override : str. Valid mayavi.mlab function name. 
		                  Overrides default function choices (plot3d,surf,quiver3d).

		returns
		-------
		-          : figure. 
		"""
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=[ynames]
		if not mlab: raise ImportError("Unable to import enthought.mayavi.mlab.")

		print("Gridding data on "+str(npts)+" point per axis grid.")
		x = []
		for n in range(self.nd):
			x += [np.linspace(self.x[n][0],self.x[n][-1],npts)]
		zs = self.interp(x,ynames=ynames)
		if self.nd==1: x += np.zeros_like(x[0])

		# the correct 3d plotting function
		plotfunc = getattr(mlab,['plot3d','surf','quiver3d'][self.nd])

        # display data
		for a,y,lbl in zip(ax.ravel(),ys,ynames):
			if 'figure' not in kwargs: mlab.figure(bgcolor=(1,1,1))
			plotfunc(*[x+zs][0],name=lbl,**kwargs)
		f.show()
		return f

	def scatter(self,ynames=None,xname=None,x1=None,x2=None,x3=None,**kwargs):
		"""
		Scatter plot display for non-regular data. Display mimics plot1d in that it displays a single independent axis. If is another independent axis exists, and is not sliced using one of the optional keyword arguments, it s represented as the color of the plotted symbols. If the data is 3D and unsliced, preference if given in order of xnames.
		args   - ynames : Dependent data to be displayed. Default is all y data.
		       - xname  : Independent axis used in scatter plots. Default is first in xnames.
		       - x1     : Use only data with first axis clossest to value.
		       - x2     : Use only data with second axis clossest to value.
		       - x3     : Use only data with thrid axis clossest to value.
		returns         : figure.
		"""
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=[ynames]
		if not xname: xname=self.xnames[-1]
		indx = self.xnames.index(xname)

		rng = np.ones_like(self.pts[:,0],dtype=bool)
		ylbl=''
		indx2 = range(self.nd)
		indx2.remove(indx)
		if x1: 
			indx2.remove(0)
			x1 = self.pts[:,0][abs(self.pts[:,0]-x1).argmin()]
			ylbl+=', '+self.xnames[0]+'={0:.3}'.format(x1)
			rng*=self.pts[:,0]==x1
		if x2: 
			indx2.remove(1)
			x2 = self.pts[:,1][abs(self.pts[:,1]-x2).argmin()]
			ylbl+=', '+self.xnames[1]+'={0:.3}'.format(x2)
			rng*=self.pts[:,1]==x2
		if x3: 
			indx2.remove(2)
			x3 = self.pts[:,2][abs(self.pts[:,2]-x3).argmin()]
			ylbl+=', '+self.xnames[2]+'={0:.3}'.format(x3)
			rng*=self.pts[:,2]==x3
		if indx2: 
			c = self.pts[:,indx2[0]][rng]
			cset = np.array(list(set(c))).round(4)
			clbl = self.xnames[indx2[0]]
		else:
			c = np.ones_like(rng[rng])
			clbl = ''

		#double number of axes for complex data
		for name in ynames:
			try:
				if np.any(np.iscomplex(self.y[name])):
					idx = ynames.index(name)
					ynames[idx] = 'imag '+name
					ynames.insert(idx,'real '+name)
			except: # looks at new names -> KeyError
				pass
		# actual plotting
		if 'figure' in kwargs:
			f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
		elif 'axes' in kwargs:
			if type(kwargs['axes'])==np.ndarray:
				f,ax = kwargs['axes'][0].get_figure(),kwargs['axes']
			else:
				f,ax = kwargs['axes'].get_figure(),np.array(kwargs['axes'])
		else:
			nrow,ncol =  min(len(ynames),2),max(1,(len(ynames)+1)/2)
			f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=(6+2*ncol,4+2*nrow))

		for a,name in zip(ax.ravel(),ynames):
			"""if 'real' in name:
				sc=a.scatter(self.pts[:,indx][rng],np.real(self.y[name.replace('real ','')][rng]),
						 c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
			elif 'imag' in name:
				sc=a.scatter(self.pts[:,indx][rng],np.imag(self.y[name.replace('imag ','')][rng]),
						 c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
			else:
				sc=a.scatter(self.pts[:,indx][rng],np.real(self.y[name][rng]),
						 c=c,cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)"""
			y = self.y[name.replace('real ','').replace('imag ','')][rng]
			if 'real' in name:
				y = np.real(y)
			elif 'imag' in name:
				y = np.imag(y)
			sc=a.scatter(self.pts[:,indx][rng],y,c=c,
						 cmap=plt.matplotlib.cm.gist_rainbow,**kwargs)
			# for printing purposes
			line = plt.matplotlib.lines.Line2D(
					self.pts[:,indx][rng],y,visible=False)
			if 'label' in kwargs: line.set_label(kwargs['label'])
			a.add_line(line)
			
			a.set_xlabel(_mathtext(self.xnames[indx]))
			a.set_ylabel(_mathtext(name+ylbl))
			if clbl:
				f.colorbar(sc,ax=a,ticks=cset)
				cb=f.get_axes()[-1]
				cb.set_ylabel(_mathtext(clbl),fontsize='large')
		f.show()
		return f

	def slice(self,x,xname=None,ynames=[]):#npts=100,**kwargs):
		"""
		*** Under construction **
		Return 1D data object from 2D.
        Method grids new data, so beware of loosing accuracy near rational surfaces.
		args    - x       : float. Lower bound of the axis.
		kwargs  - xname   : str. x axis allong which slice is performed.
		returns -         : obj. new data nd-1 dimensional data object.
		"""
		"""ynames = np.sort(self.y.keys()).tolist()
		if not xname: 
			idx = 0
		else:
			idx = self.xnames.index(xname)
		xslc,x = np.mgrid[x:x:1j,self.x[idx-1].min():self.x[idx-1].max():npts*1j]
		dt = {'names':[self.xnames[idx-1]]+ynames,
			  'formats':[type(x[0,0])]+[type(self.y[name][0]) for name in ynames]}
		print(dt)
			
		 
		dummytxt = np.array([x][0]+[griddata(self.pts,self.y[name],(x,xslc),**kwargs)[0] for name in ynames])
		dummytxt = np.array(dummytxt.T,dtype=dt)
		dummyamble = 'Slice in {} at {}'.format(xname,x)+' from '+self.__doc__
		for k,v in self.params:
			dummyamble += '{0} = {1}'.format(k,v)
		return DataBase(dummytxt,dummyamble)"""

		# housekeeping
		if not ynames: ynames=np.sort(self.y.keys()).tolist()
		if not type(ynames) in (list,tuple): ynames=(ynames,)
		if not xname: xname=self.xnames[0]

		indx = self.xnames.index(xname)
		x = self.pts[:,indx]

		# helper to sort through axes of ND data
		def maskmaker(ax,rng):
			if type(rng) in [float,int]: 
				mask = ax==ax[abs(ax-rng).argmin()]
			elif type(rng) in [tuple,list] and len(rng)==2:
				mask = (ax>=rng[0]) * (ax<=rng[1])
			elif type(rng) in [tuple,list] and len(rng)>2:
				print("Warning: mapping all axes to regular grid.")
				if self.nd==2:
					x,ax = np.mgrid[x.min():x.max():x.size,rng[0]:rng[1]:rng[2]]
					x = x[0]
					mask = range(rng[2])
				if self.nd==3:
					raise ValueError("Griding for 3D data is under developement.")
			else:
				mask = np.ones_like(ax,dtype=bool)
			return ax,mask
		
		for a,name in zip(ax.ravel(),ynames):
			if self.nd==1:
				a.plot(x,self.y[name],label=lbl,**kwargs)
			elif self.nd==2:
				if self.x:
					x2,rng2 = maskmaker(self.x[indx-1],x2rng)
					ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
					xlbl = _mathtext(self.xnames[indx-1])
					for y,x2 in zip(ys[rng2],x2[rng2]):
						a.plot(x,y,label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
				else:
					x2s,rng2 = maskmaker(self.pts[:,indx-1],x2rng)
					xlbl = _mathtext(self.xnames[indx-1])
					x2s = np.sort(list(set(x2s[rng2])))
					for x2 in x2s:
						x,rng2 = maskmaker(self.pts[:,indx-1],float(x2))
						x = self.pts[:,indx][rng2]
						y = self.y[name][rng2]
						x,y = zip(*sorted(zip(x,y)))
						a.plot(x,y,label=lbl+', '+xlbl+'={0:.3}'.format(x2),**kwargs)
			elif self.nd==3:
				indx23 = range(3)
				indx23.remove(indx)
				if self.x:
					x2,rng2 = maskmaker(self.x[indx23[0]],x2rng)
					x3,rng3 = maskmaker(self.x[indx23[1]],x3rng)
					ys = self.y[name].reshape(self.shape).swapaxes(-1,indx)
					xlbls = [_mathtext(name) for name in self.xnames if name!=xname]
					for y2,x2 in zip(ys[rng2],x2[rng2]):
						x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
						for y,x3 in zip(y2[rng3],x3[rng3]):
							a.plot(x,y,label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
				else:
					x2s,rng2 = maskmaker(self.pts[:,indx23[0]],x2rng)
					x2s = np.sort(list(set(x2s[rng2])))
					xs = self.pts[:,indx][rng2]
					ys = self.y[name][rng2]
					xlbls = [_mathtext(xlbl) for xlbl in self.xnames if xlbl!=xname]
					for x2 in x2s:
						x2lbl = lbl+', '+xlbls[0]+'='+str(x2)
						x3s,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],x3rng)
						x3s = np.sort(list(set(x3s[rng3])))
						for x3 in x3s:
							x3d,rng3 = maskmaker(self.pts[:,indx23[1]][rng2],float(x3))
							x = xs[rng3]
							y = ys[rng3]
							x,y = zip(*sorted(zip(x,y)))
							a.plot(x,y,label=x2lbl+', '+xlbls[1]+'={0:.3}'.format(x3),**kwargs)
		
	
	

class SurfaceBase:
	"""
	An empty class with standardized visualization methods
	for an axi-symmetric surface defined by attributes r and z.
	"""

	# Projection of Data class to surface
	def get_data(self,output,ynames=None):
		"""
		Use the built in interpolation of the data object
		to project onto geometric surfaces.
		args    - data   : obj. A data object from IPEC or PENT .out file.
			    - ynames : list. Projects subset of data (default is all).
		"""
		if type(output)==str:
			output = data.read(output)
		# must have data in r,z coordinates
		if output.xnames == ['r','z']:
			r,z = output.x
			n = output.params['n']
		else:
			#try:
			#	r,z = output.y['r'],output.y['z']
			#   form new interpolators here	
			raise ValueError("Must have data in (r,z).")
		if not ynames: ynames=np.sort(output.y.keys()).tolist()
		
		for key,val in self.surfs.iteritems():
			for name in ynames:
				val[name]=[]
				for rzphi in val['rzphi']:
					nphi = len(rzphi[0][1])
					dat = output.interp(rzphi[0][0],rzphi[1][0],ynames=name)
					val[name].append(np.real(dat[name].repeat(nphi).reshape(-1,nphi)
									  *np.exp(n*rzphi[3])))
		for name in ynames:
			if '_r' in name:
				zname = name.replace('_r','_z')
				val[name.replace('_r','_parallel')] = []
				val[name.replace('_r','_perp')] = []
				for dr,dz,a in zip(val[name],val[zname],val['angle']):
					val[name.replace('_r','_parallel')].append(dr*np.cos(a)+dz*np.sin(a))
					val[name.replace('_r','_perp')].append(dr*np.sin(a)+dz*np.cos(a))

		# 2D info: projection to surface
		br,bz,bp = data.point(self.r,self.z,self.phi)
		mpts,npts = np.shape(self.phi)
		a = self.Xsection.angle.repeat(npts).reshape(-1,npts)
		sina,cosa = np.sin(a),np.cos(a)
		self.bnrm = br*sina + bz*cosa
		self.bpar = br*cosa + bz*sina
		self.btor = bp

		# names and misc
		self.dataname = data.name
		self.n = data.n
		
		# 1D info: poloidal wavelength and amplitude
		xs = self.Xsection
		tempdict = {}
		real_b = data.point(xs.r,xs.z,np.zeros_like(xs.r))
		imag_b = data.point(xs.r,xs.z,np.ones_like(xs.r)*np.pi/2/data.n)
		tempdict['real_nrm'] = real_b[0]*np.sin(xs.angle) + real_b[1]*np.cos(xs.angle)
		tempdict['imag_nrm'] = imag_b[0]*np.sin(xs.angle) + imag_b[1]*np.cos(xs.angle)
		tempdict['real_par'] = real_b[0]*np.cos(xs.angle) + real_b[1]*np.sin(xs.angle)
		tempdict['imag_par'] = imag_b[0]*np.cos(xs.angle) + imag_b[1]*np.sin(xs.angle)
		tempdict['real_tor'] = real_b[2]
		tempdict['imag_tor'] = real_b[2]

		setattr(self.Xsection,'torque',np.pi*xs.length*xs.r**2.0*(tempdict['real_nrm']*tempdict['real_tor']+tempdict['imag_nrm']*tempdict['imag_tor'])/(4.0*np.pi*1e-7))

		for coord in ['nrm','par','tor']:
			phase = np.arctan2(tempdict['imag_'+coord],tempdict['real_'+coord])
			lmda  = np.abs(xs.length * 2.0*np.pi/(phase - np.roll(phase,1)))
			amp = (tempdict['imag_'+coord]**2. + tempdict['real_'+coord]**2.)**0.5
			setattr(self.Xsection,'phase_'+coord,phase)
			setattr(self.Xsection,'lambda_'+coord,lmda )
			setattr(self.Xsection,'amp_'+coord,amp)


	
	# Built-in visualizations
	def vis1d(self,**kwargs):
		"""
		Displays the poloidal cross section of an ipec Data class
		as a line on r,z plot.
		**kwargs -    : any matplotlib.pyplot plot keyword arguments
		returns : figure. Poloidal cross sections of the vessel
		"""
		if 'axes' not in kwargs:
			pyplot.figure()
			pyplot.subplot(111,aspect='equal')
		pyplot.plot(self.Xsection.r,self.Xsection.z,label=self.name,**kwargs)
		pyplot.xlabel('R (m)')
		pyplot.ylabel('z (m)')
		pyplot.legend(loc=1,prop={'size':'small'},frameon=False)

	def vis2d(self,component='normal',pos_theta=False,**kwargs):
		"""
		Shows an 'unrolled' surface, with magnetic fields
		represented as contours of theta and phi.
		kwargs - component : string. ('normal') Choice of 'normal', 'parallel', or 'toroidal'
		         pos_theta : bool. (False) Plot theta from 0 to 360 (centers HFS mid-plane)
	                                           For contour color.
		**kwargs -        : dict. any matplotlib.pyplot contour keywords.
		
		For all three in one figure, suggested settings:
		fig,ax = pyplot.subplots(3,1,sharex=True,sharey=True, figsize=(15,12))
		pyplot.subplots_adjust(hspace=0.33,right=0.95,bottom=0.05,top=0.95)
		"""
		
		if 'dataname' not in vars(self):
			raise ValueError('Must have projected data to surface. Use .get_data().')
		if 'axes' not in kwargs:
			pyplot.figure(figsize=(13,7))
			pyplot.subplots_adjust(hspace=0.33,right=0.97,bottom=0.13,top=0.87)
			
		if component in ['n','nrm','nor','norm','normal']:
			b,lmda,amp = self.bnrm,self.Xsection.lambda_nrm,self.Xsection.amp_nrm
			label = r' $\perp$'
		elif component in ['p','par','parallel']:
			b,lmda,amp = self.bpar,self.Xsection.lambda_par,self.Xsection.amp_par
			label = r' $\parallel$'
		elif component in ['t','tor','toroidal']:
			b,lmda,amp = self.bpar,self.Xsection.lambda_par,self.Xsection.amp_par
			label = r' $\hat{\phi}$'
		else:
			raise ValueError("Choose a component ('normal,'parallel','toroidal') of data to map at surface")

		if pos_theta:
			theta = np.roll(self.theta%(2.0*np.pi),np.shape(self.theta)[0]/2,axis=0)
			b = np.roll(b,np.shape(self.theta)[0]/2,axis=0)
			cntr = pyplot.contourf(self.phi*180./np.pi,theta*180./np.pi,b*1E4,40)
		else:
			cntr = pyplot.contourf(self.phi*180./np.pi,self.theta*180./np.pi,b*1E4,40)
		cb = pyplot.colorbar(cntr,pad=0.01)
		cb.set_label('Gauss')
		pyplot.title(self.dataname+label+' at '+self.name,fontsize='xx-large')
		pyplot.xlabel(r'$\phi$',fontsize='xx-large')

		ax = pyplot.gca()
		pyplot.setp(ax.get_yticklabels(),visible=False)
		divider = make_axes_locatable(ax)
		sideax = divider.append_axes('left',size='15%',pad=0.05,sharey=ax,xlabel=r'$\lambda$')
		sideax2 = divider.append_axes('left',size='15%',pad=0.05,sharey=ax,xlabel='|B|')
		if pos_theta:
			sideax.plot(lmda,(self.Xsection.theta*180./np.pi)%360)
			sideax2.plot(amp*1e4,(self.Xsection.theta*180./np.pi)%360)
			sideax2.set_ylim(0,360)
		else:
			sideax.plot(lmda,self.Xsection.theta*180./np.pi)
			sideax2.plot(amp*1e4,self.Xsection.theta*180./np.pi)
			sideax2.set_ylim(-180,180)
		sideax.set_xlim(0,3)
		sideax.set_xticks((1,2))
		#sideax.set_yticklabels('',visible=False)
		sideax2.set_ylabel(r'$\theta$',fontsize='xx-large')
		sideax2.set_xticks([np.round(np.mean(amp*1e4),2)])
		
		

	def vis3d(self,component=None,**kwargs):
		"""
		Display, in 3D, the data as a 2D mesh.
		kwargs - component: string. ('normal') Choice of 'normal', 'parallel', or 'toroidal'
		**kwargs -        : dict. Any Mayavi.mlab mesh keywords
		returns : figure. A mayavi figure if possible. mplot3d figure if not.
		"""
		if not mlab: raise ImportError('You must be operating with mayavi module in path')

		if component in ['n','nrm','nor','norm','normal']:
			b,label = self.bnrm*1E4,self.dataname+' Normal at '
		elif component in ['p','par','parallel']:
			b,label = self.bpar*1E4,self.dataname+' Parallel at '
		elif component in ['t','tor','toroidal']:
			b,label = self.btor*1E4,self.dataname+' Toroidal at '
		else:
			b,label=np.ones_like(self.phi),''
		
		if 'figure' not in kwargs: mlab.figure(bgcolor=(1,1,1))
		mlab.mesh(self.x,self.y,self.z,scalars=b,name=label+self.name,**kwargs)
		mlab.colorbar(orientation='vertical',title='Gauss')








######################################################## helper functions

##class operations

def _data_op(self,other,fun,op):
	"""
	This function generalizes the process of creating a new class
    and filling it with the appropriate dependent data, info, and
    indepenedent data from an operation on one or more data classes.
	"""
	NewData = copy.deepcopy(self) #intialize new data object.

	# check if same dependent variables
	if self.xnames != other.xnames:
		raise ArithmeticError('Two data objects do not have same dependent variables.')
	if (np.array(self.x) != np.array(other.x)).any():
		otherinterp = other.interp(NewData.x)
	else:
		otherinterp = other.y

	for key in other.y:
		if key in NewData.y:
			NewData.y[key] = op(NewData.y[key],otherinterp[key])
		else: 
			NewData.y[key] = otherinterp[key]
    
	return NewData




def _mathtext(text):
	"""
	Helper function to convert lazy ascii conventions to 
	nicer latex strings.
	args    - text : str.
	"""

	simplemaps = ['psi','omega','nabla','cdot','kappa','theta',
				  'nu','times','perp','parallel','lambda','ln',
				  'Lambda','xi','chi','phi','mu','int', 'Re','Im',
				  'delta','epsilon']
	
	trickymaps = [('par','parallel'), ('divxprp','nablacdotxi_perp'),
				  ('exb','EtimesB'), ('lam','lambda'),
				  ('wdian','omega_*n'), ('wdiat','omega_{*T}'),
				  ('wrot','omega_phi'), ('real','Re'),
				  ('imag','Im'), ('log','ln'), ('kxprp','kappacdotxi_perp'),
				  ('prp','perp'), ('int','int ')#,('eps_','epsilon_')
				 ]

	# replace some lazy data labeling with proper math
	for tm in trickymaps:
		text=text.replace(*tm)
	# deliniate start/stop of latexing
	for sm in simplemaps:
		text=text.replace(sm,"$\\"+sm+'$')
	for word in text.split(' '):
		if '_' in word:
			newword='$'+word.replace('_','_{').replace('$','')+'}'*word.count('_')+'$'
			text=text.replace(word,newword)
		if '^' in word:
			text=text.replace(word,'$'+word.replace('$','')+'$')
	#text.replace('$$','')
	
	return text#.encode('string-escape') #raw string


######################################################## Developer functions


def _write_ellipk(npts=3e2,fname='ellipk01.dat'):
	"""
	Writes complete elliptic integral of the first 
	kind to file.
	NOTE : This is not the original routine. Could not find it.
	kwargs - npts : int.
	       - fname: str.
	"""
	# ellipk(1-1e-17)=inf
	x = np.log(1./1e-15)
	kappas = 1.0-np.exp(np.linspace(0,x,npts))/np.exp(x)
	with open(fname,'w') as f:
		f.write('{:.0f} \n'.format(npts))
		for k in kappas[::-1]:
			f.write('{:.21} '.format(k))
		f.write('\n')
		for k in kappas[::-1]:
			f.write('{} '.format(ellipk(k)))
	return True

def _write_ellipe(npts=3e2,fname='ellipe01.dat'):
	"""
	Writes complete elliptic integral of the first 
	kind to file.
	NOTE : This is not the original routine. Could not find it.
	kwargs - npts : int.
	       - fname: str.
	"""
	kappas = np.linspace(0,1,npts)
	with open(fname,'w') as f:
		f.write('{:.0f} \n'.format(npts))
		for k in kappas:
			f.write('{} '.format(k))
		f.write('\n')
		for k in kappas:
			f.write('{} '.format(ellipe(k)))
	return True
		
def _factors(n):
    return list(sorted(set(reduce(list.__add__, 
					  ([i, n//i] for i in range(1, int(n**0.5) + 1) 
					   if n % i == 0)))))
