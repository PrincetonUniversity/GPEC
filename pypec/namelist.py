"""
MODULE

Read and write fortran namelists.

EXAMPLE:
import namelist
nl=namelist.read('filename.in')
namelist.write(nl,'copyoffile.in')
"""
"""
	@package ipec
	@author NC Logan
	@email nlogan@pppl.gov
"""

from types import InstanceType
from string import strip
from collections import OrderedDict

class Objectify(object):
    def __init__(self, d):
    	"""Class base to convert iterable to instance."""
        for a, b in d.iteritems():
            if isinstance(b, (list, tuple)):
               setattr(self, a, [Objectify(x) if isinstance(x, dict) else x for x in b])
            else:
               setattr(self, a, Objectify(b) if isinstance(b, dict) else b)
    def _todict(self): 
			"""Converts instance back to dictionary."""
			d = OrderedDict()
			for attr in dir(self):
				if not attr.startswith('__') and not attr=='_todict':
					val = getattr(self, attr)
					if isinstance(val, (list, tuple)):
						d[attr] = [x._todict() if isinstance(x, Objectify) else x for x in val]
					else:
						d[attr] = val._todict() if isinstance(val, Objectify) else val 
			return d



def _string_to_type(str):
	"""
	Convert a string representing a fortran namelist value
	to the appropriate python type object.
	args    - str  : str. String to convert
	returns        : unknown. Python object of appropriate type.
	"""
	try:
		pyobj=int(str)
	except ValueError:
		try:
			pyobj=float(str)
		except ValueError:
			if str.lower().startswith('.t') or str.lower().startswith('t'):
				pyobj=True
			elif str.lower().startswith('.f') or str.lower().startswith('f'):
				pyobj=False
			#complex written as tuple
			elif str.startswith('(') and str.endswith(')'):
				fpair = map(float,str.translate(None,'() ').split(','))
				pyobj=complex(*fpair)
			#lists written using spaces or repetition
			elif '*' in str or ' ' in str and str.translate(None,'.* ').isalnum():
				pyobj=[]
				for i in str.split():
					if '*' in i:
						n,v = i.split('*')
						pyobj+=int(n)*[_string_to_type(v)]
					else:
						pyobj+=[_string_to_type(i)]
			#default to string
			else:
				pyobj=str.replace("'","").replace('"','')
	return pyobj



def read(filename,out_type='dict'):
	"""
	Read in a file. Must have namelists in series, not embedded.
	Assumes name of list preceded by '&' and list ends with '\'.
	args   - filename : str. Path to input namelist.
	kwargs - out_type : str. Can be either 'dict' or 'obj'.

	returns           : dict (object). Top level keys (attributes) are the
	                    namelists with sub-key (sub-attribute) parameters and
                        their associated values.
						   - Note that dict is actuall "OrderedDict" type from 
						   collections module. Object is un-ordered.
	"""
	
	if not out_type in ['dict','obj']:
		raise ValueError("Must output 'dict' or 'obj' type.")
	nldict = OrderedDict()
	
	with open(filename,'r') as f:
		filestr = f.read()
	#remove tabs and split lines
	nlist = filestr.translate(None,'\t').replace('&','&\n').split('\n')
	#prevent errors in searching for & or /, without condensing list values
	nlist=map(strip,nlist)

	start=0
	stop=0
	it=0
	while stop < len(nlist):
		it+=1
		start = nlist.index('&',stop)+1
		stop = nlist.index('/',stop+2)-1
		name = nlist[start].strip().upper()
		nldict[name] = OrderedDict()
		#loop through each list, skiping empty lines
		for item in [line for line in nlist[start+1:stop+1] if line.strip()]:
			iname,ival = map(strip,item.split('='))
			try:
				nldict[name][iname]=_string_to_type(ival)
			except:
				raise ValueError("Failed to read "+name+" due to line: "+item)
		try:
			test=nlist.index('&',stop)
		except ValueError:
			stop=len(nlist)
		if it>100: stop = len(nlist) #prevent inf loop if bug

	if out_type=='obj':
		return Objectify(nldict)
	return nldict


def write(nl,filename):
	"""
	Write namelist object to file for fortran interface.
	args   - nl       : object. namelist object from read() 
						can be dictionary or Objectify type class.
		   - filename : str. Path of file to be written.

	returns           : bool. True on completion.
	"""
	if type(nl) == Objectify:
		print('Warning: lists from objects are un-ordered.')
		nl = nl._todict()
	
	with open(filename,'w') as f:
		f.write('\n') #vacuum code requires first line in vac.in
		for key in nl:
			f.write('&'+key+'\n')
			for subkey in nl[key]:
				if type(nl[key][subkey])==bool:
					f.write('   '+subkey+' = .'+str(nl[key][subkey]).upper()+'.\n')
				elif type(nl[key][subkey])==complex:
					fortranfmt = str(nl[key][subkey]).replace('+',',').replace('j','')
					f.write('   '+subkey+' = '+fortranfmt+'\n')
				elif type(nl[key][subkey])==str:
					f.write('   '+subkey+' = "'+str(nl[key][subkey])+'"\n')
				elif type(nl[key][subkey])==list:
					spacedlist = str(nl[key][subkey]).translate(None,',[]')
					f.write('   '+subkey+' = '+spacedlist+'\n')
				else:
					f.write('   '+subkey+' = '+str(nl[key][subkey])+'\n')
			f.write('/\n')
		
	return True
		
		
		
		
				
