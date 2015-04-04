#!/usr/local/env python

"""
EXE

GUI for gpec package.

"""
"""
	@package pypec
	@author NC Logan
	@email nlogan@pppl.gov
"""
    
from enthought.traits.api import HasTraits,Button, Instance, Bool,Str,List,Int,Float,Complex,File,Directory,List,Enum
from enthought.traits.ui.api import Item, View,NoButtons,Group,InstanceEditor,Handler,UItem

from string import join                 # string manipulation
import os

import data,gpec,_defaults_,_tooltips_
from collections import OrderedDict

################################ Handlers

class UpdateHandler(Handler):
	def object_updates_changed(self, info):
		print(info.object.update)
		if info.object.update:
			info.object.update=False
			info.object.configure_traits()
			info.ui.dispose()

################################ GUI Elements
class InputPanel(HasTraits):
	def __init__(self,inputs):
		"""
		#Display russian-doll tabs from dict of namelists.
		"""
		self.items = []
		self.layout = 'normal'
		for key,val in inputs.iteritems():
			self.add_item(key,val)
	
	def trait_view(self,name=None, view_elements=None):
		return View(Group(*self.items,layout=self.layout,
						   show_labels=self.layout=='normal'),
					resizable=True,height=600,buttons=NoButtons,
					scrollable=True)


	def add_item(self,key,val,update=False):
		style = 'simple'
		if type(val) in [dict,OrderedDict]:
			self.add_trait(key,InputPanel(val))
			initializer = InputPanel
			style='custom'
			self.layout='tabbed'
		else:
			if 'file' in key:
				initializer = File
			else:
				types = [bool,str,int,float,complex,list]
				Types = [Bool,Str,Int,Float,Complex,List]
				initializer = Types[types.index(type(val))]
		#setattr(self,key,initializer(val))
		#trait names cannot be any old dict keyword... same as class attrs
		traitname = key.translate(None,' .,()[]{}/:;*%$!@#%^&-+')
		#self.add_trait(key+'_show',Bool(True))
		test = [item.label==key for item in self.items]
		if any(test):
			setattr(self,traitname,initializer(val))
			self.items[test.index(True)] = Item(traitname,label=key,style=style)
		else:
			self.add_trait(traitname,initializer(val))
			self.items.append(Item(traitname,label=key,style=style))
			if key in _tooltips_.alltips: 
				self.items[-1].tooltip = _tooltips_.alltips[key]
   		if update:
			self.view = self.new_view()

	def _todict(self):
		d=OrderedDict()
		for item in self.items:
			if self.layout=='tabbed':
				d[item.label]=getattr(self,item.name)._todict()
			else:
				d[item.label]=getattr(self,item.name)
				if 'file' in item.label:
					d[item.label] = str(d[item.label])
		return d

class NewParameter(HasTraits):
	name = Str
	value= Str
	ok = Button()
	view = View('name','value','ok',buttons=NoButtons)
	def __init__(self,inputs):
		self.add_trait('inputs',inputs)
	def _ok_fired(self):
		self.inputs.add_item(self.name,self.value)

class ControlPanel(HasTraits):
	run = Button()
	load_all = Button()
	load_dir = Directory(_defaults_.inputdir)
	load = Button()
	load_file = File(_defaults_.inputdir+'/sol.in')
	qsub= Bool()
	roc = Bool()
	rundcon = Bool(True)
	runipec = Bool()
	runpent = Bool()
	rundir = File("/p/gpec/users/nlogan/ipec_3.00/nlogan/rundir/LinuxLahey64")
	location = Directory(".")
	mailon = Enum(_defaults_.mailon,'','a','b','e','ab','ae','be','abe')
	email = Str(_defaults_.email)
	inputpanel = Instance(InputPanel({'No Defaults Set':''}),())

	emacs = Button()

	update = Bool(False)

	view=View(Item('rundir',label='Run Directory',width=300),
			  Item('location',label='Location',width=300),
			  Group(Item('mailon',label='Mail On'),
					Item('email',label='email',width=200,visible_when="mailon!=''"),
					orientation='horizontal'),
			  Group(Item('load_all',show_label=False),
					Item('load_dir',show_label=False,width=300),
					columns=2,label='New Inputs from Directory',
					show_border=True),
			  Group(UItem('load',show_label=False),
					Item('load_file',show_label=False,width=300),
					columns=2,label='New NameList from File',
					show_border=True),
			  Group(Item('rundcon',label='DCON'),
					Item('runipec',label='IPEC'),
					Item('runpent',label='PENT'),
					orientation='horizontal'),
			  Group(Item('qsub',label='Submit to Cluster'),
			  		Item('roc',label='Wait until Complete'),
			  		orientation='horizontal'),
			  Item('inputpanel',style='custom',show_label=False,
				   editor=InstanceEditor(editable=True)),
			  Item('emacs',label='Brute Edit Namelists',show_label=False),
			  Item('run',show_label=False),
			  scrollable=True, handler=UpdateHandler(),
			  resizable=True,height=1000, width=525,
			  buttons=NoButtons,title='IPEC 3.00 Package')

	def _run_fired(self):
		print('runing program(s)')
		gpec.run(loc=str(self.location),rundir=str(self.rundir),qsub=self.qsub,
				 return_on_complete=self.roc,mailon=str(self.mailon),email=str(self.email),
				 rundcon=self.rundcon,runipec=self.runipec,runpent=self.runpent,
				 **self.inputpanel._todict())

	def _load_fired(self):
		f = str(self.load_file)
		name = f.split('/')[-1][:-3]
		self.inputpanel.add_item(name,gpec.namelist.read(f))
		self.update=True
		
	def _load_all_fired(self):
		d = str(self.load_dir)
		ins = gpec.Inputset(indir=d).infiles
		for k,v in ins:
			if k in _defaults_.inputs:
				self.inputpanel.add_item(k,v)
		self.update=True

	def _emacs_fired(self):
		gpec.run(loc=_defaults_.tempdir,rundir=str(self.rundir),qsub=False,
				 return_on_complete=False,rundcon=False,runipec=False,runpent=False,
				 mailon=str(self.mailon),email=str(self.email),
				 **self.inputpanel._todict())
		os.system('emacs '+_defaults_.tempdir+'/*.in')
		ins = gpec.InputSet(indir=_defaults_.tempdir).infiles
		ordered_ins = OrderedDict()
		for key in _defaults_.inputs:
			ordered_ins[key]=ins[key]
		for key,val in ins.iteritems():
			if key not in _defaults_.inputs:
				ordered_ins[key]=val
		print('forming new inputs')
		self.inputpanel = InputPanel(ordered_ins)
		print('updating')
		self.update=True
		

def main():
	"""
	Run GUI.
	"""
	default = OrderedDict()
	ins = gpec.InputSet(indir=_defaults_.inputdir).infiles
	for name in _defaults_.inputs:
		if name in ins: default[name] = ins[name]
	inputs = InputPanel(inputs=default)
	control = ControlPanel(inputpanel=inputs)
	control.configure_traits()
	#MainWindow().configure_traits()

if __name__ == '__main__':
	main()
