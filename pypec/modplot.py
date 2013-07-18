#!/usr/local/env python
"""
MODULE

Collection of modified matplotlib functions 
and objects.

"""
"""
	@package nikplotlib
	@author Nikolas Logan
    @email nlogan@princeton.edu
"""

#standard matplotlib modules
import matplotlib                       # standard ploting library
from mpl_toolkits.axes_grid1 import make_axes_locatable	#for splitting axes
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

# standard python modules
import numpy as np                      # math
import copy
from scipy.interpolate import interp1d  # math
from types import MethodType            # to modify instances

pyplot = matplotlib.pyplot
show = pyplot.show
draw = pyplot.draw


#override user preferences from ~/.matplotlib/matplotlibrc file.
#really should define new functions
matplotlib.rcParams['figure.autolayout']=True # Necessary for tight_layout
matplotlib.rcParams['legend.frameon']=False
matplotlib.rcParams['legend.loc']='best'

linestyles = [ '-' , '--' , '-.' , ':' , 'None' , ' ' ,'' ]


########################################### modified pyplot functions

def subplots(nrows=1, ncols=1, sharex=True, sharey=False, squeeze=True, subplot_kw=None, powerlim=(-3,3),useOffset=False,**fig_kw):
	"""
	Matplotlib subplots with default sharex=True. 
	Additional kwargs-  powlim    : tuple. (-3,3) Axis labels use scientific notion above 10^power.
                    useOffset : bool. (False) Axis labels use offset if range<<average.
        Accepts standargs args and kwargs for pyplot.subplots.
        -
        """
	f,ax = pyplot.subplots(nrows=nrows,ncols=ncols,sharex=sharex,sharey=sharey,squeeze=squeeze,subplot_kw=subplot_kw,**fig_kw)
	# customize axes and figure
	#for a in np.array(ax).ravel():
	#   a = modaxes(a,useOffset=useOffset)
	f = _modfigure(f)
	return f,ax
subplots.__doc__ = pyplot.subplots.__doc__

def figure(num=None, figsize=None, dpi=None, facecolor=None, edgecolor=None, frameon=True, FigureClass=matplotlib.figure.Figure, **kwargs):
	f = pyplot.figure(num=num, figsize=figsize, dpi=dpi, facecolor=facecolor, edgecolor=edgecolor, frameon=frameon, FigureClass=FigureClass, **kwargs)
	f = _modfigure(f)
	return f

def colorbar(mappable=None, cax=None, ax=None, use_gridspec=True,**kw):
    pyplot.colorbar(mappable=mappable, cax=cax, ax=ax, use_gridspec=use_gridspec,**kw)
colorbar.__doc__='Modified for default use_gridspec=True. \n - \n'+pyplot.colorbar.__doc__





########################################### Figure object modifications

# internal function to apply customizations post initialization
def _modfigure(f):
	"""
	Modplot figure object modifications. These include:
	    1) canvas connection to modplot.key_press event funtion
	args   - f : obj. Figure object.
	"""
	f.canvas.mpl_connect('key_press_event',onkey) #custom hotkeys
	#pppl Figure base doesn't have show from nomachine for some reason
	#f._orig_show = MethodType(copy.deepcopy(f.show),f)
	#f.show = MethodType(_figure_show,f)
	
	return f

# customize axes when added	
def _figure_add_subplot(self,*args,**kwargs):
    ax = self._orig_add_subplot(*args,**kwargs)
    ax = _modaxes(ax)
    return ax
if not hasattr(pyplot.Figure,"_orig_add_subplot"): #prevent recursion when reloaded
    pyplot.Figure._orig_add_subplot = copy.deepcopy(pyplot.Figure.add_subplot)
pyplot.Figure.add_subplot = _figure_add_subplot

# tight_layout when show
def _figure_show(self,*args):
    self.tight_layout()
    self._orig_show(*args)
_figure_show.__doc__ = pyplot.show.__doc__
#pppl Figure base doesn't have show from nomachine for some reason
if not hasattr(pyplot.Figure,"_orig_show"): #prevent recursion when reloaded
    pyplot.Figure._orig_show = copy.deepcopy(pyplot.Figure.show)#pyplot.show)#pyplot.Figure.show)
pyplot.Figure.show = _figure_show

#colorbars are axes (for tight_layout)
def _figure_colorbar(self, mappable, cax=None, ax=None, use_gridspec=True,**kw): 
    self._orig_colorbar(mappable, cax=cax, ax=ax, use_gridspec=use_gridspec,**kw)
_figure_colorbar.__doc__='Modified default use_gridspec=True. \n - \n'+pyplot.colorbar.__doc__
if not hasattr(pyplot.Figure,"_orig_colorbar"): #prevent recursion when reloaded
    pyplot.Figure._orig_colorbar = copy.deepcopy(pyplot.Figure.colorbar)
pyplot.Figure.colorbar = _figure_colorbar


# new method to print text file of lines in figure
def printlines(self,filename,squeeze=False):
	"""
	Print all data in line plots to text file.
	The x values will be taken from the first line 
	in each axis, and other lines are interpolated
	if their x values do not match. Column lables 
	are the line labels and xlabel.
	args   - filename : str. File to print to. 
	"""
	with open(filename,'w') as f:
		data=[]
		for a in self.get_axes():
			if not squeeze: 
				data = []
				f.write(a.get_title()+'\n\n')
				f.write('  ylabel = '+a.get_ylabel().translate(None,' \${}')+'\n\n')
			for line in a.get_lines():
				x,y = line.get_xdata(),line.get_ydata()
				if not data: #get x-axis
					data.append(x)
					f.write('{0:>25s}'.format(a.get_xlabel().translate(None,' \${}')))
				label = line.get_label().translate(None,' \${}')
				if squeeze: label = a.get_ylabel().translate(None,' \${}')+label
				f.write('{0:>25s}'.format(label))
				# standard axis
				if np.any(x!=data[0]):
					fit = interp1d(x,y,bounds_error=False,fill_value=np.nan)
					data.append(fit(data[0]))
				else:
					data.append(y)
			if not squeeze:
				f.write('\n ')
				data=np.array(data).T
				for row in data:
					row.tofile(f,sep=' ',format='%24.9E')
					f.write('\n ')
				f.write('\n')
		if squeeze:
			f.write('\n ')
			data=np.array(data).T
			for row in data:
				row.tofile(f,sep=' ',format='%24.9E')
				f.write('\n ')
			f.write('\n')
	print("Wrote lines to "+filename)
	return True
pyplot.Figure.printlines = printlines


########################################### Axes object modifications


# internal function to apply customizations post initialization
def _modaxes(ax,useOffset=False):
	"""
	Modplot axis object modifications. These include:
	    1) Optional offset
	    2) Ticker format power lims (-3,3)
	    3) New downsample instance, connected to x-axis rescale
	    4) New plot instance with complex arg handling
	    
	args   - ax        : obj. Axes object.
	kwargs - useOffset : bool. Offset axis zero.
	"""
	afmt = matplotlib.ticker.ScalarFormatter(useOffset=useOffset)
	afmt.set_powerlimits((-3,3))
	ax.callbacks.connect('xlim_changed', ax.downsample)
	ax.yaxis.set_major_formatter(afmt)
	ax.xaxis.set_major_formatter(afmt)
	return ax
	
# down sample all lines in axes
def _axes_downsample(self,ax,npts=1000):
	"""
	Attempt to down sample original data in
	each Line2D object in lines.
	args -   ax   : obj. Initialized Axes object containing lines
	kwargs - npts : int. Number of points plotted in each line
	"""
	for l in ax.lines:
		try:
			l.downsample(npts=npts)
		except:
			pass
	#draw()
pyplot.Axes.downsample = _axes_downsample

# complex value plotting
def _axes_plot(self, *args, **kwargs):
	"""
	Modified axes line plotting, with suport for complex plotting.
	Original matplotlib Axes plot method now located at _orig_plot.
	All original plot args and kwargs accepted.
	-
	"""
	linestyles = matplotlib.lines.lineStyles.keys()
	newargs = map(np.real,args)
	lines1=self._orig_plot(*newargs,**kwargs)
	lines2 = []
	if any([any(np.ravel(np.iscomplex(arg))) for arg in args]):
		if len(args)==1:
			newargs = np.imag(args)
		else:
			newargs[1]=np.imag(args[1])
		lines2=self._orig_plot(*newargs,**kwargs)
		for l1,l2 in zip(lines1,lines2): 
			l2.set_color(l1.get_color())
			ils = linestyles.index(l1.get_linestyle())
			l2.set_linestyle(linestyles[ils+1])
			l2.set_label(r'$\Im$ '+l1.get_label())
			l1.set_label(r'$\Re$ '+l1.get_label())
		pyplot.gca().legend(ncol=max(1,len(self.lines)/6))
	#draw()
	return lines1+lines2
_axes_plot.__doc__+=pyplot.plot.__doc__
if not hasattr(pyplot.Axes,"_orig_plot"): #prevent recursion when reloaded
    pyplot.Axes._orig_plot = copy.deepcopy(pyplot.Axes.plot)
pyplot.Axes.plot = _axes_plot


def _center_axes(ax):
	ax.spines['left'].set_position('center')
	ax.spines['left'].set_linewidth(2)
	ax.spines['right'].set_color('none')
	ax.spines['bottom'].set_position('center')
	ax.spines['bottom'].set_linewidth(2)
	ax.spines['top'].set_color('none')
	ax.spines['left'].set_smart_bounds(True)
	ax.spines['bottom'].set_smart_bounds(True)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
pyplot.Axes.center_axes = _center_axes


########################################### modified pyplot functions

# downsample data from the full set used to initialize the line
def line_downsample(line,npts=1000):
        """
        Downsample line data for increased plotting speed.
        kwargs - npts : int. Number of data points within axes xlims.
        """
	lims = line.axes.get_xlim()
	if not hasattr(line,'_xinit'):
		line._xinit = line.get_xdata()
		line._yinit = line.get_ydata()
	window = (line._xinit>=lims[0]) & (line._xinit<=lims[1])
	#extend out edge
	window[window.argmax()-1]=True
	window[window.argmax()+window[window.argmax():].argmin()]=True
	step = int(line._xinit[window].size/npts) + 1
	line.set_xdata(line._xinit[window][::step])
	line.set_ydata(line._yinit[window][::step])
pyplot.matplotlib.lines.Line2D.downsample = line_downsample





########################################### my own functions

	
def plot_axes(ax, fig=None, geometry=(1,1,1)):
	"""
	Re-create a given axis in a new figure. This allows, for instance,
	a subplot to be moved to its own figure where it can be manipulated
	and/or saved independent of the original.
	args   - ax : obj. An initialized Axes object  
	kwargs - fig: obj. A figure in which to re-create the axis
  	       - geometry : tuple. Axes geometry of re-created axis
    returns           : figure.
	"""
	if fig is None:
		fig = figure()
	a2 = copy.copy(ax)
	a2.set_figure(fig)
	a2 = fig.add_axes(a2)
#	for k,v in vars(ax).iteritems():
#		if '_' not in k:
#			try:
#				setattr(ax2,k,v)
#			except:
#				pass
	if ax.get_geometry() != geometry :
		a2.change_geometry(*geometry)
#	draw()
	fig.show()
	return fig


def onkey(event):
	"""
    Function to connect key_press_event events from matplotlib
    to custom functions.

    Matplotlib defaults (may be changed in matplotlibrc):
    keymap.fullscreen : f               # toggling
    keymap.home : h, r, home            # home or reset mnemonic
    keymap.back : left, c, backspace    # forward / backward keys to enable
    keymap.forward : right, v           #   left handed quick navigation
    keymap.pan : p                      # pan mnemonic
    keymap.zoom : o                     # zoom mnemonic
    keymap.save : s                     # saving current figure
    keymap.quit : ctrl+w                # close the current figure
    keymap.grid : g                     # switching on/off a grid in current axes
    keymap.yscale : l                   # toggle scaling of y-axes ('log'/'linear')
    keymap.xscale : L, k                # toggle scaling of x-axes ('log'/'linear')
    keymap.all_axes : a                 # enable all axes

    My custom function mapping:
    popaxes : n                         #Creates new figure with current axes
    tighten : t                         #Call tight_layout bound method for figure
	"""
	if event.key=='n':
		if event.inaxes:
			f = plot_axes(event.inaxes)
	if event.key=='t':
		f = pyplot.gcf()
		f.tight_layout()
		f.show()
	return



	
