#!/usr/local/env python
"""
MODULE

Collection synthetic diagnostics for analysis of data objects.

Base classes define line and surface diagnostic operations and visualization.

Specific diagnostics are read from local private ascii files.



EXAMPLES
--------

>>> mydata, = read('_example_data.dat')
Casting table 1 into Data object.
#>>> d3d.Vessel.set_data(mydata)
....
"""
# standard module imports
import os
import numpy as np
from scipy.interpolate import interp1d,interpn
from string import join                 # string manipulation
from mpl_toolkits.axes_grid1 import make_axes_locatable    #for splitting axes

#in this package
import data
#import _syntheticdata_


######################################################## Variables

plt = data.plt
_mathtext = data._mathtext
MaxNLocator = data.plt.matplotlib.ticker.MaxNLocator
packagedir = join(os.path.abspath(__file__).split('/')[:-2],'/')+'/' # one directory up

dtor = np.pi/180

######################################################## 1D classes

class LinearBase:
    """
    General base type for 1D synthetic diagnostics, providing
    data retreaval and visualization techniques.
    
    Specifc coordinate classes should provide the following
    standard attributes:
    pts    - dict. All relavent coordinates.
    xnames - list. Default coordinates.

    """
    def set_data(self,data,quiet=True,**kwargs):
        """
        Interpolate data from Data object to diagnostic points.

        args
        ----
        - data     : obj. Data object

        **kwargs   : dict. Passed to Data.interp
        """
        pts = np.array([self.pts[name] for name in self.xnames])
        
        # Don't spend time interpolating auto-generated phases and amplitudes
        keys = data.y.keys()
        for k in data.y.keys():
            if (('angle ' in k and k.replace('angle ','') in keys) or
                ('|' in k and k.rstrip('|').lstrip('|') in keys)):
                keys.remove(k)
        if not quiet: print("Interpolating {:} to {:}".format(keys,self.xnames[:2]))

        pts = zip(*[self.pts[name] for name in self.xnames[:2]])
        self.data = data.interp(pts,quiet=quiet,**kwargs).y
        
        # mode number
        if 'n' in data.params:
            self.data_n = data.params['n']
        else:
            self.data_n = 0
 
        # some quick default manipulation
        keys = self.data.keys()
        for k in keys:
            if np.any(np.iscomplex(self.data[k])):
                self.data['angle '+k] = np.angle(self.data[k]) #phase
                self.data['| '+k+' |'] = np.abs(self.data[k])   #amplitude

        return

    def mean(self,datakeys=None):
        """
        Mean value of data at diagnostic.

        kwargs
        ------
        - datakeys  : list. Any number of keys from data attribute.

        returns
        -------
        -           : dict. Means of each data key.
        """
        # housekeeping
        if not datakeys: datakeys=np.sort(self.data.keys()).tolist()
        if not type(datakeys) in (list,tuple): datakeys=(datakeys,)
        # calculation
        datamean = {}
        for k in datakeys:
            datamean[k] = np.mean(self.data[k])
        return datamean

    def plot(self,coord,datakeys=None,**kwargs):
        """
        Visualize data in 1D plots.

        args
        ----
        - coord     : str. x coordinate name (must be key of pts attribute).
        
        kwargs
        ------
        - datakeys      : list. Any number of keys from data attribute.
        **kwargs    : dict. Valid matplotlib pyplot.plot keyword arguments.
        """
        # housekeeping
        if not datakeys: datakeys=np.sort(self.data.keys()).tolist()
        if not type(datakeys) in (list,tuple): datakeys=(datakeys,)
        x = self.pts[coord]
        
        # display data
        if 'figure' in kwargs:
            f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
        elif 'axes' in kwargs:
            ax = np.array(kwargs['axes']).ravel()
            f = ax[0].get_figure()
        else:
            nrow,ncol =  min(len(datakeys),2),max(1,(len(datakeys)+1)/2)
            f,ax = plt.subplots(nrow,ncol,squeeze=False,figsize=plt.rcp_size*(ncol,nrow))
        ax = np.array(ax)
        for kwarg in ['figure','axes']:
            if kwarg in kwargs: del kwargs[kwarg]

        if 'label' in kwargs:
            ulbl = kwargs['label']+', '
            del kwargs['label']
        else:
            ulbl = ''
        
        for a,name in zip(ax.ravel(),datakeys):
            lbl = ulbl + _mathtext(name)
            a.plot(x,self.data[name],label=lbl,**kwargs)
            a.set_xlabel(_mathtext(coord))
            a.set_title(self.name)
            a.legend(ncol=max(1,len(a.lines)/6))
            f.show()
        return f

class RZCrossSection(LinearBase):
    def __init__(self,r,z,mtheta=540,r0=1.69,z0=0,t0=-np.pi/2,close=False,name=''):
        """
        Fill in a profile in poloidal angle from (r0,z0) by linearly 
        extrapolating r(theta) andf z(theta).

        Known Issue: all points must have a line of sight to the center,
        preventing overhanging features such as a diverter shelf.
        
        args
        ----
        - r      : list. Major radius points.
        - z      : list. Vertical points.

        kwargs
        ------
        - mtheta : int. Number of points in theta
        - r0     : float. Major radius of centroid
        - z0     : float. Vertical position of centroid
        - name   : str. Name of diagnostic.
        - close  : bool. Connected poloidal profile (0->2pi)

        returns
        --------
        -        : obj. Class containing attributes r,z,theta,angle,length
        """
        self.name = name

        t = (np.arctan2(z-z0,r-r0)-t0)%(2*np.pi)
        t,r,z = np.array( zip( *sorted(zip(t,r,z)) ) ) #sorts by theta

        if close:
            t = np.hstack((t,t[0]+2*np.pi))
            r = np.hstack((r,r[0]))
            z = np.hstack((z,z[0]))
        fun_r = interp1d(t,r)
        fun_z = interp1d(t,z)
    
        theta = np.linspace(t[0],t[-1],mtheta)
        # general points dictionary attribute
        self.xnames = ['r','z']
        self.pts = {}
        self.pts['r'] = fun_r(theta)
        self.pts['z'] = fun_z(theta)
        self.pts['theta'] = (np.arctan2(self.pts['z']-z0,self.pts['r']-r0)-t0)%(2*np.pi)
        self.pts['theta (deg)'] = self.pts['theta']*180/np.pi -90

        # closer was causing plot issues... resort by theta
        t,r,z = self.pts['theta'],self.pts['r'],self.pts['z']
        t,r,z = np.array( zip( *sorted(zip(t,r,z)) ) ) #sorts by theta
        if close:
            t = np.hstack((t,t[0]+2*np.pi))
            r = np.hstack((r,r[0]))
            z = np.hstack((z,z[0]))
        self.pts['theta'],self.pts['r'],self.pts['z'] = t,r,z
        self.pts['theta (deg)'] = self.pts['theta']*180/np.pi# -90

        # check extrapolations
        #f,ax = plt.subplots()
        #ax.plot(theta)
        #ax.plot(self.pts['theta'])
        #ax.plot(self.pts['r'])
        #ax.plot(self.pts['z'])
        #ax.set_title(name)
        
        # attributes specific to r,z coords
        r,z = self.pts['r'],self.pts['z']
        #dz,dr = (z-np.roll(z,1)),(r-np.roll(r,1))
        dz,dr = np.gradient(z),np.gradient(r)
        self.angle = np.arctan2(dz,dr)
        #self.length = ( (z-np.roll(z,1))**2. + (r-np.roll(r,1))**2. )**0.5
        self.length = np.sqrt(dz**2 + dr**2)
        self.pts['length'] = np.cumsum(self.length)
        
    def set_rzcharacteristics(self):
        """
        Converts r,z vector data identified by suffixes '_r' and '_z'
        to perpendicular and parallel vector data. New data
        is added to the data dictionary with '_prp' and '_par' suffixes.

        Also calculates local wavelength of data along the linear diagnostic
        using the unique RZ length attribute.
        """
        # perp and parallel vector data
        sina,cosa = np.sin(self.angle),np.cos(self.angle)
        keys = [k.replace('_r','') for k in self.data.keys() if k.endswith('_r')]
        for k in keys:
            self.data[k+'_prp'] =-1*(self.data[k+'_r']*sina-self.data[k+'_z']*cosa)
            self.data[k+'_par'] = 1*(self.data[k+'_r']*cosa+self.data[k+'_z']*sina)
            for suffix in ['_prp','_par']:
                if np.any(np.iscomplex(self.data[k+suffix])):
                    self.data['angle '+k+suffix] = np.angle(self.data[k+suffix]) #phase
                    self.data['| '+k+suffix+' |'] = np.abs(self.data[k+suffix])   #amplitude

        # wavelength
        keys = self.data.keys()
        for k in keys:
            if k.startswith('angle '):
                newk = k.replace('angle ','lambda ')
                dphase = np.gradient(self.data[k])#-np.roll(self.data[k],1)
                for i in range(len(dphase)):
                    if np.abs(dphase[i])>np.pi/2:
                        dphase[i] = dphase[i-1]
                self.data[newk] = np.abs(self.length*2*np.pi/dphase)
                # clean up ends
                #self.data[newk][0] = self.data[newk][1]
                #self.data[newk][-1] = self.data[newk][-2]

######################################################## 2D classes

class SurfaceBase:
    """
    General base type for 2D synthetic diagnostics, providing
    data retreaval and visualization techniques.
    
    Specifc coordinate classes should provide the following
    standard attributes:
    pts    - dict. All relavent coordinates.
    xnames - list. Default coordinates.
    """

    def set_data(self,data,quiet=True,**kwargs):
        """
        Interpolate data from Data object to diagnostic points.

        args
        ----
        data     - obj. Data object

        **kwargs - dict. Passed to Data.interp
        """
        pts = zip([self.pts[name] for name in self.xnames])
        self.data = data.interp(pts,quiet=quiet,**kwargs).y

        # some quick default manipulation
        keys = self.data.keys()
        for k in keys:
            if np.any(np.iscomplex(self.data[k])):
                self.data['angle '+k] = np.angle(self.data[k]) #phase
                self.data['|'+k+' |'] = np.abs(self.data[k])   #amplitude

        return

    def plot(self,coords,datakeys=None,aspect='auto',plot_type='mesh',cbar=True,
             title='',units='',amplitude=False,phase=False,wavelength=False,**kwargs):
        """
        Matplotlib pcolormesh or contour plots of 2D data.
        
        kwargs
        ------
        - coords   : list. 2 keys from pts attribute for x,y axes.
        - datakeys : list. Strings specifying dependent data displayed.
        - aspect : str. Choose from 'auto' or 'equal' (see matplotlib axes.set_aspect).
                   float. Stretches a circle to hight=aspect*width. 
        - plot_type: str. Choose one of
                          - "mesh" : Use matplotlib pcolormesh to plot.
                          - "line" : Use matplotlib contour to plot.
                          - "fill" : Use matplotlib contourf to plot.
                          - "surf" : Use matplotlib Axes3D plot_surface to plot (forces x,y,z coords).
        - cbar   : bool. Show colorbar. 
        - title  : str. Added to front of the axes title(s).
        - units  : str. Labels the colorbar.
        - amplitude  : bool. Appends 1D plot of amplitude to y axis.
        - phase  : bool. Appends 1D plot of pahse to y axis.
        - wavelength  : bool. Appends 1D plot of wavelength to y axis.
        **kwargs : dict. Valid pcolormesh/contour/contourf key arguments.
        
        returns
        -------
        -        : figure. 
        """
        if not datakeys: datakeys=np.sort(self.data.keys()).tolist()
        if not type(datakeys) in (list,tuple): datakeys=[datakeys]
            
        x1,x2 = [self.pts[k] for k in coords]

        # display data
        if 'figure' in kwargs:
            f,ax = kwargs['figure'],np.array(kwargs['figure'].get_axes())
            del kwargs['figure']
        elif 'axes' in kwargs:
            ax = np.array(kwargs['axes']).ravel()
            f = ax[0].get_figure()
            del kwargs['axes']
        else:
            nrow,ncol =  min(len(datakeys),2),max(1,(len(datakeys)+1)/2)
            subplot_kw = {}
            if plot_type=='surf': subplot_kw={'projection':'3d'}
            f,ax = plt.subplots(nrow,ncol,squeeze=False,subplot_kw=subplot_kw,
                                figsize=plt.rcp_size*(ncol,nrow))#(4+2*ncol,3+2*nrow))
        ax = np.array(ax)
        for kwarg in ['figure','axes']:
            if kwarg in kwargs: del kwargs[kwarg]

        for a,name in zip(ax.ravel(),datakeys):
            print("Plotting "+name+".")
            y = self.data[name]
            if(plot_type=='line'):
                plotter = a.contour
            elif(plot_type=='fill'):
                plotter=a.contourf
            elif plot_type=='surf':
                a.grid()
                plotter=a.plot_surface
                cmap = plt.matplotlib.cm.jet
                kwargs['facecolors'] = cmap(self.data[name])#cmap(self.data[name])
                defaults = {'linewidth':0,'rstride':1,'cstride':1,'antialiased':False,'cmap':cmap}
                for k in defaults:
                    if not k in kwargs: kwargs[k]=defaults[k]
                x1,x2,y = self.pts['x'],self.pts['y'],self.pts['z']
            else:
                plotter = a.pcolormesh
            pcm = plotter(x1,x2,y,**kwargs)

            a.set_xlabel(_mathtext(coords[0]))
            a.set_ylabel(_mathtext(coords[1]))
            a.set_xlim(x1.min(),x1.max())
            a.set_ylim(x2.min(),x2.max())
            ttl = _mathtext(name+bool(self.name)*(' on '+self.name))
            if title: ttl = title+' '+ttl
            a.set_title(ttl)
            a.set_aspect(aspect)
            if cbar:
                cb = f.colorbar(pcm,ax=a)
                if units: cb.set_label('('+units+')')
            
            if amplitude or phase or wavelength:
                ylim = a.get_ylim()
                ylbl = a.get_ylabel()
                divider = make_axes_locatable(a)
                a.set_ylabel('')
                f.show()
                #a.yaxis.set_ticklabels([])
            if amplitude:
                k = '| '+name+' |'
                a1 = divider.append_axes('left',size='25%',pad=0.05,sharey=a)
                y = self.CrossSection.pts[coords[1]]
                x  = self.CrossSection.data[k]
                a1.plot(x,y)
                a1.set_xlabel(_mathtext(k))
                a1.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='lower'))
                a1.set_ylim(ylim)
                if not (phase or wavelength):
                    a1.set_ylabel(ylbl)
            if phase:
                k = 'angle '+name
                a1 = divider.append_axes('left',size='25%',pad=0.05,sharey=a)
                y = self.CrossSection.pts[coords[1]]
                x  = self.CrossSection.data[k]
                a1.plot(x,y)
                a1.set_xlabel(_mathtext(k))
                a1.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='lower'))
                a1.set_ylim(ylim)
                if not wavelength:
                    a1.set_ylabel(ylbl)
            if wavelength:
                k = 'lambda '+name
                a1 = divider.append_axes('left',size='25%',pad=0.05,sharey=a)
                y = self.CrossSection.pts[coords[1]]
                x  = self.CrossSection.data[k]
                a1.plot(x,y)
                a1.set_xlabel(_mathtext(k))
                a1.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='lower'))
                a1.set_ylim(ylim)
                a1.set_ylabel(ylbl)
            # may need to do a.yaxis.set_ticklabels([]) for inner axes
            
        f.show()
        return f
    
    def plot3d(self,datakeys=[],**kwargs):
        """
        Plot 3D structure.
        
        **Key Word Arguments:**
          datakeys : list. 
            Strings specifying dependent data displayed.
        
        **Returns:**
          figure.
        
        """
        # Just the surface: use color and alpha kwargs
        if not datakeys:
            data.mmlab.figure(bgcolor=(1,1,1))
            data.mmlab.mesh(self.pts['x'],self.pts['y'],self.pts['z'],**kwargs)
        # the interpolated data
        for key in datakeys:
            data.mmlab.figure(bgcolor=(1,1,1))
            data.mmlab.mesh(self.pts['x'],self.pts['y'],self.pts['z'],scalars=self.data[key],**kwargs)
        
        return
        
class AxisymSurface(SurfaceBase):
    """
    Two dimensional axisymmetric surface in (r,z,phi) created
    by extending an RZCrossSection in phi.
    
    Data manipulation and visualization are handled by base type
    SurfaceBase.
    
    """
    def __init__(self,RZCrossSection,phi_extent=(0,2*np.pi),nphi=180):
        """
        Fill a toroidally symmetric surface.

        args
        ----
        - r       : list. Major radius nodes of poloidal cross section.
        - z       : list. Vertical position nodes of poloidal cross section.
        
        kwargs
        ------
        - nphi    : int. Number of toroidal points on surface.
        - name    : str. Name of surface object.
        **kwargs  : dict. Passed to RZCrossSection.

        returns
        -------
        -         : obj. SurfaceBase type object.
        """        
        # Associate CrossSection
        self.name = RZCrossSection.name
        self.CrossSection = RZCrossSection
        
        # Toroidal extent
        self.extent = phi_extent
        
        # Fill out 2D surface in variety of coords
        if nphi%2: nphi+=1     # even number for fft
        self.pts = {}
        self.pts['r'] = self.CrossSection.pts['r'].repeat(nphi).reshape(-1,nphi)
        self.pts['z'] = self.CrossSection.pts['z'].repeat(nphi).reshape(-1,nphi)
        self.pts['theta'] = self.CrossSection.pts['theta'].repeat(nphi).reshape(-1,nphi)
        self.pts['theta (deg)'] = self.pts['theta']*180/np.pi - 90
        phi = np.linspace(phi_extent[0],phi_extent[1],nphi)
        self.pts['phi'] = np.array([phi]).repeat(len(self.CrossSection.pts['theta']),axis=0)
        self.pts['phi (deg)'] = self.pts['phi']*180/np.pi
        self.pts['x'] = self.pts['r']*np.cos(self.pts['phi'])
        self.pts['y'] = self.pts['r']*np.sin(self.pts['phi'])

    def mean(self,datakeys=None):
        """
        Mean value of data at diagnostic.

        kwargs
        ------
        - datakeys  : list. Any number of keys from data attribute.

        returns
        -------
        -           : dict. Means of each data key.
        """
        # housekeeping
        if not datakeys: datakeys=np.sort(self.data.keys()).tolist()
        if not type(datakeys) in (list,tuple): datakeys=(datakeys,)
        # calculation
        datamean = {}
        n = self.CrossSection.data_n
        for k in datakeys:
            if n==0: # no toroidal averaging
                datamean[k] = np.mean(self.CrossSection.data[k])
            else: # average toroidal mode
                datamean[k] = np.mean(np.real(self.data[k]
                                              *(np.exp(1j*n*self.extent[0])
                                                -np.exp(1j*n*self.extent[1]))/(1j*n)))
        return datamean
        

    def set_data(self,data,quiet=True,**kwargs):
        """
        Interpolate data from Data object to diagnostic points.

        args
        ----
        data     - obj. Data object

        **kwargs - dict. Passed to Data.interp
        """
        # get 1D data
        self.CrossSection.set_data(data,quiet=quiet,**kwargs)
        self.CrossSection.set_rzcharacteristics()
        
        # extent toroidally
        self.data = {}
        n = self.CrossSection.data_n
        keys = self.CrossSection.data.keys()
        for k in keys:
            self.data[k] = np.real(self.CrossSection.data[k].reshape(-1,1)
                                   *np.exp(1j*n*self.pts['phi']))

        return


######################################################## Form Diagnostics

# d3d RZCrossSections

_d3d_nodes = {'DIII-D Vessel Wall': [np.array([95.748,112.210,168.285,216.931,242.824,
                                     242.824,216.931,183.657,112.006,95.748,
                                     95.748,95.748,95.748])/100,
                           np.array([126.083,142.545,142.545,103.195,40.310,
                                     -39.916,-102.799,-142.342,-142.342,-126.083,
                                     -1e-10,1e-10,126.083])/100],

         'DIII-D Center Stack': [np.array([95.748,95.748,95.748,95.748])/100,
                                 np.array([-126.083,-1e-10,1e-10,126.083])/100],

         'DIII-D Low Field Side': [np.array([168.285,216.931,242.824,242.824,
                                             216.931,183.657])/100,
                                   np.array([142.545,103.195,40.310,-39.916,
                                             -102.799,-142.342])/100],

         'DIII-D Lower Baffle': [np.array([2.134,1.786,1.768,1.768,1.682,1.372]),
                                 np.array([-0.9730,-1.174,-1.211,-1.250,-1.250,-1.250])],

         'DIII-D Upper Baffle': [np.array([1.372,1.608,1.647,1.785,2.070,2.128]),
                                 np.array([1.292,1.095,1.077,1.077,1.040,0.9930])],

         'DIII-D First Wall': [np.array([1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.012,1.001,1.029,1.042,1.046,1.056,1.097,1.108,1.116,1.134,1.148,1.162,1.181,1.182,1.185,1.190,1.195,1.201,1.209,1.215,1.222,1.228,1.234,1.239,1.242,1.248,1.258,1.263,1.280,1.280,1.280,1.310,1.328,1.361,1.380,1.419,1.419,1.372,1.372,1.608,1.647,1.785,2.070,2.128,2.245,2.323,2.377,2.352,2.351,2.351,2.354,2.354,2.354,2.353,2.351,2.353,2.377,2.134,1.786,1.768,1.768,1.682,1.372,1.372,1.420,1.420,1.273,1.153,1.016,1.016,1.016,1.016,1.016,1.016,1.016,1.016]),
                               np.array([0.000E+00,0.9640,0.9680,1.001,1.019,1.077,1.070,1.096,1.113,1.138,1.147,1.165,1.217,1.217,1.1624,1.16238,1.1626,1.1645,1.16594,1.16591,1.16896,1.17175,1.17556,1.183,1.1835,1.185,1.188,1.191,1.196,1.202,1.208,1.214,1.221,1.231,1.238,1.244,1.254,1.278,1.290,1.331,1.347,1.348,1.348,1.348,1.348,1.348,1.348,1.310,1.310,1.292,1.095,1.077,1.077,1.040,0.9930,0.7090,0.5190,0.3890,0.4000,0.3370,0.2050,6.800E-02,0.000E+00,-6.800E-02,-0.2050,-0.3370,-0.4000,-0.3890,-0.9730,-1.174,-1.211,-1.250,-1.250,-1.250,-1.329,-1.329,-1.363,-1.363,-1.363,-1.223,-1.223,-0.8300,-0.8000,-0.4150,-0.4000,-1.000E-03,0.000E+00])]

         }

d3d_crosssections = {}
d3d_surfaces = {}
for k,v in _d3d_nodes.iteritems():
    close = 'Vessel' in k
    d3d_crosssections[k] = RZCrossSection(v[0],v[1],r0=1.69,z0=0,name=k,close=close)    
    d3d_surfaces[k] = AxisymSurface(d3d_crosssections[k])






######################################################## Input Visualization "Diagnostics"



def coils(coil,dim=3,cmap='RdBu_r',curlim=None,exclude=[],**kwargs):
    """
    Plot each coil referenced in the given coil input namelist.
    Coloring of the coils will be determined from its current normalized to +/-curlim.
    
    Colormaps can be seen at http://matplotlib.org/examples/color/colormaps_reference.html
    
    *Arguments:*
        - coil : str (obj).
            Coil input namelist file path or namelist.read dictionary object.
    
    *Key Word Arguments:*
        - dim : int.
            Choose from: 1) Plots R,Z cross-sections 2) Plots phi,theta projections 3) Plots 3D x,y,z using mayavi.
        - cmap : str.
            Valid matplotlib colormap key.
        - curlim : float.
            Normalizes all currents to range (-curlim,+curlim) for color mapping.
            Default is to use the maximum absolute current in the namelist.
        - exclude : list.
            Strings of coil names that are ignored (example: exclude=['il','iu']).
        
    Additional kwargs are passed to pyplot.plot or mayavi.mlab.plot3d functions.
    
    """
    # set up
    if type(coil)==str:
        coil = data.namelist.read(coil)
    f = kwargs.pop('figure',None)    
    axes = kwargs.pop('axes',None)
    if not f and not axes:
        if dim==3:
            f = data.mmlab.figure()#bgcolor=(1,1,1))
        else:
            f,axes = data.plt.subplots()
    if not axes and hasattr(f,'axes'):
        axes = np.array(f.axes).flat[0]
            
    colormap = data.plt.pyplot.get_cmap(cmap)
    to_rgb =  data.plt.matplotlib.colors.ColorConverter().to_rgb
    r0 = {'d3d':1.69,'nstx':1.00,'iter':6.42,'kstar':1.84,'rfxmod':0.46}
    
    # set normalization of currents using max absolute current
    if not curlim:
        curs = []
        for i in range(10):
            for j in range(50):
                if 'coil_cur({:},{:})'.format(i+1,j+1) in coil['COIL_CONTROL']:
                    curs.append(coil['COIL_CONTROL']['coil_cur({:},{:})'.format(i+1,j+1)])
        curlim = np.max(np.abs(curs)) #-np.min(curs)
    def current_rgb(cur):
        normcur = (cur+curlim)/(2*curlim) #*256
        rgb = to_rgb(colormap(normcur))
        return rgb
    if 'color' in kwargs:
        def current_rgb(cur):
            return kwargs['color']
    
    # get x,y,z and plot each coil
    coils = {}
    for i in range(1,10):
        if 'coil_name({:})'.format(i) in coil['COIL_CONTROL']:
            name = coil['COIL_CONTROL']['coil_name({:})'.format(i)]
            machine = coil['COIL_CONTROL']['machine']
            if name in exclude: continue
            cfile = '{:}_{:}.dat'.format(machine,name)
            with open(packagedir+'coil/'+cfile,'r') as f:
                line1 = f.readline()
                ncoil,s,nsec,nw = map(int,map(float,line1.split()))
            x,y,z = np.genfromtxt(packagedir+'coil/'+cfile,skip_header=1).T.reshape(3,ncoil,s,nsec)
            cur = [coil['COIL_CONTROL']['coil_cur({:},{:})'.format(i,j+1)] for j in range(ncoil)]
            # store info
            coils[name]={}
            coils[name]['x'] = x
            coils[name]['y'] = y
            coils[name]['z'] = z
            coils[name]['r']     = np.sqrt(x**2+y**2)
            coils[name]['phi']   = np.angle(x+1j*y)
            if machine in r0:
                coils[name]['theta'] = np.angle(coils[name]['r']-r0[machine]+1j*z)
            else:
                print('WARNING: {:} major radius unknown... no 2D plots available.'.format(machine))
            coils[name]['cur']   = cur
            coils[name]['ncoil'] = ncoil
            coils[name]['s']     = s
            coils[name]['nsec']  = nsec
            coils[name]['nw']    = nw
            #print(cur)
            #print([current_rgb(cur[j]) for j in range(ncoil)])
            # Plot
            for j in range(ncoil):
                if curlim: kwargs['color']=current_rgb(cur[j])
                if dim==3:
                    data.mmlab.plot3d(x[j].ravel(),y[j].ravel(),z[j].ravel(),#cur[j]*np.ones_like(z[j].ravel()),
                                      name=name+str(j),**kwargs)
                elif dim==2:
                    axes.plot(coils[name]['phi'][j].ravel()/dtor,coils[name]['theta'][j].ravel()/dtor,
                            label=name+str(j),**kwargs)
                else:
                    axes.plot(coils[name]['r'][j].ravel(),coils[name]['z'][j].ravel(),
                            label=name+str(j),**kwargs)
    
    return coils


############################################################ Magnetic Diagnostics

class _magnetics(object):
    """
    Class for getting synthetic magnetic diagnostics from IPEC results.
    """
    def __init__(self):
        """
        Just stores hardcoded geometries for now.
        Could eventually pass key words like machine='DIII-D' to make seperate
        classes for each machine.
        """
        self.toroidal_arrays = {
            # DIII-D
            'MPID66M' :{'rz':(2.413,  0.000),'length':1.40e-01,'angle':-np.pi/2},
            'MPID1A'  :{'rz':(0.977,  0.070),'length':1.40e-01,'angle': np.pi/2},  #note wall at r=0.95748
            'MPID1B'  :{'rz':(0.977, -0.070),'length':1.40e-01,'angle': np.pi/2},   #note wall at r=0.95748
            'MPID67A' :{'rz':(2.268,  0.759),'length':1.55e0-1,'angle':-1.17766},
            'MPID67B' :{'rz':(2.261, -0.753),'length':1.55e0-1,'angle':-1.96524},
            'MPID79A' :{'rz':(1.762,  1.312),'length':1.53e0-1,'angle':-0.68548},
            'MPID79B' :{'rz':(2.043, -1.110),'length':1.53e0-1,'angle':-2.25758},
        }
        
        # extend points along length
        for k,v in self.toroidal_arrays.iteritems():
            r,z = v['rz']
            l,a = v['length'],v['angle']
            rs = np.linspace(r-np.cos(a)*l/2,r+np.cos(a)*l/2,30)
            zs = np.linspace(z-np.sin(a)*l/2,z+np.sin(a)*l/2,30)
            v['rzs'] = np.array(zip(rs,zs))
    
    def get_signal(self,key,ds,op='mean'):
        """
        Get synthetic signal for array of magentic sensors.
        
        *Arguments:*
            - key : str. Sensor array from mags.
            - ds : Dataset. Dataset from ipec cylindrical output netcdf.
    
        *Key Word Arguments:*
            - op : str. Operation applied to the points along the surface
            spaned by the sensors (i.e. 'mean','max',etc.)
    
        *Returns:*
            - complex. Sensor array signal.
        
        """
        v = self.toroidal_arrays[key]
        
        # project over length of sensor cross section
        bz = interpn((ds['R'].values,ds['z'].values),ds['b_z-plas'].values.T,v['rzs'])
        br = interpn((ds['R'].values,ds['z'].values),ds['b_r-plas'].values.T,v['rzs'])
        bprp =-1*(br*np.sin(v['angle'])-bz*np.cos(v['angle']))
        bpar = 1*(br*np.cos(v['angle'])+bz*np.sin(v['angle']))
        if 'SL' in key:
            sigs = bprp
        else:
            sigs = bpar
        if op=='max':
            result = np.abs(sigs).max()
        else:
            result = np.mean(sigs)
        
        return result
magnetics = _magnetics()