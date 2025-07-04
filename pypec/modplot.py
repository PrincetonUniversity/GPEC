#!/usr/local/env python
"""
:mod:`pypec.modplot` -- matplotlib Customizations
=================================================

Collection of modified matplotlib functions and objects.

:Examples:

Plot complex arguments.

>>> f,ax = subplots()
>>> lines = ax.plot(np.arange(10)*(1+0.5j),label='test')
>>> pyplot.show()

Automatically resize plots to fit labels.

>>> xlbl = ax.set_xlabel('X AXIS')

Plot huge data sets quickly.

>>> x = np.linspace(0,9,1e5)
>>> data = np.arange(1e5)/1.5e4+(np.random.rand(1e5)-0.5)
>>> newline = ax.plot(x,data,label='large data set')

This plots a line capped at 1000-pts by default. The maximum number
of points is maintained as you manipulate the axis, so zooming
in will provide you with new points and increased detail until the
window samples fewer than that many points in the orginal data. The
first two lines, for instance, contain only their original 10 points
(not 1000 interpolated points).

"""

# standard matplotlib modules
import matplotlib  # standard plotting library
from matplotlib import pyplot
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable  # for splitting axes
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

# standard python modules
import sys, copy, os
import numpy as np  # math
from scipy.interpolate import interp1d  # math

# Want separate ID pyplot module for my hacks
# so this module does not change global matplotlib,pyplot,pylab,etc.
# Issue is import matplotlib imports many modules... all intertwined.
# METHOD 1
# import imp #built in method for assigning IDs
# mfile,mpath,mdesc = imp.find_module('matplotlib')
# pyplot = imp.load_module('modpyplot', *imp.find_module('pyplot',[mpath]))
# Method 2
# save = sys.modules.pop('matplotlib.pyplot', None)
# import matplotlib.pyplot as pyplot
# sys.modules['matplotlib.pyplot'] = save

# matplotlib modifications
try:
    interactive = matplotlib.is_interactive()
    import seaborn

    seaborn.reset_orig()
    matplotlib.interactive(interactive)
except:
    seaborn = None
# override user preferences from ~/.matplotlib/matplotlibrc file.
# really should define new functions
# matplotlib.rcParams['figure.autolayout']=True # Necessary for tight_layout
# matplotlib.rcParams['legend.frameon']=False
# matplotlib.rcParams['legend.loc']='best'

linestyles = ["-", "--", "-.", ":", "None", " ", ""]
show = pyplot.show
draw = pyplot.draw


########################################### some figure sizes (journal defualts)

rcp_size = np.array(matplotlib.rcParams["figure.figsize"])
pop_size = np.array([3.346, 3.346])  # single column. Do *2 for double

########################################### default colormaps

# this got the new matplotlib colormaps before their were made standard
from . import colormaps as cmaps

if not hasattr(matplotlib.cm, "viridis"):
    for k in ["magma", "inferno", "plasma", "viridis"]:
        pyplot.register_cmap(name=k, cmap=cmaps.cmaps[k])
        pyplot.register_cmap(name=k + "_r", cmap=cmaps.cmaps[k + "_r"])
    f, ax = pyplot.subplots()
    pyplot.set_cmap(cmaps.cmaps["viridis"])
    pyplot.close(f)

########################################### default colormaps


def png_to_gif(files, gif_file, delay=20, loop=0, clean=False):
    """gif_file should end in .gif"""
    gif_file = gif_file.rstrip(".gif")
    os.system(
        "convert -delay {d} -loop {l} {i} {o}".format(
            d=delay, l=loop, i=" ".join(files), o=gif_file + ".gif"
        )
    )
    if clean:
        os.system(
            "zip -j {zipfile} {files}".format(
                zipfile=gif_file + ".zip", files=" ".join(files)
            )
        )
        os.system("rm {files}".format(files=" ".join(files)))

    return


########################################### Custom Styles


def set_style(style=None, rc={}):
    """
    Seaborn set_style with additional 'thick',

    Thick (thin) style multiplies rcParams axes.linewidth, lines.linewidth,
    xtick.major.width, xtick.minor.width, ytick.major.width, and
    ytick.minor.width by 2 (0.5).
    """
    if style in ["thick", "thin"]:
        if pyplot.rcParams["font.weight"] == "normal":
            pyplot.rcParams["font.weight"] = "400"
        elif pyplot.rcParams["font.weight"] == "bold":
            pyplot.rcParams["font.weight"] = "700"
        if style == "thick":
            fac = 3 / 2.0
            pyplot.rcParams["font.weight"] = str(
                min((900, int(pyplot.rcParams["font.weight"]) + 200))
            )
        elif style == "thin":
            fac = 2.0 / 3
            pyplot.rcParams["font.weight"] = str(
                max((100, int(pyplot.rcParams["font.weight"]) - 200))
            )

        for k, v in pyplot.rcParams.items():
            if "weight" in k:
                pyplot.rcParams[k] = pyplot.rcParams["font.weight"]

        for k, v in dict(rc).items():
            pyplot.rcParams[k] = v

        for k in [
            "axes.linewidth",
            "lines.linewidth",
            "xtick.major.width",
            "xtick.minor.width",
            "ytick.major.width",
            "ytick.minor.width",
        ]:
            pyplot.rcParams[k] *= fac

    else:
        seaborn.set_style(style, rc)

    return


########################################### modified pyplot functions


def subplots(
    nrows=1,
    ncols=1,
    sharex=False,
    sharey=False,
    squeeze=True,
    autosize=True,
    subplot_kw=None,
    powerlim=(-3, 3),
    useOffset=False,
    nybins=None,
    nxbins=None,
    **fig_kw
):
    """
    Matplotlib subplots with default figsize = rcp_size*[ncols,nrows]
    fig_kw.setdefault('figsize',figsize).

    Additional Key Word Arguments:
      powlim    : tuple.
        Axis labels use scientific notion above 10^power.
      useOffset : bool.
        Axis labels use offset if range<<average.

    Accepts standargs args and kwargs for pyplot.subplots.


    """
    figsize = np.array(matplotlib.rcParams["figure.figsize"]) * [ncols, nrows]
    if autosize:
        fig_kw.setdefault("figsize", figsize)
    f, ax = pyplot.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=sharex,
        sharey=sharey,
        squeeze=squeeze,
        subplot_kw=subplot_kw,
        **fig_kw
    )
    # customize axes and figure
    if nybins is not None:
        for a in ax.flat:
            a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nybins))
    if nxbins is not None:
        for a in ax.flat:
            a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nxbins))
    f = _modfigure(f)

    return f, ax


subplots.__doc__ += "ORIGINAL DOCUMENTATION \n\n" + str(pyplot.subplots.__doc__)


def figure(
    num=None,
    figsize=None,
    dpi=None,
    facecolor=None,
    edgecolor=None,
    frameon=True,
    FigureClass=matplotlib.figure.Figure,
    **kwargs
):
    f = pyplot.figure(
        num=num,
        figsize=figsize,
        dpi=dpi,
        facecolor=facecolor,
        edgecolor=edgecolor,
        frameon=frameon,
        FigureClass=FigureClass,
        **kwargs
    )
    f = _modfigure(f)
    return f


def colorbar(mappable=None, cax=None, ax=None, use_gridspec=True, **kw):
    """
    Modified pyplot colorbar for default use_gridspec=True.

    """
    cb = pyplot.colorbar(
        mappable=mappable, cax=cax, ax=ax, use_gridspec=use_gridspec, **kw
    )
    return cb


colorbar.__doc__ += "ORIGINAL DOCUMENTATION \n\n" + str(pyplot.colorbar.__doc__)


########################################### Figure object modifications


# internal function to apply customizations post initialization
def _modfigure(f):
    """
    Modplot figure object modifications. These include:
        1) canvas connection to modplot.key_press event funtion

    Arguments:
      f : Figure object.

    """
    f.canvas.mpl_connect("key_press_event", onkey)  # custom hotkeys
    return f


# customize axes when added
def _figure_add_subplot(self, *args, **kwargs):
    ax = self._orig_add_subplot(*args, **kwargs)
    ax = _modaxes(ax)
    return ax


if not hasattr(pyplot.Figure, "_orig_add_subplot"):  # prevent recursion when reloaded
    pyplot.Figure._orig_add_subplot = copy.deepcopy(pyplot.Figure.add_subplot)
pyplot.Figure.add_subplot = _figure_add_subplot


# tight_layout when show
def _figure_show(self, *args):
    self.tight_layout()
    self._orig_show(*args)


_figure_show.__doc__ = pyplot.show.__doc__
# pppl Figure base doesn't have show from nomachine or 2.7.2 for some reason
if not hasattr(pyplot.Figure, "_orig_show"):  # prevent recursion when reloaded
    if sys.version_info < (2, 7, 3):
        pyplot.Figure._orig_show = copy.deepcopy(pyplot.show)
    else:
        pyplot.Figure._orig_show = copy.deepcopy(pyplot.Figure.show)
pyplot.Figure.show = _figure_show


# colorbars are axes (for tight_layout)
def _figure_colorbar(self, mappable, cax=None, ax=None, use_gridspec=True, **kw):
    cb = self._orig_colorbar(mappable, cax=cax, ax=ax, use_gridspec=use_gridspec, **kw)
    return cb


_figure_colorbar.__doc__ = "Modified default use_gridspec=True. \n - \n" + str(
    pyplot.colorbar.__doc__
)
if not hasattr(pyplot.Figure, "_orig_colorbar"):  # prevent recursion when reloaded
    pyplot.Figure._orig_colorbar = copy.deepcopy(pyplot.Figure.colorbar)
pyplot.Figure.colorbar = _figure_colorbar


# new method to print text file of lines in figure
def printlines(self, filename, labeled_only=False, squeeze=False):
    """
    Print all data in line plot(s) to text file.The x values
    will be taken from the line with the greatest number of
    points in the (first) axis, and other lines are interpolated
    if their x values do not match. Column lables are the
    line labels and xlabel.

    Arguments:
      filename : str.
        Path to print to.

    Key Word Arguments:
      labeled_only : bool
        Only print labeled lines to file.
      squeeze : bool.
        Print lines from all axes in a single table.

    returns
      bool.
        True.

    """
    with open(filename, "w") as f:
        data = []

        # Try to avoid clutter in squeeze
        samettle = True
        sameylbl = True
        if squeeze:
            ylbl = self.get_axes()[0].get_ylabel()
            ttle = self.get_axes()[0].get_title()
            for a in self.get_axes():
                if a.get_ylabel() != ylbl:
                    sameylbl = False
                if a.get_title() != ttle:
                    samettle = False
            if samettle:
                f.write(ttle + "\n\n")
            if sameylbl:
                f.write(
                    "  ylabel = "
                    + ylbl.encode().translate(str.maketrans(dict.fromkeys(" \${}")))
                    + "\n\n"
                )

        for a in self.get_axes():
            if not a.lines:
                continue
            if not squeeze:
                data = []
                f.write("\n" + a.get_title() + "\n\n")
                f.write(
                    "  ylabel = "
                    + a.get_ylabel()
                    .encode()
                    .translate(str.maketrans(dict.fromkeys(" \${}")))
                    + "\n\n"
                )
            # use x-axis from line with greatest number of pts
            xs = []
            for line in a.lines:
                label = (
                    line.get_label()
                    .encode()
                    .translate(str.maketrans(dict.fromkeys(" \${}")))
                )
                if labeled_only and (label.startswith("_") or not label):
                    continue
                xs.append(line.get_xdata())
            longest = np.array([len(x) for x in xs]).argmax()
            if not data:
                data.append(xs[longest])
                f.write(
                    "{0:>25s}".format(
                        a.get_xlabel()
                        .encode()
                        .translate(str.maketrans(dict.fromkeys(" \${}")))
                    )
                )
            # label and extrapolate each line
            for line in a.lines:
                label = (
                    line.get_label()
                    .encode()
                    .translate(str.maketrans(dict.fromkeys(" \${}")))
                )
                if labeled_only and (label.startswith("_") or not label):
                    continue
                if squeeze and not sameylbl:
                    label = (
                        a.get_ylabel()
                        .encode()
                        .translate(str.maketrans(dict.fromkeys(" \${}")))
                        + label
                    )
                f.write("{0:>25s}".format(label))
                # standard axis
                x, y = line.get_xdata(), line.get_ydata()
                if np.any(x != data[0]):
                    fit = interp1d(x, y, bounds_error=False, fill_value=np.nan)
                    data.append(fit(data[0]))
                else:
                    data.append(y)
            if not squeeze:
                f.write("\n ")
                if np.any(np.iscomplex(data)):
                    print("WARNING: Complex data may be lost.")
                data = np.array(data).real.T
                for row in data:
                    row.tofile(f, sep=" ", format="%24.9E")
                    f.write("\n ")
                f.write("\n")
        if squeeze:
            f.write("\n ")
            if np.any(np.iscomplex(data)):
                print("WARNING: Complex data may be lost.")
            data = np.array(data).real.T
            for row in data:
                row.tofile(f, sep=" ", format="%24.9E")
                f.write("\n ")
            f.write("\n")
    print("Wrote lines to " + filename)
    return True


pyplot.Figure.printlines = printlines


########################################### Axes object modifications


# internal function to apply customizations post initialization
def _modaxes(ax, useOffset=False):
    """
    Modplot axis object modifications. These include:
        1) Optional offset
        2) Ticker format power lims (-3,3)
        3) New downsample instance, connected to x-axis rescale
        4) New plot instance with complex arg handling

    Arguments:
      ax        : obj.
        Axes object.
    Key Word Arguments:
      useOffset : bool.
        Offset axis zero.

    returns
      Axis object.

    """
    afmt = matplotlib.ticker.ScalarFormatter(useOffset=useOffset)
    afmt.set_powerlimits((-3, 3))
    ax.callbacks.connect("xlim_changed", ax.downsample)
    ax.yaxis.set_major_formatter(afmt)
    ax.xaxis.set_major_formatter(afmt)
    return ax


# down sample all lines in axes
def _axes_downsample(self, ax, npts=10000):
    """
    Attempt to down sample original data in
    each Line2D object in lines.

    Arguments:
      ax   : obj.
        Initialized Axes object containing lines
    Key Word Arguments:
      npts : int.
        Number of points plotted in each line

    """
    for l in ax.lines:
        try:
            l.downsample(npts=npts)
        except:
            pass
    # draw()


pyplot.Axes.downsample = _axes_downsample


# complex value plotting
def _axes_plot(self, *args, **kwargs):
    """
    Modified axes line plotting, with suport for complex plotting.
    Original matplotlib Axes plot method now located at _orig_plot.
    All original plot args and kwargs accepted.

    """
    linestyles = matplotlib.lines.lineStyles.keys()
    newargs = map(_real_arg, args)
    lines1 = self._orig_plot(*newargs, **kwargs)
    lines2 = []
    if any([any(np.ravel(np.iscomplex(arg))) for arg in args]):
        newargs = [args[0]]
        if len(args) > 1:
            newargs += map(_imag_arg, args[1:])
        lines2 = self._orig_plot(*newargs, **kwargs)
        for l1, l2 in zip(lines1, lines2):
            l2.set_color(l1.get_color())
            ils = linestyles.index(l1.get_linestyle())
            l2.set_linestyle(linestyles[ils + 1])
            l2.set_label(r"$\Im$ " + l1.get_label())
            l1.set_label(r"$\Re$ " + l1.get_label())
        pyplot.gca().legend(ncol=max(1, len(self.lines) / 6))
    # draw()
    return lines1 + lines2


_axes_plot.__doc__ += str(pyplot.plot.__doc__)
if not hasattr(pyplot.Axes, "_orig_plot"):  # prevent recursion when reloaded
    pyplot.Axes._orig_plot = copy.deepcopy(pyplot.Axes.plot)
pyplot.Axes.plot = _axes_plot


def _imag_arg(arg):
    # if type(arg)==str:
    #    return arg
    # else:
    if np.any(np.iscomplex(arg)):
        return np.imag(arg)
    else:
        return arg


def _real_arg(arg):
    # if type(arg)==str:
    #    return arg
    # else:
    if np.any(np.iscomplex(arg)):
        return np.real(arg)
    else:
        return arg


def _center_axes(ax):
    ax.spines["left"].set_position("center")
    ax.spines["left"].set_linewidth(2)
    ax.spines["right"].set_color("none")
    ax.spines["bottom"].set_position("center")
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_smart_bounds(True)
    ax.spines["bottom"].set_smart_bounds(True)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")


pyplot.Axes.center_axes = _center_axes


########################################### modified pyplot functions


# downsample data from the full set used to initialize the line
def line_downsample(line, npts=1000):
    """
    Downsample line data for increased plotting speed.

    Key Word Arguments:
      npts : int.
        Number of data points within axes xlims.
    """
    lims = line.axes.get_xlim()
    if not hasattr(line, "_xinit"):
        line._xinit = line.get_xdata()
        line._yinit = line.get_ydata()
    window = (line._xinit >= lims[0]) & (line._xinit <= lims[1])
    # extend out edge
    window[window.argmax() - 1] = True
    window[window.argmax() + window[window.argmax() :].argmin()] = True
    step = int(line._xinit[window].size / npts) + 1
    line.set_xdata(line._xinit[window][::step])
    line.set_ydata(line._yinit[window][::step])


pyplot.matplotlib.lines.Line2D.downsample = line_downsample


########################################### my own functions


def add_inset_axes(ax, rect, axisbg="w", polar=False, **kwargs):
    """
    Add a axes inset in another axes.
    Unlike the native matplotlib inset_axes, here the inset axes can be a polar axes.

    Note that the parent axes should be drawn first if using autolayout (changes postion).
    The inset axes will not track the original if it is moved.

    From http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib/17479417#17479417

    :param ax: Axes. Parent axes.
    :param rect: list. Inset axes location and dimensions [left, bottom, width, height]
    :param axisbg: str. Background color.
    :param polar: bool. Polar plot.
    :param kwargs: dict. Passed to add_subplot.
    :return: Axes. A new inset axes.

    """
    fig = ax.figure
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x, y, width, height], axisbg=axisbg, polar=polar, **kwargs)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2] ** 0.5
    y_labelsize *= rect[3] ** 0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name="shiftedcmap"):
    """
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Taken from http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Args:
        cmap: The matplotlib colormap to be altered
        start: Offset from lowest point in the colormap's range. Should be between 0 and 'midpoint'.
        midpoint: The new center of the colormap. Should be between 0.0 and 1.0. Recommended is  1 - vmax/(vmax + abs(vmin)).
        stop: Offset from highest point in the colormap's range. Should be between 'midpoint' and 1.0.
        name: str. Added to the name of cmap when creating the new colormap.

    Example:

    The following creates a shifted colormap that can be used in matplotlib arguments,

        >>> orig_cmap = matplotlib.cm.coolwarm
        >>> shifted_cmap = shiftedColorMap(orig_cmap, midpoint=0.75, name='shifted')

    """
    cdict = {"red": [], "green": [], "blue": [], "alpha": []}

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack(
        [
            np.linspace(0.0, midpoint, 128, endpoint=False),
            np.linspace(midpoint, 1.0, 129, endpoint=True),
        ]
    )

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict["red"].append((si, r, r))
        cdict["green"].append((si, g, g))
        cdict["blue"].append((si, b, b))
        cdict["alpha"].append((si, a, a))
    name = "{:}_{:}".format(cmap.name, name)
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    pyplot.register_cmap(cmap=newcmap)

    return newcmap


class MidPointNorm(Normalize):
    """
    Taken from http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    """

    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        """
        A Normalize subclass that, when used, will stretch a colormap out such that the
        dynamic range has (vmin, midpoint, vmax) even if midpoint is not half way between the ends.

        Note, a side-product of this approach is that the colorbar labeling is similarly stretched.

        Args:
            midpoint: numeric. Sets the data value that will be the center of the colormap dynamic range.
            vmin: numeric. Sets the data value that will be the min of the colormap dynamic range.
            vmax: numeric. Sets the data value that will be the max of the colormap dynamic range.
            clip: numeric.

        Returns:
            class. To be used as a norm with matplotlib.

        Example:
            >>> norm = MidPointNorm(midpoint=0.8)
            >>> f,ax = subplots()
            >>> im = ax.imshow(np.random.random((10,10)), norm=norm, cmap='seismic')
            >>> f.colorbar(im,ax=ax)

        """
        Normalize.__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0)  # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = np.ma.getmask(result)
                result = np.ma.array(
                    np.clip(result.filled(vmax), vmin, vmax), mask=mask
                )

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            # First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat > 0] /= abs(vmax - midpoint)
            resdat[resdat < 0] /= abs(vmin - midpoint)

            resdat /= 2.0
            resdat += 0.5
            result = np.ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if matplotlib.cbook.iterable(value):
            val = np.ma.asarray(value)
            val = 2 * (val - 0.5)
            val[val > 0] *= abs(vmax - midpoint)
            val[val < 0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0:
                return val * abs(vmin - midpoint) + midpoint
            else:
                return val * abs(vmax - midpoint) + midpoint


def plot_axes(ax, fig=None, geometry=(1, 1, 1)):
    """
    Re-create a given axis in a new figure. This allows, for instance,
    a subplot to be moved to its own figure where it can be manipulated
    and/or saved independent of the original.

    Arguments:
      ax : obj.
        An initialized Axes object
    Key Word Arguments:
      fig: obj.
        A figure in which to re-create the axis
      geometry : tuple.
        Axes geometry of re-created axis

    returns
      Figure.

    """
    if fig is None:
        fig = figure()
    a2 = copy.copy(ax)
    a2.set_figure(fig)
    a2 = fig.add_axes(a2)
    #    for k,v in vars(ax).items():
    #        if '_' not in k:
    #            try:
    #                setattr(ax2,k,v)
    #            except:
    #                pass
    if ax.get_geometry() != geometry:
        a2.change_geometry(*geometry)
    #    draw()
    fig.show()
    return fig


def xyaxes(axes, x=True, y=True, **kwargs):
    """
    Plot hline and vline through (0,0).

    **Arguments:**
        - axes : obj. Matplotlib Axes object.

    **key Word Arguments:**
        - x : bool. Draw hline
        - y : bool. Draw vline

    **Returns:**
        - bool.

    All other kwargs passed to hlines and vlines functions.

    """
    if y:
        ylim = axes.get_ylim()
        axes.vlines(0, *ylim, **kwargs)
        axes.set_ylim(*ylim)
    if x:
        xlim = axes.get_xlim()
        axes.hlines(0, *xlim, **kwargs)
        axes.set_xlim(*xlim)
    return True


def set_linearray(lines, values=None, cmap=None, vmin=None, vmax=None):
    """
    Set colors of lines to colormaping of values.

    Other good sequential colormaps are YlOrBr and autumn.
    A good diverging colormap is bwr.

    **Arguments:**
        - lines : list. Lines to get set colors.

    **Key Word Arguments:**
        - values : array like. Values corresponding to each line. Default is indexing.
        - cmap : str. Valid matplotlib colormap name.
        - vmax : float. Upper bound of colormaping.
        - vmin : float. Lower bound of colormaping.

    **Returns:**
        - ScalarMappable. A mapping object used for colorbars.

    """
    if values == None:
        values = range(len(lines))
    nm = matplotlib.colors.Normalize(vmin, vmax)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=nm)
    sm.set_array(values)
    colors = sm.cmap(sm.norm(values))
    for l, c in zip(lines, colors):
        l.set_color(c)

    return sm


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
    if event.key == "n":
        if event.inaxes:
            f = plot_axes(event.inaxes)
    if event.key == "t":
        f = pyplot.gcf()
        f.tight_layout()
        f.show()
    if event.key == "c":
        event.inaxes._center_axis()
    return
