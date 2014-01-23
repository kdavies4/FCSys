#!/usr/bin/env python
"""Load and analyze results from FCSys_.

.. _FCSys: http://kdavies4.github.com/FCSys/
"""
__author__ = "Kevin Davies"
__email__ = "kdavies4@gmail.com"
__copyright__ = "Copyright 2012-2013, Georgia Tech Research Corporation"
__license__ = "BSD-compatible (see LICENSE.txt)"


import os
import numpy as np
import matplotlib.pyplot as plt

from collections import namedtuple
from matplotlib import rcParams
from matplotlib.cbook import iterable
from modelicares import (figure, unit2tex, label_number, label_quantity,
                         setup_subplots)
from modelicares.base import flatten_list, get_pow10, convert
from simres import SimRes

class FCSimRes(SimRes):
    """Fuel cell simulation results from FCSys_ and methods to analyze those
    results
    """

    # Global constants
    LAYERS = ['anFP', 'anGDL', 'anCL', 'PEM', 'caCL', 'caGDL', 'caFP']
    LAYER_INFO = {'anFP': "Anode flow plate",
                  'anGDL': "Anode GDL",
                  'anCL': "Anode catalyst layer",
                  'PEM': "PEM",
                  'caCL': "Cathode catalyst layer",
                  'caGDL': "Cathode GDL",
                  'caFP': "Cathode flow plate"}

    def animate_quiverfig_current(self, times, fname='current_eminus_xy'):
        """Create an animation of the x-direction current in the x-y plane.

        **Arguments:**

        - *times*: Times at which the frames should be generated

        - *fname*: Filename for the movie

            ".mpg" will be appended if necessary.
        """
        from res import animate, saveall
        # TODO 10/28/11: Use the new movie utilities in matplotlib.

        # Set up.
        #[start_time, stop_time] = self.get_times('Time', [0,-1])
        #if not times:
        #    times = self.get_values
        #movie_dir = os.path.join(self.dir, movie_dir)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        # Plot and animate.
        for i, time in enumerate(times):
            figure(os.path.join(output_dir, '_tmp%03d' % i))
            #self.currentx_at_times(time=time)
            plt.plot(time, np.sin(time), '*') # Temporary; for debug
        saveall('png')
        animate(fname)

# Example to add to documentation somewhere:
     #   **Example:**
#
      #     >>> from modelicares import saveall
      #     >>> from fcres import FCSimRes
#
     #      >>> sim = FCSimRes('examples/Polarization.mat')
     #      >>> sim.plotfig(xname='cell.I',
     #      ...             ynames1='cell.v', ylabel1="Potential",
     #      ...             legends1="Average voltage",
    #       ...             ynames2='cell.Wdot', ylabel2="Power",
    #       ...             legends2="Power output",
    #       ...             title="Cell Polarization",
    #       ...             label='examples/Polarization') # doctest: +ELLIPSIS
    #       (<matplotlib.axes.AxesSubplot object at 0x...>, <matplotlib.axes.AxesSubplot object at 0x...>)
    #       >>> aveall()
     #      Saved examples/Polarization.pdf
     #      Saved examples/Polarization.png
#
     #   .. only:: html
#
     #      .. image:: examples/Polarization.png
     #         :scale: 50 %
     #         :alt: plot of cell profile

     #   .. only:: latex
#
      #     .. figure:: examples/Polarization.pdf
      #        :scale: 70 %"""



    def colorfig(self, suffix, times=[0], n_rows=1, title="", subtitles=[],
                 label="color", slice_axis='z', slice_index=0,
                 xlabel="", xticklabels=[], xticks=[],
                 ylabel="", yticklabels=[], yticks=[],
                 clabel="", cbar_orientation='vertical',
                 margin_left=rcParams['figure.subplot.left'],
                 margin_right=1-rcParams['figure.subplot.right'],
                 margin_bottom=rcParams['figure.subplot.bottom'],
                 margin_top=1-rcParams['figure.subplot.top'],
                 margin_cbar=0.2,
                 wspace=0.1, hspace=0.25,
                 cbar_space=0.1, cbar_width=0.05,
                 **kwargs):
        """Create a figure with 2D scalar data at given time(s) on a color axis
        in 2D Cartesian coordinates.

        **Arguments:**

        - *suffix*: Name of the variable to be plotted (relative to the names
          of the subregions)

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *n_rows*: Number of rows of (sub)plots

        - *title*: Title for the figure

        - *subtitles*: List of subtitles (i.e., titles for each subplot)

             If not provided, "t = xx s" will be used, where xx is the time of
             each entry.  "(initial)" or "(final)" are appended if appropriate.

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *slice_axis*: Axis normal to the screen or page ('x', 'y', or 'z')

        - *slice_index*: Position along slice_axis

        - *xlabel*: Label for the x-axes (only shown for the subplots on the
          bottom row)

        - *xticklabels*: Labels for the x-axis ticks (only shown for the
          subplots on the bottom row)

        - *xticks*: Positions of the x-axis ticks

        - *ylabel*: Label for the y axis (only shown for the subplots on the
          left column)

        - *yticklabels*: Labels for the y-axis ticks (only shown for the
          subplots on the left column)

        - *yticks*: Positions of the y-axis ticks

        - *clabel*: Label for the color- or c-bar axis

        - *cbar_orientation*: Orientation of the colorbar ("vertical" or
          "horizontal")

        - *margin_left*: Left margin

        - *margin_right*: Right margin (ignored if ``cbar_orientation ==
          'vertical'``)

        - *margin_bottom*: Bottom margin (ignored if ``cbar_orientation ==
          'horizontal'``)

        - *margin_top*: Top margin

        - *margin_cbar*: Margin reserved for the colorbar (right margin if
          ``cbar_orientation == 'vertical'`` and bottom margin if
          ``cbar_orientation == 'horizontal'``)

        - *wspace*: The amount of width reserved for blank space between
          subplots

        - *hspace*: The amount of height reserved for white space between
          subplots

        - *cbar_space*: Space between the subplot rectangles and the colorbar

        - *cbar_width*: Width of the colorbar if vertical (or height if
          horizontal)

        - *\*\*kwargs*: Additional arguments for  :meth:`modelicares.base.color`
        """
        # The procedure for coordinating the images with a single colorbar was
        # copied and modified from
        # http://matplotlib.sourceforge.net/examples/pylab_examples/multi_image.html,
        # accessed 11/8/10.

        # 5/26/11: TODO: This method needs to be updated. Use quiverfig() as a
        # reference.

        from res import color

        # Get the data.
        n_plots = len(times)
        if slice_axis == 'x':
            names = presuffix(self.subregions[slice_index], suffix=suffix)
        elif slice_axis == 'y':
            names = presuffix(self.subregions[:][slice_index], suffix=suffix)
        else:
            names = presuffix(self.subregions[:][:][slice_index],
                              suffix=suffix)
        c = self.get_values_at_times(names, times)

        # Cast the data into a list of matrices (one at each sample time).
        if slice_axis == 'x':
            c = [np.array([[c[i_y + self.n_y*i_z][i]
                 for i_y in range(self.n_y-1,-1,-1)]
                 for i_z in range(self.n_z)]) for i in range(n_plots)]
        elif slice_axis == 'y':
            c = [np.array([[c[i_z + self.n_z*i_x][i]
                 for i_z in range(self.n_z-1,-1,-1)]
                 for i_x in range(self.n_x)]) for i in range(n_plots)]
        else:
            c = [np.array([[c[i_y + self.n_y*i_x][i]
                 for i_z in range(self.n_y-1,-1,-1)]
                 for i_x in range(self.n_x)]) for i in range(n_plots)]
        [start_time, stop_time] = self.get_times('Time', [0,-1])

        # Generate xlabel, xticks, and xticklabels.
        if not xlabel:
            if slice_axis == 'x':
                xlabel = 'z-axis index'
            elif slice_axis == 'y':
                xlabel = 'x-axis index'
            elif slice_axis == 'z':
                xlabel = 'y-axis index'
        if not xticklabels:
            if slice_axis == 'x':
                xticklabels = [str(i+1) for i in range(self.n_z)]
            elif slice_axis == 'y':
                xticklabels = [str(i+1) for i in range(self.n_x)]
            else:
                xticklabels = [str(i+1) for i in range(self.n_x)]
        if not xticks:
            if slice_axis == 'x':
                xticks = range(self.n_z)
            elif slice_axis == 'y':
                xticks = range(self.n_x)
            else:
                xticks = range(self.n_x)

        # Generate ylabel, yticks, and yticklabels.
        if not ylabel:
            if slice_axis == 'x':
                ylabel = 'y-axis index'
            elif slice_axis == 'y':
                ylabel = 'z-axis index'
            else:
                ylabel = 'x-axis index'
        if not yticklabels:
            if slice_axis == 'x':
                yticklabels = [str(i+1) for i in range(self.n_y)]
            elif slice_axis == 'y':
                yticklabels = [str(i+1) for i in range(self.n_z)]
            else:
                yticklabels = [str(i+1) for i in range(self.n_y)]
        if not yticks:
            if slice_axis == 'x':
                yticks = range(self.n_y)
            elif slice_axis == 'y':
                yticks = range(self.n_z)
            else:
                yticks = range(self.n_y)

        # Generate clabel.
        if not clabel:
            clabels = self.get_description(names)
            # If all of the descriptions are the same, use the first one.
            if len(set(clabels)) == 1:
                clabel = clabels[0]
            else:
                clabel = "Value"
        #units = self.get_unit(names)
        units = self.get_unit(names)
        if len(set(units)) == 1:
            clabel += label_number("", units[0])
        else:
            raise UserWarning("The variables have inconsistent units.  The "
                "colorbar unit will not match the units of all of the "
                "variables.")

        # Set up the subplots.
        if not subtitles:
            #unit = unit2tex(self.get_unit('Time'))
            unit = unit2tex(self.get_unit('Time'))
            subtitles = ["t = " + label_quantity(times[i], unit)
                for i in n_plots]
            for i, time in enumerate(times):
                if time == start_time:
                    subtitle += " (initial)"
                elif time == stop_time:
                    subtitle += " (final)"
        ax, cax = setup_subplots(n_plots=n_plots, n_rows=n_rows,
            title=title, subtitles=subtitles, label=label,
            xlabel=xlabel, xticklabels=xticklabels, xticks=xticks,
            ylabel=ylabel, yticklabels=yticklabels, yticks=yticks,
            ctype=cbar_orientation, clabel=clabel,
            margin_left=margin_left, margin_right=margin_right,
            margin_bottom=margin_bottom, margin_top=margin_top,
            margin_cbar=margin_cbar, wspace=wspace, hspace=hspace,
            cbar_space=cbar_space, cbar_width=cbar_width)

        # Create the plots.
        profiles = []
        c_min = np.inf
        c_max = -np.inf
        for i, time in enumerate(times):
            profiles.append(color(ax[i], c[i], **kwargs))

            # Find the minimum and maximum of the color-axis data.
            c_min = min(c_min, np.amin(c[i]))
            c_max = max(c_max, np.amax(c[i]))

        # Set the first image as the master, with all the others observing it
        # for changes in norm.
        class ImageFollower:
            """Update image in response to changes in clim on another image.
            """
            def __init__(self, follower):
                self.follower = follower
            def __call__(self, leader):
                self.follower.set_clim(leader.get_clim())
        norm = Normalize(c_min=c_min, c_max=c_max)
        for i, profile in enumerate(profiles):
            profile.set_norm(norm)
            if i > 0:
                profiles[0].callbacksSM.connect('changed',
                                                ImageFollower(profile))

        # Add the colorbar. (It is also based on the master image.)
        cbar = fig.colorbar(images[0], cax, orientation=cbar_orientation)

        # Scale the colorbar.
        if c_min == c_max:
            yticks = [c_min]
            for i in range(len(cbar.ax.get_yticks()) - 1):
                yticks.append(None)
            cbar.set_ticks(yticks)

    def colorfig_pressure(self, title="Pressure Profile",
                          label='pressureprofile', times=None):
        """Plot the pressure profiles of the subregions.

        **Arguments:**

        - *title*: Title for the figure

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.
        """
        if times is None:
            [start_time, stop_time] = self.get_values('Time', [0,-1])
            Delta_time = (stop_time - start_time)/5
            times = np.arange(start_time, stop_time + Delta_time/2, Delta_time)
        self.gen_profiles_fig(fname=fname, title=title,
            names=presuffix(self.storage_names, suffix='.p'),
            n_y=self.N_storagey, times=times,
            n_rows=int(np.ceil(np.sqrt(len(times)))),
            xlabel='', xticklabels=['AnFP','AnCL','CaCL','CaFP'],
            xticks=range(self.N_storagex),
            n_x_layers=self.N_storagex_layers,
            ylabel="y index", yticklabels=range(self.N_storagey-1,-1,-1),
            yticks=range(self.N_storagey), clabel=None)

        # TODO: Fix these later and make the code in gen_profiles_fig() more
        # general.
        #xticks = range(self.N_storagex)
        #xticklabels = ['AnFP','AnCL','CaCL','CaFP']

        # TODO: Create similar methods for temperature, density, etc.

    def currdenfig(self, title=None, times=[0], z_index=0, leg_kwargs={},
                   **kwargs):
        """Plot current densities of the cell segments at times.

        **Arguments:**

        - *title*: Title for the figure

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *z_index*: z-index at which the current densities are taken

        - *leg_kwargs*: Dictionary of keyword arguments for
          :meth:`matplotlib.pyplot.legend`

             If *leg_kwargs* is *None*, then no legend will be shown.

        - *\*\*kwargs*: Additional arguments for  :meth:`barfig`
        """
        # Process the arguments.
        xlabel = kwargs.pop('xlabel', 'y-axis index (inlet to outlet)')
        name_template = self.cell + ('iprimeprime_seg[%i, ' + '%i]'
            % (z_index+1))
        #description = self.get_description(name_template % 1)
        #unit = self.get_unit(name_template % 1)
        unit = self.get_unit(name_template % 1)
        #ylabel = kwargs.pop('ylabel', label_number(description, unit))
        ylabel = kwargs.pop('ylabel', label_number("Current density", unit))

        # Create the plot.
        ax = self.barfig(names=[name_template % (i_y+1) for i_y in range(1)],
                         times=times, xlabel=xlabel, ylabel=ylabel,
                         leg_kwargs=None, **kwargs)
        for a, time in zip(ax, times):
            a.axhline(self.get_values('cell.iprimeprime', time),
                      linestyle='--', color='k', label='Entire cell')

        # Decorate.
        if title is None:
            if self.n_z == 1:
                plt.title("Current Distribution of Cell Segments")
            else:
                plt.title("Current Distribution of Cell Segments\n"
                          "z-axis index %i (of %i)" % (z_index+1, self.n_z))
        if leg_kwargs is not None:
            loc = leg_kwargs.pop('loc', 'best')
            if len(ax) == 1:
                ax[0].legend(loc=loc, **leg_kwargs)
            else:
                plt.figlegend(ax[0].lines, **leg_kwargs)

    def gen_subtitles_time(self, times):
        """Generate titles for subplots at a list of times.
        """
        time_unit = unit2tex(self.get_unit('Time'))
        subtitles = ["t = %s" % (time if time == 0 else
            label_quantity(time, time_unit)) for time in times]
        start_time, stop_time = self.get_times('Time', [0, -1])
        for i, time in enumerate(times):
            if time == start_time:
                subtitles[i] += " (initial)"
            elif time == stop_time:
                subtitles[i] += " (final)"
        return subtitles

    def get_vector_comps2D(self, axis, index, vect, times):
        """Retrieve 2-dimensional vector components at specified times.

        **Arguments:**

        - *axis*: Slice axis ('x', 'y', or 'z')

        - *index*: Position along *axis*

        - *vect*: Name of the vector to be extracted (relative to the
          subregions)

             *vect* should contain '%i' where the vector indices should be
             inserted; the vector will be indexed according to the choice of
             *slice_axis* (see below).  *vect* will be appended to the names of
             the subregions (separated by '.') in order to determine the
             full names of the vectors.

        - *times*: Times at which to sample the vectors

        **Returns:**

        1. Tuple of the x-axis vector component (u) and the y-axis vector
           component (v)

              Here, the x axis is in screen coordinates (left/right)---not
              model coordinates. The y axis is in also screen coordinates
              (bottom/top).

        2. Tuple of the number of subregions along the x- and y- axes

              This will not be the same as (*self.n_x*, *self.n_y*) if a
              chemical component is not modeled within one or more of the
              subregions.

        3. Unit of the vector components

        4. 2D boolean list indicating if a subregion actually contains the
           species)
        """
        # Create a 2D slice of the 3D space of subregions.
        names = self.subregions_slice2D(axis, index)

        # Numbers of subregions along the cross axes
        if axis == 'x':
            # Cross axes
            axes = ('y', 'z')
            n_x, n_y = self.n_y, self.n_z
        elif axis =='y':
            # Cross axes (wrap around through Cartesion space)
            axes = ('z', 'x')
            n_x, n_y = self.n_z, self.n_x
        else:
            axes = ('x', 'y')
            n_x, n_y = self.n_x, self.n_y

        # Generate the names of the variables for the u and v components of
        # the vector.  It is necessary to look up the model's indices for
        # those components.
        names_u = [[names[i_x][i_y] + '.' + vect % self.get_IV(names[i_x][i_y]
            + '.i_' + axes[0]) for i_y in range(n_y)] for i_x in range(n_x)]
        names_v = [[names[i_x][i_y] + '.' + vect % self.get_IV(names[i_x][i_y]
            + '.i_' + axes[1]) for i_y in range(n_y)] for i_x in range(n_x)]
        isvalid = [[(self.get_IV(names[i_x][i_y] + '.%sOpen' % axes[0]) == 1,
            self.get_IV(names[i_x][i_y] + '.%sOpen' % axes[1]) == 1)
            for i_y in range(n_y)]  for i_x in range(n_x)]

        # Extract the time sequences of the components at each position.
        u_traj = self.get_values_at_times(names_u, times)
        v_traj = self.get_values_at_times(names_v, times)

        # If both of the vector components are missing, it's a good indication
        # that the species isn't included in the subregion.  Mark it as blank.
        blanks = [[u_traj[i_x][i_y] is None and v_traj[i_x][i_y] is None
                  for i_y in range(n_y)] for i_x in range(n_x)]

        # Retrieve the vector components.
        u = []
        v = []
        zeros = np.zeros((n_y, n_x)) # y is 1st index (row)
        for i in range(len(times)):
            u.append(zeros)
            v.append(zeros)
            for i_x in range(n_x):
                for i_y in range(n_y):
                    if u_traj[i_x][i_y] is not None:
                        u[-1][i_y, i_x] = u_traj[i_x][i_y][i]
                    if v_traj[i_x][i_y] is not None:
                        v[-1][i_y, i_x] = v_traj[i_x][i_y][i]

        # Determine the unit of the data.  Assume that the unit of all the
        # variables are the same; there don't seem to be any cases where this
        # wouldn't be true.
        #unit = set(self.get_unit(names_u[0][0]))
        unit = self.get_unit(names_u[0][0])

        return (u, v), (n_x, n_y), unit, blanks

    def label_layers(self, axs, y=-0.11, ygap=0.03, shrink=0.2):
        """Label the layers along the x axis of plot(s).

        **Arguments:**

        - *axs*: List of axes which should be labeled

        - *y*: Vertical position of the horizontal grouping bar

        - *ygap*: Vertical gap between the grouping bar and the top of the
          labels

        - *shrink*: Fraction that each grouping bar should shrink away from
          touching its neighbor

             A value of 0.5 causes the bar to disappear.
        """
        # TODO: An attempt to make y automatic...
        #bboxes = []
        #for label in axs[-1].get_xticklabels():
        #    bbox = label.get_window_extent()
        #    # The figure transform goes from relative coords->pixels and we
        #    # want the inverse of that.
        #    bboxi = bbox.inverse_transformed(fig.transFigure)
        #    bboxes.append(bboxi)
        ## This is the bbox that bounds all the bboxes, again in relative
        ## figure coords.
        #bbox = mtransforms.Bbox.union(bboxes)
        #y = bbox.y0
        #

        # Generate a list of the layer names.
        n_x = self.n_x
        layers = [self.locate_layer(i)[0] for i in range(n_x)]
        layers.append('') # Add one more entry; be sure it is different.

        # Label the layers.
        for ax in axs:
            i_start = 0
            for i in range(self.n_x):
                if layers[i+1] != layers[i]:
                    ax.text(0.5*(i_start+i)/(n_x-1), y-ygap,
                            layers[i], ha='center', va='top',
                            transform=ax.transAxes)
                    ax.annotate('', xy=((i_start-0.5+shrink)/(n_x-1), y),
                                xytext=((i+0.5-shrink)/(n_x-1), y),
                                xycoords='axes fraction',
                                arrowprops=dict(arrowstyle='-', shrinkA=0,
                                shrinkB=0))
                    i_start = i + 1

    def layer(self, i):
        """Return the (abbreviated) name of the layer at a given x index
        (within the entire cell) and the x index within that layer.

        **Arguments:**

        - *i*: Index of the x axis (zero-based) within the entire cell---anode
          to cathode
        """
        n = 0
        for layer in self.LAYERS:
            if i < n + self.n_x_layers[layer]:
                return layer, i-n
            n += self.n_x_layers[layer]
        return None

    def plotfig_subregions(self, prop, **kwargs):
        """Plot a property within all subregions (by default, vs. time).

        **Arguments:**

        - *prop*: Name of the property

             This will be appended to the names of the subregions (separated by
             '.') in order to determine the names of the variables to plot on
             the primary y axis.

        - *\*\*kwargs*: Additional arguments for  :meth:`modelicares.plot` (and thus to
          :meth:`modelicares.base.plot` and finally to
          :meth:`matplotlib.pyplot.plot`)
        """
        # Generate reasonable defaults.
        title = kwargs.pop('title', "%s within Subregions" % prop)
        label = kwargs.pop('label', "subregion_prop_%s" % prop)
        ynames1 = kwargs.pop('ynames1', flatten_list(self.subregions_w_prop(prop)))
        default_legends1 = [subregion.replace('.subregions', "")
                            for subregion in flatten_list(self.rel_subregions)]
        legends1 = kwargs.pop('legends1', default_legends1)

        # Create the plot.
        self.plotfig(title=title, label=label, ynames1=ynames1,
                     legends1=legends1, **kwargs)

    def quiverfig_subregions(self, vect, times=[0], n_rows=1,
                             title="", subtitles=None, label="quiver",
                             slice_axis='z', slice_index=0,
                             xlabel=None, xticklabels=None,
                             ylabel=None, yticklabels=None,
                             margin_left=rcParams['figure.subplot.left'],
                             margin_right=1-rcParams['figure.subplot.right'],
                             margin_bottom=rcParams['figure.subplot.bottom'],
                             margin_top=1-rcParams['figure.subplot.top'],
                             wspace=0.1, hspace=0.25,
                             **kwargs):
        """Create a figure with 2D vector data at given time(s) as arrows in 2D
        Cartesian coordinates.

        **Arguments:**

        - *vect*: Name of the vector to be plotted

             *vect* should contain "%i" where the vector indices should be
             inserted; the vector will be indexed according to the choice of
             *slice_axis* (see below).  *vect* will be appended to the names of
             the subregions (separated by '.') in order to determine the full
             names of the vectors.

        - *times*: List of times at which the data should be sampled

             If multiple times are given, then subfigures will be generated.

        - *n_rows*: Number of rows of (sub)plots

        - *title*: Title for the figure

        - *subtitles*: List of subtitles (i.e., titles for each subplot)

             If not provided, "t = *xx* s" will be used, where *xx* is  the
             time of each entry.  "(initial)" or "(final)" is appended if
             appropriate.

        - *label*: Label for the figure

             This will be used as a base filename if the figure is saved.

        - *slice_axis*: Axis normal to the screen or page ('x', 'y', or 'z')

        - *slice_index*: Position along *slice_axis*

        - *xlabel*: Label for the x-axes (only shown for the subplots on the
          bottom row)

             If *xlabel* is *None*, then the axis is labeled with "X Axis",
             "Y Axis", or "Z Axis" (as appropriate).

        - *xticklabels*: Labels for the x-axis ticks

             X-axis ticks are only shown for the subplots on the bottom row.
             If *xticklabels* is *None*, then the tick labels are "1", "2", ...

        - *ylabel*: Label for the y axis (only shown for the subplots on the
          left column)

             If *ylabel* is *None*, then the axis is labeled with "X Axis",
             "Y Axis", or "Z Axis" (as appropriate).

        - *yticklabels*: Labels for the y-axis ticks

             Y-axis ticks are only shown for the subplots on the left column.
             If *None*, then the tick labels are "1", "2", ...

        - *margin_left*: Left margin

        - *margin_right*: Right margin

        - *margin_bottom*: Bottom margin

        - *margin_top*: Top margin

        - *wspace*: The amount of width reserved for blank space between
          subplots

        - *hspace*: The amount of height reserved for white space between
          subplots

        - *\*\*kwargs*: Additional arguments for  :meth:`modelicares.base.quiver` (and
          thus to :meth:`matplotlib.pyplot.quiver`)
        """
        # The quiver() scaling is static and manual; therefore, it shouldn't be
        # necessary to coordinate the scaling among subplots like it is for
        # colorfig().

        # TODO 10/31/11:
        # 1) Move the dots and the pivots to the centers of masses of the
        #    subregions.
        # 2) Use and show a grid that is in scale with the geometry of the
        #    subregions (create a separate method for this it so that it can be
        #    used for colorfig too; x and y scales would need to be different).
        # 3) Overlay the velocities at the faces with the velocities of the
        #    subregions (by calling this method again with hold on?).

        from res import quiver

        # Change the default quiver color.
        color = kwargs.pop('color', 'b')

        # Slice the subregions and retrieve the data.
        (us, vs), (n_x, n_y), unit, blanks = self.get_vector_comps2D(
            slice_axis, slice_index, vect, times)

        # Generate xlabel, xticks, and xticklabels.
        if xlabel is None:
            xlabel = ('y-axis index' if slice_axis == 'x' else 'z-axis index'
                      if slice_axis == 'y' else 'x-axis index')
        if xticklabels is None:
            xticklabels = [str(i+1) for i in range(n_x)]
        xticks = range(n_x)
        #xticks = [0,1,2,10] # TODO temp

        # Generate ylabel, yticks, and yticklabels.
        if ylabel is None:
            ylabel = ('z-axis index' if slice_axis == 'x' else 'x-axis index'
                      if slice_axis == 'y' else 'y-axis index')
        if yticklabels is None:
            yticklabels = [str(i+1) for i in range(n_y)]
        yticks = range(n_y)

        # Set up the subplots.
        n_plots = len(times) # Number of plots
        if not subtitles:
            subtitles = self.gen_subtitles_time(times)
        axs, n_cols = setup_subplots(n_plots=n_plots, n_rows=n_rows,
            title=title, subtitles=subtitles, label=label,
            xlabel=xlabel, xticklabels=xticklabels, xticks=xticks,
            ylabel=ylabel, yticklabels=yticklabels, yticks=yticks,
            margin_left=margin_left, margin_right=margin_right,
            margin_bottom=margin_bottom, margin_top=margin_top,
            wspace=wspace, hspace=hspace)

        # Create the plots.
        # scale = 1e-6 # 1 um/s; Zero may cause trouble wth quiverkey() TODO ?
        quivers = []
        scale = 0
        for ax, u, v in zip(axs, us, vs):
            quivers.append(quiver(ax, u, v, pad=2*0.5, color=color, **kwargs))
            # TODO test moving the locations
            #quivers.append(quiver(ax, x = [0,1,2,10], y = [0], u=u, v=v,
            #    **kwargs))

            # Add dots at the pivots.
            for i_x in range(n_x):
                for i_y in range(n_y):
                    ax.plot(i_x, i_y, marker='o', color='k',
                            markerfacecolor='w' if blanks[i_x][i_y] else 'k')

            # Find the maximum magnitude of the data.
            scale = max(scale, np.amax(np.abs(u)), np.amax(np.abs(v)))
        # TODO: Hack: For some reason, the last plot is badly scaled.  Update
        # it to match the one before it.
        if n_plots > 1:
            axs[-1].axis(axs[-2].axis())

        # Add the key.
        # Round the scale to one significant digit.
        pow10 = get_pow10(scale)
        scale = np.round(scale/10**pow10)*10**pow10
        axs[n_cols-1].quiverkey(Q=quivers[n_cols-1], X=1, Y=1.15, U=scale,
                                label=label_quantity(scale, unit),
                                labelpos='W')

        # Label the layers.
        self.label_layers(axs[len(axs)-n_cols:])

    def sankeyfig_energy(self):
        """TODO
        """
        pass

    def subregions_slice2D(self, axis, index):
        """Return a 2-dimensional slice of the names of the subregions.

        **Arguments:**

        - *axis*: The slice axis ('x', 'y', or 'z')

        - *index*: The index of the slice axis
        """
        # It would be possible to create a class for subregions that has a
        # __getiter__ method to allow direct slicing (e.g.,
        # http://stackoverflow.com/questions/2936863/python-implementing-slicing-in-getitem).
        # However, it is simpler for now to use the brute-force method.
        if axis == 'x':
            return [[self.subregions[index][i_y][i_z]
                    for i_z in range(self.n_z)] for i_y in range(self.n_y)]
        elif axis == 'y':
            return [[self.subregions[i_x][index][i_z]
                    for i_x in range(self.n_x)] for i_z in range(self.n_z)]
        else:
            return [[self.subregions[i_x][i_y][index]
                    for i_y in range(self.n_y)] for i_x in range(self.n_x)]

    def subregions_w_prop(self, prop):
        """Return a list of the names of the subregions appended with the name
        of a property.

        **Arguments:**

        - *prop*: Name of the property to be appended

             The '.' separator is automatically inserted (between the names of
             the subregion and property).
        """
        return presuffix(self.subregions, suffix='.' + prop)

# Define the reactions and species.
#REACTIONS = {'AqHO':'Aqeous H Oxidation', 'AqOR':'Aqueous O Reduction',
#             'HO':'H Oxidation', 'OR':'O Reduction',
#             'H2OCond':'H_2O Condensation'}
#SPECIES = {'eminus':'e^-', 'H2':'H_2', 'H2Og':'H_2O_{(g)}',
#           'H2Ol':'H_2O_{(l)}', 'H3Oplus':'H_3O^+', 'Hplus':'H^+', 'N2':'N_2',
#           'O2':'O_2'}

## Create the temporal plots.

# Plot the species concentration(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.N'
#for i in range(1):
#for i in range(len(SPECIES)):
#    N = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                        SPECIES[i].name))
#    L_x = getElementData(s, n, t, storageElement, '.face.x.L')
#    L_y = getElementData(s, n, t, storageElement, '.face.y.L')
#    L_z = getElementData(s, n, t, storageElement, '.face.z.L')
#    conc = N/(L_x*L_y*L_z)
#    plotTrends(curFig, reshape(conc,size(conc,1),[]), t, storageElement,
#               label_number("Substance Concentration", 'mol/m3'),
#               ['Concentration of %s'SPECIES[i].text],
#               [fullfile(cd,'timeConc_'), SPECIES[i].name])

# Plot the species storage/consumption rate(s).
#if createAll or ishandle(curFig):
#    dataStr = '.face.%s.Ndot'
#for i in range(len(SPECIES)):
#    data = getElementData(s, n, t, storageElement,
#                          sprintf(dataStr, SPECIES[i].name))
#    plotTrends(curFig, reshape(data,size(data,1),[]), t, storageElement,
#               label_number("Substance Rate", 'mol/s'),
#               ['Storage/consumption Rate of %s' % SPECIES[i].text],
#               [fullfile(cd,'timeStorage_'), SPECIES[i].name])

# Plot the reaction rate(s).
#if createAll or ishandle(curFig):
#for i in range(len(reaction)):
#    dataStr = '.%s.Ndot_react'
#    data = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                           reaction[i].name))
#    plotTrends(curFig, reshape(data,size(data,1),[]), t, storageElement,
#               label_number("Substance Rate", 'mol/s'),
#               [reaction[i].text, ' Rate'], [fullfile(cd,'timeRate_'),
#               reaction[i].name])


# Plot the species transfer rate(s).
#if createAll or ishandle(curFig):
#dataStr1 = '.face_n.%s.Ndot'
#dataStr2 = '.face_p.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    Ndot_n = getElementData(s, n, t, transferElement_x, sprintf(dataStr1,
#                            SPECIES[i].name))
#    Ndot_p = getElementData(s, n, t, transferElement_x, sprintf(dataStr2,
#                            SPECIES[i].name))
#    Ndot = Ndot_n - Ndot_p
#    plotTrends(curFig, reshape(Ndot,size(Ndot,1),[]), t, transferElement_x,
#               label_number("Substance Rate", 'mol/s'),
#               ['Transfer Rate of %s'%SPECIES[i].text, ' in the x direction'],
#               [fullfile(cd,'timeTransferx_'), SPECIES[i].name])

## Create the spatial plots.

# Subsample the time vector.
#t = [0:5]*t(end)/5
# The dimensions of the time vector determine the number of horizontal and
# vertical subplots.
#t = [t(1:3); t(4:6)];


# Plot the species concentration(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.N'
#for i in range(1):
#for i in range(len(SPECIES)):
#    N = getElementData(s, n, t, storageElement, sprintf(dataStr,
#                                                        SPECIES[i].name))
#    L_x = getElementData(s, n, t, storageElement, '.face.x.L')
#    L_y = getElementData(s, n, t, storageElement, '.face.y.L')
#    L_z = getElementData(s, n, t, storageElement, '.face.z.L')
#    conc = N./(L_x.*L_y.*L_z)
#    plotProfiles(curFig, conc, t, n_xStorage, 'Substance Concentration /mol '
#        'm^{-3}',['Concentration of ', SPECIES[i].text],
#        [fullfile(cd,'spaceConc_'),SPECIES[i].name])

# Plot the species storage/consumption rate(s).
#if createAll or ishandle(curFig):
#dataStr = '.face.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    data = getElementData(s, n, t, storageElement,
#        sprintf(dataStr,SPECIES[i].name))
#    plotProfiles(curFig, data, t, n_xStorage, 'Substance Rate /mol s^{-1}',
#        ['Storage/consumption Rate of ', SPECIES[i].text],
#        [fullfile(cd,'spaceStorage_'),SPECIES[i].name])

# Plot the reaction rate(s).
#if createAll or ishandle(curFig):
#for i in 3:3:
#for i = 1:length(reaction)
#    dataStr = '.%s.Ndot_react'
#    data = getElementData(s, n, t, storageElement,
#        sprintf(dataStr,reaction[i].name))
#    plotProfiles(curFig, data, t, n_xStorage, 'Substance Rate /mol s^{-1}',
#        [reaction[i].text, ' Rate'],[fullfile(cd,'spaceRate_'),
#        reaction[i].name])

# Plot the species transfer rate(s).
#if createAll or ishandle(curFig):
#dataStr1 = '.face_n.%s.Ndot'
#dataStr2 = '.face_p.%s.Ndot'
#for i in range(1):
#for i in range(len(SPECIES)):
#    Ndot_n = getElementData(s, n, t, transferElement_x,
#        sprintf(dataStr1,SPECIES[i].name))
#    Ndot_p = getElementData(s, n, t, transferElement_x,
#        sprintf(dataStr2,SPECIES[i].name))
#    Ndot = Ndot_n - Ndot_p
#    plotProfiles(curFig, Ndot, t, n_xTransfer_x, 'Substance Rate /mol s^{-1}',
#        ['Transfer Rate of ', SPECIES[i].text, ' in the x direction'],
#        [fullfile(cd,'spaceTransferx_'),SPECIES[i].name])


    def __init__(self, fname="Cell.mat", cell='cell'):
        """On initialization, load and preprocess the data.

        **Arguments:**

        - *fname*: Name of the Dymosim results trajectory file (\*.mat)

        - *cell*: Name of the cell model in the trajectory file.

             This should be relative to the top level of the simulated model
             and expressed in Modelica_ dot notation.  If the cell **is** the
             simulated model, then this should be an empty string.
        """
        self._load(fname)
        self._set_constants()

        # Save the base filename and the directory.
        self.dir, self.fbase = os.path.split(fname)
        self.fbase = os.path.splitext(self.fbase)[0]

        # Read the number of subregions in the y and z directions.
        if cell != "": cell += '.'
        self.cell = cell
        self.n_y = int(self.get_IV(cell + 'n_y'))
        self.n_z = int(self.get_IV(cell + 'n_z'))

        # Read the number of subregions in the x direction in each layer.
        self.n_x_layers = {}

        # Total number of subregions along the x axis (through-the-cell)
        self.n_x = 0
        for layer in self.LAYERS:
            n_x = 1
            while cell + layer + '.L_x[%i]' % (n_x+1) in self._traj:
                n_x += 1
            self.n_x_layers.update({layer:n_x})
            self.n_x += self.n_x_layers[layer]

        # Names of the subregions as a 3D list
        # Relative to the cell model
        self.rel_subregions = []
        for layer in self.LAYERS:
            if self.n_x_layers[layer] > 0:
                self.rel_subregions.extend([[[layer + ".subregions[%i, %i, %i]"
                    % (i_x+1, i_y+1, i_z+1) for i_z in range(self.n_z)]
                    for i_y in range(self.n_y)]
                    for i_x in range(self.n_x_layers[layer])])
        # Relative to the root of the simulated model
        self.subregions = presuffix(self.rel_subregions, prefix=cell)


def annotate_conditions(conditions, temp=True, composition=True,
                        pressure=True, humidity=True, flow=True):
    r"""Create a description of the operating conditions (e.g., to be used as a
    subtitle).

    **Arguments:**

    - *conditions*: Dictionary of operating conditions

         The dictionary may include the keys *T*, *p_an_out*, *p_ca_out*,
         *n_O2*, *p_an_out*, *p_ca_out*, *T_an_in*, *T_ca_in*, *anInletRH*,
         *caInletRH*, *anStoich*, and *caStoich*.  Any entry is missing, it will
         not be exclued.  Each entry is a tuple of (number, unit string) that
         describes a condition.

    - *temp*: *True*, if temperature should be included in the description

    - *composition*: *True*, if composition should be included in the
      description

    - *pressure*: *True*, if pressure should be included in the description

    - *humidity*: *True*, if humidity should be included in the description

    - *flow*: *True*, if the flow conditions (stoichiometry) should be included
      in the description

    **Example:**

       >>> conditions = dict(T = (60, 'degC'),
       ...                   p = (150, 'kPa'),
       ...                   n_O2 = (21, '%'),
       ...                   anInletRH = (80, '%'),
       ...                   caInletRH = (50, '%'),
       ...                   anStoich = (1.5, '1'),
       ...                   caStoich = (2.0, '1'))
       >>> annotate_conditions(conditions)
       '60$\\,^{\\circ}\\!C$, 150$\\,kPa$; An|Ca: $H_2$|Air, 80|50$\\!$$\\,\\%$$\\,$RH, 1.5|2.0$\\,$stoich'
    """
    # Main description
    main_desc = []
    if temp:
        try:
            main_desc += [label_quantity(*conditions['T'])]
            temp = False
            if 'T_an_in' in conditions:
                print("T_an_in will be ignored.")
            if 'T_ca_in' in conditions:
                print("T_ca_in will be ignored.")
        except KeyError:
            pass
    if pressure:
        try:
            main_desc += [label_quantity(*conditions['p'])]
            pressure = False
            if 'p_an_out' in conditions:
                print("p_an_out will be ignored.")
            if 'p_ca_out' in conditions:
                print("p_ca_out will be ignored.")
        except KeyError:
            pass
    main_desc = ", ".join(main_desc)

    # Anode/cathode description
    anca_desc = []
    if flow:
        try:
            anca_desc += ['%s|%sstoich' % (label_quantity(conditions['anStoich'][0],
                                                          format='%.1f'),
                                           label_quantity(*conditions['caStoich'],
                                                          format='%.1f'))]
        except KeyError:
            pass
    if composition:
        try:
            if conditions['n_O2'][0] == 1.0:
                anca_desc += ['$H_2$|$O_2$']
            else:
                anca_desc += ['$H_2$|Air']
        except KeyError:
            pass
    if temp:
        try:
            anca_desc += ['%s|%s' % (label_quantity(conditions['T_an_in'][0],
                                                    format='%.0f'),
                                     label_quantity(*conditions['T_ca_in'],
                                                    format='%.0f'))]
        except KeyError:
            pass
    if pressure:
        try:
            anca_desc += ['%s|%s' % (label_quantity(conditions['p_an_out'][0],
                                                    format='%.1f'),
                                     label_quantity(*conditions['p_ca_out'],
                                                    format='%.1f'))]
        except KeyError:
            pass
    if humidity:
        try:
            anca_desc += ['%s|%s$\,$RH' % (label_quantity(conditions['anInletRH'][0],
                                                       format='%.0f'),
                                        label_quantity(*conditions['caInletRH'],
                                                       format='%.0f'))]
        except KeyError:
            pass
    anca_desc = ", ".join(anca_desc)
    if anca_desc:
        anca_desc = "An|Ca: " + anca_desc

    # Finish.
    if not anca_desc:
        return main_desc
    if not main_desc:
        return anca_desc
    return main_desc + "; " + anca_desc

def presuffix(items, prefix='', suffix=''):
    """Add a prefix and/or suffix to every entry in a list.

    The list may be multi-dimensional.  This is useful to join names of
    subregions with the name of a property.

    **Arguments:**

    - *items*: The list of items

    - *prefix*: The prefix to add

    - *suffix*: The suffix to add

    **Example:**

       >>> presuffix(['fix', 'historic', 'date', 'view'], 'pre')
       ['prefix', 'prehistoric', 'predate', 'preview']
    """
    if not prefix and not suffix:
        return items
    if iterable(items):
        if isinstance(items[0], basestring):
            return [prefix + item + suffix for item in items]
        else:
            return [presuffix(item, prefix, suffix) for item in items]

if __name__ == '__main__':
    """Test the contents of this file."""
    import doctest
    doctest.testmod()
    exit()
