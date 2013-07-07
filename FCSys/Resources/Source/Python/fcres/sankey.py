#!/usr/bin/env python
"""Module for creating Sankey diagrams using matplotlib
"""
__author__ = "Kevin Davies"
__credits__ = ["Yannick Copin"]
__license__ = "BSD"
__version__ = "0.1"
# Original version by Yannick Copin (ycopin@ipnl.in2p3.fr) 10/2/2010, available
# at:
#     http://matplotlib.sourceforge.net/examples/api/sankey_demo_old.html
# Modifications by Kevin Davies (kdavies4@gmail.com) 6/3/2011:
#   --Used arcs for the curves (so that the widths of the paths are uniform)
#   --Converted the function to a class and created methods to join multiple
#     simple Sankey diagrams
#   --Provided handling for cases where the total of the inputs isn't 100
#     Now, the default layout is based on the assumption that the inputs sum to
#     1.  A scaling parameter can be used in other cases.
#   --The call structure was changed to be more explicit about layout,
#     including the length of the trunk, length of the paths, gap between the
#     paths, and the margin around the diagram.
#   --Allowed the lengths of paths to be adjusted individually, with an option
#     to automatically justify them
#   --The call structure was changed to make the specification of path
#     orientation more flexible.  Flows are passed through one array, with
#     inputs being positive and outputs being negative.  An orientation
#     argument specifies the direction of the arrows.  The "main" inputs/
#     outputs are now specified via an orientation of 0, and there may be
#     several of each.
#   --Added assertions to catch common calling errors
#   --Added the physical unit as a string argument to be used in the labels, so
#     that the values of the flows can usually be applied automatically
#   --Added an argument for a minimum magnitude below which flows are not shown
#   --Added a tapered trunk in the case that the flows do not sum to 0
#   --Allowed the diagram to be rotated

import numpy as np
import warnings

from matplotlib.sankey import Sankey
from matplotlib.cbook import iterable, Bunch
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from matplotlib import verbose
from texunit import label_quantity

# Angles [deg/90]
RIGHT = 0
UP = 1
# LEFT = 2
DOWN = 3


class Sankey(Sankey):
    """Sankey diagram in matplotlib

    "Sankey diagrams are a specific type of flow diagram, in which the width of
    the arrows is shown proportionally to the flow quantity.  They are
    typically used to visualize energy or material or cost transfers between
    processes."
    --http://en.wikipedia.org/wiki/Sankey_diagram, accessed 6/1/2011
    """
    def add(self, patchlabel='', flows=np.array([1.0,-1.0]),
            orientations=[0,0], labels='', trunklength=1.0, pathlengths=0.25,
             prior=None, connect=(0,0), rotation=0, **kwargs):
        """
        call signature::

          add(patchlabel='', flows=np.array([1.0,-1.0]), orientations=[0,0],
              labels='', trunklength=1.0, pathlengths=0.25, prior=None,
              connect=(0,0), rotation=0, **kwargs)

        Add a simple Sankey diagram with flows at the same hierarchical level.

        Return value is the instance of :class:`Sankey`.

        Optional keyword arguments:

          ===============   ===================================================
          Keyword           Description
          ===============   ===================================================
          *patchlabel*      label to be placed at the center of the diagram
                            Note: *label* (not *patchlabel*) will be passed to
                            create the patch through **kwargs and can be used
                            to an entry in the legend.
          *flows*           array of flow values
                            By convention, inputs are positive and outputs are
                            negative.
          *orientations*    list of orientations of the paths
                            Valid values are 1 (from/to the top), 0 (from/to
                            the left or right), or -1 (from/to the bottom).  If
                            *orientations* == 0, inputs will break in from the
                            left and outputs will break away to the right.
          *labels*          list of specifications of the labels for the flows
                            Each value may be *None* (no labels), '' (just
                            label the quantities), or a labeling string.  If a
                            single value is provided, it will be applied to all
                            flows.  If an entry is a non-empty string, then the
                            quantity for the corresponding flow will be shown
                            below the string.  However, if the *unit* of the
                            main diagram is *None*, then quantities are never
                            shown, regardless of the value of this argument.
          *trunklength*     length between the bases of the input and output
                            groups
          *pathlengths*     list of lengths of the arrows before break-in or
                            after break-away
                            If a single value is given, then it will be applied
                            to the first (inside) paths on the top and bottom,
                            and the length of all other arrows will be
                            justified accordingly.  The *pathlengths* are not
                            applied to the horizontal inputs and outputs.
          *prior*           index of the prior diagram to which this diagram
                            should be connected
          *connect*         a (prior, this) tuple indexing the flow of the
                            prior diagram and the flow of this diagram which
                            should be connected
                            If this is the first diagram or *prior* is *None*,
                            *connect* will be ignored.
          *rotation*        angle of rotation of the diagram [deg]
                            *rotation* is ignored if this diagram is connected
                            to an existing one (using *prior* and *connect*).
                            The interpretation of the *orientations* argument
                            will be rotated accordingly (e.g., if *rotation*
                            == 90, an *orientations* entry of 1 means to/from
                            the left).
          ===============   ===================================================

        Valid kwargs are :meth:`~matplotlib.patches.PathPatch` arguments:
        %(PathPatch)s
        As examples, *fill*=False and *label*="A legend entry".  By default,
        *facecolor*='#bfd1d4' (light blue) and *lineweight*=0.5.

        The indexing parameters (*prior* and *connect*) are zero-based.

        The flows are placed along the top of the diagram from the inside out
        order of their index within the *flows* list or array.  They are placed
        in along the sides of the diagram from the top down and along the
        bottom from the outside in.

        If the the sum of the inputs and outputs is nonzero, the discrepancy
        will appear as a cubic Bezier curve along the top and bottom edges of
        the trunk.

        .. seealso::

            :meth:`finish`
        """
        # Check and preprocess the arguments.
        flows = np.array(flows)
        n = flows.shape[0] # Number of flows
        if rotation == None:
            rotation = 0
        else:
            # In the code below, angles are expressed in deg/90.
            rotation /= 90.0
        assert len(orientations) == n, ("orientations and flows must have the "
                                        "same length.\norientations has "
                                        "length %d, but flows has length %d."\
                                        %(len(orientations), n))
        if getattr(labels, '__iter__', False):
        # iterable() isn't used because it would give True if labels is a
        # string.
            assert len(labels) == n, ("If labels is a list, then labels and "
                                      "flows must have the same length.\n"
                                      "labels has length %d, but flows has "
                                      "length %d."%(len(labels), n))
        else:
            labels = [labels]*n
        assert trunklength >= 0, ("trunklength is negative.\nThis isn't "
                                  "allowed, because it would cause poor "
                                  "layout.")
        if np.absolute(np.sum(flows)) > self.tolerance:
            verbose.report("The sum of the flows is nonzero (%f).\nIs the "
                           "system not at steady state?" % np.sum(flows),
                           'helpful')
        scaled_flows = self.scale*flows
        gain = sum(max(flow, 0) for flow in scaled_flows)
        loss = sum(min(flow, 0) for flow in scaled_flows)
        if not (0.5 <= gain <= 2.0):
            verbose.report("The scaled sum of the inputs is %f.\nThis may "
                           "cause poor layout.\nConsider changing the scale "
                           "so that the scaled sum is approximately 1.0."%gain,
                           'helpful')
        if not (-2.0 <= loss <= -0.5):
            verbose.report("The scaled sum of the outputs is %f.\nThis may "
                           "cause poor layout.\nConsider changing the scale "
                           "so that the scaled sum is approximately 1.0."%gain,
                           'helpful')
        if prior is not None:
            assert prior >= 0, "The index of the prior diagram is negative."
            assert min(connect) >= 0, ("At least one of the connection "
                                       "indices is negative.")
            assert prior < len(self.diagrams), ("The index of the prior "
                                                "diagram is %d, but there are "
                                                "only %d other diagrams.\nThe "
                                                "index is zero-based."%(prior,
                                                len(self.diagrams)))
            assert connect[0] < len(self.diagrams[prior].flows), \
                   ("The connection index to the source diagram is %d, but "
                    "that diagram has only %d flows.\nThe index is zero-based."
                    % (connect[0], len(self.diagrams[prior].flows)))
            assert connect[1] < n, ("The connection index to this diagram is "
                                    "%d, but this diagram has only %d flows.\n"
                                    "The index is zero-based."%(connect[1], n))
            assert self.diagrams[prior].angles[connect[0]] is not None, \
                   ("The connection cannot be made.  Check that the magnitude "
                    "of flow %d of diagram %d is greater than or equal to the "
                    "specified tolerance."%(connect[0], prior))
            flow_error = self.diagrams[prior].flows[connect[0]] \
                         + flows[connect[1]]
            assert abs(flow_error) < self.tolerance, \
                  ("The scaled sum of the connected flows is %f, which is not "
                   "within the tolerance (%f)."%(flow_error, self.tolerance))

        # Determine if the flows are inputs.
        are_inputs = [None]*n
        for i, flow in enumerate(flows):
            if flow >= self.tolerance:
                are_inputs[i] = True
            elif flow <= -self.tolerance:
                are_inputs[i] = False
            else:
                verbose.report("The magnitude of flow %d (%f) is below the "
                               "tolerance (%f).\nIt will not be shown, and it "
                               "cannot be used in a connection." % (i, flow,
                               self.tolerance), 'helpful')

        # Determine the angles of the arrows (before rotation).
        angles = [None]*n
        for i, (orient, is_input) in enumerate(zip(orientations, are_inputs)):
            if orient == 1:
                if is_input:
                    angles[i] = DOWN
                elif is_input == False:
                    # Be specific since is_input can be None.
                    angles[i] = UP
            elif orient == 0:
                if is_input is not None:
                    angles[i] = RIGHT
            else:
                assert orient == -1, ("The value of orientations[%d] is %d, "
                                     "but it must be -1, 0, or 1."%(i, orient))
                if is_input:
                    angles[i] = UP
                elif is_input == False:
                    angles[i] = DOWN

        # Justify the lengths of the paths.
        if iterable(pathlengths):
            assert len(pathlengths) == n, ("If pathlengths is a list, then "
                                           "pathlengths and flows must have "
                                           "the same length.\npathlengths has "
                                           "length %d, but flows has length "
                                           "%d."%(len(pathlengths), n))
        else:  # Make pathlengths into a list.
            urlength = pathlengths
            ullength = pathlengths
            lrlength = pathlengths
            lllength = pathlengths
            d = dict(RIGHT=pathlengths)
            pathlengths = [d.get(angle, 0) for angle in angles]
            # Determine the lengths of the top-side arrows
            # from the middle outwards.
            for i, (angle, is_input, flow) \
                in enumerate(zip(angles, are_inputs, scaled_flows)):
                if angle == DOWN and is_input:
                    pathlengths[i] = ullength
                    ullength += flow
                elif angle == UP and not is_input:
                    pathlengths[i] = urlength
                    urlength -= flow # Flow is negative for outputs.
            # Determine the lengths of the bottom-side arrows
            # from the middle outwards.
            for i, (angle, is_input, flow) \
                in enumerate(zip(angles, are_inputs, scaled_flows)[::-1]):
                if angle == UP and is_input:
                    pathlengths[n-i-1] = lllength
                    lllength += flow
                elif angle == DOWN and not is_input:
                    pathlengths[n-i-1] = lrlength
                    lrlength -= flow
            # Determine the lengths of the left-side arrows
            # from the bottom upwards.
            has_left_input = False
            for i, (angle, is_input, spec) \
                in enumerate(zip(angles, are_inputs, zip(scaled_flows,
                                                         pathlengths))[::-1]):
                if angle == RIGHT:
                    if is_input:
                        if has_left_input:
                            pathlengths[n-i-1] = 0
                        else:
                            has_left_input = True
            # Determine the lengths of the right-side arrows
            # from the top downwards.
            has_right_output = False
            for i, (angle, is_input, spec) \
                in enumerate(zip(angles, are_inputs, zip(scaled_flows,
                                                         pathlengths))):
                if angle == RIGHT:
                    if not is_input:
                        if has_right_output:
                            pathlengths[i] = 0
                        else:
                            has_right_output = True

        # Begin the subpaths, and smooth the transition if the sum of the flows
        # is nonzero.
        urpath = [(Path.MOVETO, [(self.gap - trunklength / 2.0), # Upper right
                                 gain / 2.0]),
                  (Path.LINETO, [(self.gap - trunklength / 2.0) / 2.0,
                                 gain / 2.0]),
                  (Path.CURVE4, [(self.gap - trunklength / 2.0) / 8.0,
                                 gain / 2.0]),
                  (Path.CURVE4, [(trunklength / 2.0 - self.gap) / 8.0,
                                 -loss / 2.0]),
                  (Path.LINETO, [(trunklength / 2.0 - self.gap) / 2.0,
                                 -loss / 2.0]),
                  (Path.LINETO, [(trunklength / 2.0 - self.gap),
                                 -loss / 2.0])]
        llpath = [(Path.LINETO, [(trunklength / 2.0 - self.gap), # Lower left
                                 loss / 2.0]),
                  (Path.LINETO, [(trunklength / 2.0 - self.gap) / 2.0,
                                 loss / 2.0]),
                  (Path.CURVE4, [(trunklength / 2.0 - self.gap) / 8.0,
                                 loss / 2.0]),
                  (Path.CURVE4, [(self.gap - trunklength / 2.0) / 8.0,
                                 -gain / 2.0]),
                  (Path.LINETO, [(self.gap - trunklength / 2.0) / 2.0,
                                 -gain / 2.0]),
                  (Path.LINETO, [(self.gap - trunklength / 2.0),
                                 -gain / 2.0])]
        lrpath = [(Path.LINETO, [(trunklength / 2.0 - self.gap), # Lower right
                                 loss / 2.0])]
        ulpath = [(Path.LINETO, [self.gap - trunklength / 2.0, # Upper left
                                 gain / 2.0])]

        # Add the subpaths and assign the locations of the tips and labels.
        tips = np.zeros((n,2))
        label_locations = np.zeros((n,2))
        # Add the top-side inputs and outputs from the middle outwards.
        for i, (angle, is_input, spec) \
            in enumerate(zip(angles, are_inputs,
                             zip(scaled_flows, pathlengths))):
            if angle == DOWN and is_input:
                tips[i,:], label_locations[i,:] = self._add_input(ulpath,
                                                                  angle, *spec)
            elif angle == UP and not is_input:
                tips[i,:], label_locations[i,:] = self._add_output(urpath,
                                                                  angle, *spec)
        # Add the bottom-side inputs and outputs from the middle outwards.
        for i, (angle, is_input, spec) \
            in enumerate(zip(angles, are_inputs,
                             zip(scaled_flows, pathlengths))[::-1]):
            if angle == UP and is_input:
                (tips[n-i-1,:],
                 label_locations[n-i-1,:]) = self._add_input(llpath, angle,
                                                             *spec)
            elif angle == DOWN and not is_input:
                (tips[n-i-1,:],
                 label_locations[n-i-1,:]) = self._add_output(lrpath, angle,
                                                              *spec)
        # Add the left-side inputs from the bottom upwards.
        has_left_input = False
        for i, (angle, is_input, spec) \
            in enumerate(zip(angles, are_inputs,
                             zip(scaled_flows, pathlengths))[::-1]):
            if angle == RIGHT and is_input:
                if not has_left_input:
                    # Make sure the lower path extends
                    # at least as far as the upper one.
                    if llpath[-1][1][0] > ulpath[-1][1][0]:
                        llpath.append((Path.LINETO, [ulpath[-1][1][0],
                                                     llpath[-1][1][1]]))
                    has_left_input = True
                (tips[n-i-1,:],
                 label_locations[n-i-1,:]) = self._add_input(llpath, angle,
                                                             *spec)
        # Add the right-side outputs from the top downwards.
        has_right_output = False
        for i, (angle, is_input, spec) \
            in enumerate(zip(angles, are_inputs,
                             zip(scaled_flows, pathlengths))):
            if angle == RIGHT and not is_input:
                if not has_right_output:
                    # Make sure the upper path extends
                    # at least as far as the lower one.
                    if urpath[-1][1][0] < lrpath[-1][1][0]:
                        urpath.append((Path.LINETO, [lrpath[-1][1][0],
                                                     urpath[-1][1][1]]))
                    has_right_output = True
                (tips[i,:],
                 label_locations[i,:]) = self._add_output(urpath, angle, *spec)
        # Trim any hanging vertices.
        if not has_left_input:
            ulpath.pop()
            llpath.pop()
        if not has_right_output:
            lrpath.pop()
            urpath.pop()

        # Concatenate the subpaths in the correct order (clockwise from top).
        path = (urpath + self._revert(lrpath) + llpath + self._revert(ulpath) +
                [(Path.CLOSEPOLY, urpath[0][1])])

        # Create a patch with the Sankey outline.
        codes, vertices = zip(*path)
        vertices = np.array(vertices)
        def _get_angle(a, r):
            if a is None: return None
            else: return a + r

        if prior is None:
            if rotation != 0: # By default, none of this is needed.
                angles = [_get_angle(angle, rotation) for angle in angles]
                rotate = Affine2D().rotate_deg(rotation*90).transform_point
                tips = rotate(tips)
                label_locations = rotate(label_locations)
                vertices = rotate(vertices)
            text = self.ax.text(0, 0, s=patchlabel, ha='center', va='center')
        else:
            rotation = (self.diagrams[prior].angles[connect[0]] -
                        angles[connect[1]])
            angles = [_get_angle(angle, rotation) for angle in angles]
            rotate = Affine2D().rotate_deg(rotation*90).transform_point
            tips = rotate(tips)
            offset = self.diagrams[prior].tips[connect[0]] - tips[connect[1]]
            translate = Affine2D().translate(*offset).transform_point
            tips = translate(tips)
            label_locations = translate(rotate(label_locations))
            vertices = translate(rotate(vertices))
            kwds = dict(s=patchlabel, ha='center', va='center')
            text = self.ax.text(*offset, **kwds)
        if False: # Debug
            print("llpath\n", llpath)
            print("ulpath\n", self._revert(ulpath))
            print("urpath\n", urpath)
            print("lrpath\n", self._revert(lrpath))
            xs, ys = zip(*vertices)
            self.ax.plot(xs, ys, 'go-')
        patch = PathPatch(Path(vertices, codes),
                          fc=kwargs.pop('fc', kwargs.pop('facecolor',
                                        '#bfd1d4')), # Custom defaults
                          lw=kwargs.pop('lw', kwargs.pop('linewidth',
                                        0.5)),
                          **kwargs)
        self.ax.add_patch(patch)

        # Add the path labels.
        texts = []
        for i, (number, angle, label, location) in enumerate(zip(flows, angles,
                                                     labels, label_locations)):
            if label is None or angle is None:
                label = ''
            elif self.unit is not None:
                quantity = label_quantity(abs(number), self.unit, self.format)
                if label != '':
                    label += "\n"
                label += quantity
            texts.append(self.ax.text(x=location[0], y=location[1],
                                      s=label,
                                      ha='center', va='center'))
        # Text objects are placed even they are empty (as long as the magnitude
        # of the corresponding flow is larger than the tolerance) in case the
        # user wants to provide labels later.

        # Expand the size of the diagram if necessary.
        self.extent = (min(np.min(vertices[:,0]), np.min(label_locations[:,0]),
                           self.extent[0]),
                       max(np.max(vertices[:,0]), np.max(label_locations[:,0]),
                           self.extent[1]),
                       min(np.min(vertices[:,1]), np.min(label_locations[:,1]),
                           self.extent[2]),
                       max(np.max(vertices[:,1]), np.max(label_locations[:,1]),
                           self.extent[3]))
        # Include both vertices _and_ label locations in the extents; there are
        # where either could determine the margins (e.g., arrow shoulders).

        # Add this diagram as a subdiagram.
        self.diagrams.append(Bunch(patch=patch, flows=flows, angles=angles,
                                        tips=tips, text=text, texts=texts))

        # Allow a daisy-chained call structure (see docstring for the class).
        return self
