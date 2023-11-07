
# =============================================================================='
# sividl.devices
# ==============================================================================
# This module contains phidl wrapper classes used for quick generation of
# nanophotonic structures for e-beam lithography.
#
# Implemented devices are so far:
#   - Write field with alignment marker
#   - Etch slabs
#   - Generator of rectangular labelled parameter sweeps
#   - Representation of images as pixel array
#
# Note that all dimensions are in micrometer.
# ==============================================================================

from __future__ import absolute_import, division, print_function

import copy
import itertools as it
import string
from collections.abc import Iterable

import gdspy
from matplotlib.font_manager import FontProperties
import numpy as np
from phidl import CrossSection, Device
import phidl.geometry as pg
import phidl.path as pp
from sividl.sividl_utils import image_to_binary_bitmap, render_text


class SividdleDevice(Device):
    """Global device class holding all common class functions.

    Parameters
    ----------
    name: string
        Name of top-level cell.
    """

    def __init__(self, name):
        Device.__init__(self, name=name)

    def invert(self, layer, padding=None):
        """Inverts a pattern from positive to negative or vice versa.

        Parameters
        ----------
        layer: string
            Layer of new, inverted device.
        padding: np.array
            Array containing padding of bounding box used to substract device
            [left, right, top, bottom]
        """
        bounding_box = Device('interim_bounding_box')

        if padding is None:
            padding = [0, 0, 0, 0]

        bounding_box.add_polygon(
            [
                (self.xmin - padding[0], self.ymin - padding[3]),
                (self.xmin - padding[0], self.ymax + padding[2]),
                (self.xmax + padding[1], self.ymax + padding[2]),
                (self.xmax + padding[1], self.ymin - padding[3])
            ],
            layer=layer
        )

        # perform subtraction for positive e-beam resist
        inverse = pg.boolean(
            A=bounding_box,
            B=self,
            operation='A-B',
            layer=layer
        )

        return inverse

    def add_device_label(self, params):
        """Adding a label to device.

        Parameters
        ----------
        params: dict
            Dictionary containing the settings.
        params['orientation']: string
            Indicators 'r', 'l', 't', 'b' to indicate
            label position with reference to device
            (right, left, top, bottom).
        params['layer']: int
            Layer of text.
        params['text']: string
            Text to render.
        params['style']: string
            Font style,
            from here:
            https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/fonts_demo.html
        params['distance']: float
            Distance between label and device.
        """
        assert params['orientation'] in ['r', 'l', 't', 'b'], \
            "Orientation must be one of the following: ['r', 'l', 't', 'b']"

        # Add text
        params['name'] = 'label_{}'.format(params['text'])
        text_device = RenderedText(params)

        # Shift center of bounding box to origin,
        # text center should now overlapp with center of device.
        text_device.center = [0, 0]

        # Move text.
        if params['orientation'] == 'l':
            text_device.movex(
                -(text_device.xsize + self.xsize) / 2 - params['distance']
            )
        elif params['orientation'] == 'r':
            text_device.movex(
                (text_device.xsize + self.xsize) / 2 + params['distance']
            )
        elif params['orientation'] == 't':
            text_device.movey(
                (text_device.ysize + self.ysize) / 2 + params['distance']
            )
        elif params['orientation'] == 'b':
            text_device.movey(
                -(text_device.ysize + self.ysize) / 2 - params['distance']
            )

        self << text_device


class BoundingBox(SividdleDevice):
    """Contains a writefield and all plot functions.

    Parameters
    ----------
    layer: int
        Layer of bounding box.
    wf_size: int
        Write filed size in um.
    """

    def __init__(self, layer, wf_size):
        # Properly instantiate Device
        SividdleDevice.__init__(self, name='writefield_boundingbox')
        self.add_polygon([(-wf_size * 0.5, -wf_size * 0.5),
                          (-wf_size * 0.5, wf_size * 0.5),
                          (wf_size * 0.5, wf_size * 0.5),
                          (wf_size * 0.5, -wf_size * 0.5)], layer=layer)


class Ellipse(SividdleDevice):
    """Device wrapper for gdspy ellipse

    Parameters
    ----------
    params: dict
        Dictionary containing the following
        configuration parameter.
    params['layer']: int
        Layer of ellipse.
    params['rx']: float
        ellipse axis1 radius
    params['ry']: float
        ellipse axis2 radius
    """
    def __init__(self, params):
        SividdleDevice.__init__(self, name='ellipse')

        # Center device
        self.center = [0, 0]

        ellipse = gdspy.Round(
            (0, 0),
            [params['ry'], params['rx']],
            tolerance=1e-4,
            layer=params['layer']
        )

        self.add(ellipse)


class CrossAligmentMark(SividdleDevice):
    """Write alignment marker.

    Parameters
    ----------
    params: dict
        Dictionary containing the following
        configuration parameter.
    params['layer']: int
        Layer of aligment mark.
    params['invert']: Boolean
        If true, invert pattern.
    params['d_small']: int
        Width of small rectangle.
    params['d_large']:  int
        Width of large rectangle.
    params['sep']: int
        Gap between rectangles.
    params['exposure_box']: Boolean
        If True, add circle around alignment markers,
        used for marker freeing in aligned writes.
    params['exposure_box_dx']: float
        Distance between edge of alignment mark and
        edge of box.
    params['exposure_box_layer']:
        Layer of exposure box.
    """

    def __init__(self, params):
        SividdleDevice.__init__(self, name='aligment_mark')

        interim_alignment_mark = SividdleDevice('interim_alignment_mark')

        # Add four squares defining the alignment mark
        interim_alignment_mark << pg.rectangle(
            size=(
                params['d_large'],
                params['d_large']
            ),
            layer=params['layer']
        )

        interim_alignment_mark << pg.rectangle(
            size=(
                params['d_small'],
                params['d_small']
            ),
            layer=params['layer']
        ).movex(params['d_large'] + params['sep'])

        interim_alignment_mark << pg.rectangle(
            size=(
                params['d_small'],
                params['d_small']
            ),
            layer=params['layer']
        ).movey(params['d_large'] + params['sep'])

        interim_alignment_mark << pg.rectangle(
            size=(
                params['d_large'],
                params['d_large']
            ),
            layer=params['layer']).move(
                [
                    params['d_small'] + params['sep'],
                    params['d_small'] + params['sep']
                ]
        )

        # Add interim device to self, invert if choosen.
        if params['invert']:
            self << interim_alignment_mark.invert(params['layer'])
        else:
            self << interim_alignment_mark

        # Center device
        # Shift center of bounding box to origin
        self.center = [0, 0]
        # Make marker clear window.
        if params['exposure_box']:
            exposure_box_width = self.xsize \
                + 2 * params['exposure_box_dx']

            # exposure box is really a circle b/c corners make cracks
            exposure_box = gdspy.Round(
                (0, 0),
                [exposure_box_width * 0.5, exposure_box_width * 0.5],
                tolerance=1e-4,
                layer=params['exposure_box_layer']
            )

            self.add(exposure_box)

            # Add dot at alignment point
            if params['make_dot']:
                dot = pg.rectangle(
                    size=(params['dot_size'], params['dot_size']),
                    layer=params['dot_layer']
                )
                dot.move(
                    (
                        -params['dot_size'] * 0.5,
                        - params['dot_size'] * 0.5
                    )
                )

                self << dot

    def record_dot_position(self, textlayer):
        """Read center position and record it in layer."""
        center = (self.x, self.y)
        self.add_label(
            text='Alignment mark center = ({:.2f}, {:.3f}) '.format(
                center[0],
                center[1]
            ),
            position=center,
            layer=textlayer
        )


class WriteFieldCrossAligmentMark(SividdleDevice):
    """Writefiled with four cross-type alignment markers.

    Parameters
    ----------
    params: dict
        Contains the writefield parameters:
    params['bounding_box_size']: float
        Dimension of write field.
    params['bounding_box_layer']: int
        Layer of bounding box.
    params['alignment_layer']: int:
        Layer of alignment markers.
    params['alignment_offset_dx']: int
        Offset of alignment markers
        from edge of writefield in x-direction.
    params['alignment_offset_dy']: int
        Offset of alignment markers
        from edge of writefield in x-direction.
    params['alignment_mark_params']:
        Dictionary containing the parameters for the
        alignment marker, see docstring of CrossAligmentMark.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='writefield')

        # make bounding box
        bounding_box = BoundingBox(
            params['bounding_box_layer'], params['bounding_box_size']
        )

        self << bounding_box

        # Shift center of bounding box to origin
        self.center = [0, 0]

        # Position alignment marks on writefield
        delta_x = params['alignment_offset_dx']
        delta_y = params['alignment_offset_dy']

        for i in range(4):

            alignment_mark = CrossAligmentMark(
                params['alignment_mark_params']).move(
                (
                    - delta_x * (i < 2) + delta_x * (i >= 2) ,
                    - delta_y * (i in (0, 2)) + delta_x * (i in (1, 3))
                )
            )

            if params['add_text_label']:
                alignment_mark.record_dot_position(params['text_label_layer'])

            self << alignment_mark


class EtchSlap(SividdleDevice):
    """Generate two etching strip for isotropic etching tests.

    Parameters
    ----------
    params: dict
        Dictionary containing the following parameters:

    params['expose_layer']: int
        Layer for exposed regions.
    params['id_string']: string
        Identifier of device, will be printed
        in label layer.
    params['label_layer']: int
        Layer for labelling which details
        slit/slap dimensions.
    params['length_slab']: float
        Length of the slits and the resulting slab.
    params['width_slit']: float
        Width of the two slits.
    params['width_slab']: float:
        Separation of two slits which will result
        in width of slab.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='etchslab')

        # retrieve parameters
        self.id_string = params['id_string']
        self.expose_layer = params['expose_layer']
        self.label_layer = params['label_layer']
        self.length_slab = params['length_slab']
        self.width_slit = params['width_slit']
        self.width_slab = params['width_slab']

        slit = pg.rectangle(
            size=(self.width_slit, self.length_slab),
            layer=self.expose_layer
        ).rotate(90)

        self << pg.copy(slit).movey((self.width_slit + self.width_slab) * 0.5)
        self << pg.copy(slit).movey(-(self.width_slit + self.width_slab) * 0.5)
        self.add_label(
            text='{} \n slab_width = {:.2f} \
                \n slit_width = {:.2f} \
                \n slab_length = {:.2f}'.format(
                self.id_string,
                self.width_slab,
                self.width_slit,
                self.length_slab
            ),
            position=(self.xmin, self.ymax),
            layer=self.label_layer
        )

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class EtchSlabCluster(SividdleDevice):
    """Generates multiple pairs of etching strips for isotropic etching tests.

    This is meant for making extremely long slabs and/or slabs with tethers to
    the substrate. Note that for ease of placing photonic crystals and other
    features on top, the slab cluster will be placed such that (0,0) is located
    at the leftmost edge of the leftmost slab (and centered in the
    y-direction).

    Parameters
    ----------
    params: dict
        Dictionary containing the following parameters:

    params['expose_layer']: int
        Layer for exposed regions.
    params['id_string']: string
        Identifier of device, will be printed
        in label layer.
    params['label_layer']: int
        Layer for labelling which details
        slit/slap dimensions.
    params['num_slabs']: int
        Number of adjacent slabs
    params['lengths_slab']: array of floats
        Lengths of the slits and the resulting slabs.
    params['widths_slit']: array of floats
        Widths of the pairs of slits.
    params['widths_slab']: array of floats:
        Separations of each pair of slits which will result
        in widths of slabs.
    params['slab_dx']: array of floats:
        Distances in the x direction between pairs of slits.
        Note that this array will have one fewer element than
        the other arrays.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='etchslab')

        # retrieve parameters
        self.id_string = params['id_string']
        self.expose_layer = params['expose_layer']
        self.label_layer = params['label_layer']
        self.num_slabs = params['num_slabs']
        self.lengths_slab = params['lengths_slab']
        self.widths_slit = params['widths_slit']
        self.widths_slab = params['widths_slab']
        self.slab_dx = params['slab_dx']

        for i in range(self.num_slabs):
            length_slab = self.lengths_slab[i]
            width_slit = self.widths_slit[i]
            width_slab = self.widths_slab[i]
            if i > 0:
                slab_i_dx = self.slab_dx[i - 1] + np.sum(self.lengths_slab[:i])
            else:
                slab_i_dx = 0
            slit = pg.rectangle(
                size=(width_slit, length_slab),
                layer=self.expose_layer
            ).rotate(90)

            self << pg.copy(slit).move((
                slab_i_dx,
                (width_slit + width_slab) * 0.5
            ))
            self << pg.copy(slit).move((
                slab_i_dx,
                -(width_slit + width_slab) * 0.5
            ))
            self.add_label(
                text='{} \n slab_width = {:.2f} \
                    \n slit_width = {:.2f} \
                    \n slab_length = {:.2f}'.format(
                    self.id_string,
                    width_slab,
                    width_slit,
                    length_slab
                ),
                position=(self.xmin, self.ymax),
                layer=self.label_layer
            )

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class Taper(SividdleDevice):
    """Device describing a tapering section of a waveguide.

    This device will have two ports associated left and right ends
    of the tapered sections, named 'tpport1' and 'tpport2'.

    Parameters
    ----------
    layer: int
        Layer of waveguide.
    length: float
        length of taper.
    dy_min: float
        Minimum height of taper.
    dy_max: float
        MAximum height of taper.
    """

    def __init__(self, layer, name, length, dy_min, dy_max):

        SividdleDevice.__init__(self, name=name)

        self.add_polygon(
            [
                (0, 0),
                (0, dy_min),
                (length, (dy_max - dy_min) / 2 + dy_min),
                (length, -(dy_max - dy_min) / 2)
            ],
            layer=layer
        )

        self.add_port(
            name='tpport1',
            midpoint=[0, dy_min / 2],
            width=dy_min,
            orientation=180
        )

        self.add_port(
            name='tpport2',
            midpoint=[length, dy_min / 2],
            width=dy_max,
            orientation=0
        )


class TaperedCouplerWSupport(SividdleDevice):
    """Device describing a tapering section of a waveguide with a support.

    The tapering section of a waveguide is combined with a
    house-pentagon to support the tip. This device will have
    two ports associated with its left and right ends named
    'tpport1' and 'tpport2'.

    Parameters
    ----------
    layer: int
        Layer of waveguide.
    length: float
        length of taper.
    dy_min: float
        Minimum width of taper.
    dy_max: float
        Maximum width of taper.
    beam_length: float
        length of the beam that connects the house
        to the taper
    beam_width: float
        width of the support beam that connects
        the house to the taper
    roof_height: float
    house_width: float
        w
    house_length: float

    """

    def __init__(self, layer, name, length, dy_min, dy_max, beam_length,
                 beam_width, roof_height, house_width, house_length):

        SividdleDevice.__init__(self, name=name)

        self.add_polygon(
            [
                (0, -dy_min / 2),
                (0, dy_min / 2),
                (length, dy_max / 2),
                (length, -dy_max / 2)
            ],
            layer=layer
        )

        taper_slope = (dy_min - dy_max) / (2 * length)
        beam_overlap_x = (beam_width / 2 - dy_max / 2) / taper_slope

        self.add_polygon(
            [
                (beam_overlap_x, -beam_width / 2),
                (beam_overlap_x, beam_width / 2),
                (-beam_length, beam_width / 2),
                (-beam_length, -beam_width / 2)
            ],
            layer=layer
        )

        # add triangular transition to rectangular support
        self.add_polygon(
            [
                (-beam_length, -beam_width / 2),
                (-beam_length, beam_width / 2),
                (-beam_length - roof_height, house_width / 2),
                (-beam_length - roof_height, -house_width / 2)
            ],
            layer=layer
        )

        # add rectangular support
        self.add_polygon(
            [
                (-beam_length - roof_height, -house_width / 2),
                (-beam_length - roof_height, house_width / 2),
                (-beam_length - roof_height - house_length, house_width / 2),
                (-beam_length - roof_height - house_length, -house_width / 2)
            ],
            layer=layer
        )

        self.add_port(
            name='tpport1',
            midpoint=[-beam_length - roof_height - house_length, 0],
            width=house_width,
            orientation=180
        )

        self.add_port(
            name='tpport2',
            midpoint=[length, 0],
            width=dy_max,
            orientation=0
        )


class TaperedSupport(SividdleDevice):
    """Device describing a tapered support section of a waveguide.

    This tapered support is in the style of Mike and Bart. There
    are several segments in a tapered support. In addition to the straight
    section at the support, the tapering sections can be broken
    down into 2 parts. That is, there is a concave section
    followed by a convex section as you go from waveguide width
    to support width.

    This device has two ports associated left and right ends

    Parameters
    ----------
    params: dict
        Dictionary containing the following parameters:
    params['layer']: int
        Layer of waveguide
    params['taper_length_1']: float
        tapered length on the LHS of the tapered support
    params['straight_length_center']: float
        straight length at center of the tapered support
    params['taper_length_2']: float
        tapered length on the RHS of the tapered support
    params['width_1']: float
        waveguide width at LHS of the tapered support
    params['width_center']: float
        waveguide width at center of the tapered support
    params['width_2']: float
        waveguide width at RHS of the tapered support
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        self.taper_length_1 = params['taper_length_1']
        self.taper_length_2 = params['taper_length_2']
        self.width_1 = params['width_1']
        self.width_center = params['width_center']
        self.width_2 = params['width_2']

        straight_center_path = pp.straight(
            length=params['straight_length_center'])

        taper_1_conc_path = pp.straight(length=params['taper_length_1'] / 2.0)
        taper_1_conv_path = pp.straight(length=params['taper_length_1'] / 2.0)
        taper_2_conc_path = pp.straight(length=params['taper_length_2'] / 2.0)
        taper_2_conv_path = pp.straight(length=params['taper_length_2'] / 2.0)

        # Create blank CrossSection objects to be used for each path
        w_taper1_conc = CrossSection()
        w_taper1_conv = CrossSection()
        w_center = CrossSection()
        w_taper2_conc = CrossSection()
        w_taper2_conv = CrossSection()

        # Add a single "section" to each of the cross-sections
        w_taper1_conc.add(
            width=self.tapered_width_1_conc, offset=0,
            layer=params['layer'],
            ports=('in_taper1_conc', 'out_taper1_conc')
        )
        w_taper1_conv.add(
            width=self.tapered_width_1_conv, offset=0,
            layer=params['layer'],
            ports=('in_taper1_conv', 'out_taper1_conv')
        )
        w_center.add(
            width=params['width_center'], offset=0,
            layer=params['layer'],
            ports=('in_center', 'out_center')
        )
        w_taper2_conv.add(
            width=self.tapered_width_2_conv, offset=0,
            layer=params['layer'],
            ports=('in_taper2_conv', 'out_taper2_conv')
        )
        w_taper2_conc.add(
            width=self.tapered_width_2_conc,
            offset=0, layer=params['layer'],
            ports=('in_taper2_conc', 'out_taper2_conc')
        )

        # Combine the Path and the CrossSection
        taper_1_conc_dev = taper_1_conc_path.extrude(
            w_taper1_conc
        )
        taper_1_conv_dev = taper_1_conv_path.extrude(
            w_taper1_conv
        )
        straight_center_dev = straight_center_path.extrude(
            w_center
        )
        taper_2_conv_dev = taper_2_conv_path.extrude(
            w_taper2_conv
        )
        taper_2_conc_dev = taper_2_conc_path.extrude(
            w_taper2_conc
        )

        taper_1_conc_ref = self.add_ref(taper_1_conc_dev)
        taper_1_conv_ref = self.add_ref(taper_1_conv_dev)
        straight_center_ref = self.add_ref(straight_center_dev)
        taper_2_conc_ref = self.add_ref(taper_2_conc_dev)
        taper_2_conv_ref = self.add_ref(taper_2_conv_dev)

        taper_1_conv_ref.connect('in_taper1_conv',
                                 taper_1_conc_ref.ports['out_taper1_conc'])
        straight_center_ref.connect('in_center',
                                    taper_1_conv_ref.ports['out_taper1_conv'])
        taper_2_conv_ref.connect('in_taper2_conv',
                                 straight_center_ref.ports['out_center'])
        taper_2_conc_ref.connect('in_taper2_conc',
                                 taper_2_conv_ref.ports['out_taper2_conv'])

        self.add_port(name='tpport1',
                      port=taper_1_conc_ref.ports['in_taper1_conc'])
        self.add_port(name='tpport2',
                      port=taper_2_conc_ref.ports['out_taper2_conc'])

    def tapered_width_1_conc(self, x):
        tapered_support_width1 = self.width_1 + 0.5 * (
            self.width_center - self.width_1) * (x)**2
        return tapered_support_width1

    def tapered_width_1_conv(self, x):
        tapered_support_width1 = self.width_center - 0.5 * (
            self.width_center - self.width_1) * ((x - 1))**2
        return tapered_support_width1

    def tapered_width_2_conv(self, x):
        tapered_support_width2 = self.width_center - 0.5 * (
            self.width_center - self.width_2) * ((x))**2
        return tapered_support_width2

    def tapered_width_2_conc(self, x):
        tapered_support_width2 = self.width_2 + 0.5 * (
            self.width_center - self.width_2) * (x - 1)**2
        return tapered_support_width2


class WaveGuide(SividdleDevice):
    """Device describing a rectangular waveguide.

    This device will have two ports associated with the axis defined
    by the 'height' dimension, which are named ''wgport1' and 'wgport2'.
    Parameters
    ----------
    layer: int
        Layer of waveguide.
    length: float
        length of waveguide.
    height: float
        Height of waveguide.
    photonic_cristal_params: TBD
    photonic_crystal_params
    """

    def __init__(self, layer, length, height, photonic_crystal_params=None):

        SividdleDevice.__init__(self, name='waveguide')

        self.add_polygon(
            [(0, 0), (length, 0), (length, height), (0, height)],
            layer=layer
        )
        self.add_port(
            name='wgport1',
            midpoint=[0, height / 2],
            width=height,
            orientation=180
        )

        self.add_port(
            name='wgport2',
            midpoint=[length, height / 2],
            width=height,
            orientation=0
        )

        # Add photonic crytstal
        if photonic_crystal_params is not None:

            holes = AdiabaticTaperedEllipseArray(
                photonic_crystal_params['holes_layer'],
                photonic_crystal_params['hx_init'],
                photonic_crystal_params['hy_init'],
                photonic_crystal_params['hy_final'],
                photonic_crystal_params['a_const'],
                photonic_crystal_params['num_taper'],
                photonic_crystal_params['num_cells'],
                both=photonic_crystal_params['both'],
                flip=True
            )

            # wg
            self << holes.move(
                (
                    holes.xsize * 0.5 + photonic_crystal_params['dx_holes'],
                    self.ysize * 0.5
                )
            )

        # Store Layer
        self.layer = layer

    def add_anchors(self, dx_anchor, width_anchor,
                    widthmax_anchor, length_anchor, which):
        """Add anchors to the waveguide.

        Parameters
        ----------
        dx_anchor: float
            Distance of anchoring point to end of waveguide.
        width_anchor: float
            Width of anchor at anchoring point.
        widthmax_anchor: float
            Width of anchor opposite of anchoring point.
        length_anchor: float
            Length of anchor.
        which: list
            Which anchor to add, labelled:
            1 ..... 2
            3 ..... 4
            Default: All anchors.
        """
        if which is None:
            which = [1, 2, 3, 4]

        # add ports
        self.add_port(
            name='anchorport1',
            midpoint=[dx_anchor, self.ysize],
            width=width_anchor,
            orientation=90
        )

        self.add_port(
            name='anchorport2',
            midpoint=[self.xsize - dx_anchor, self.ysize],
            width=width_anchor,
            orientation=90
        )

        self.add_port(
            name='anchorport3',
            midpoint=[dx_anchor, 0],
            width=width_anchor,
            orientation=-90
        )

        self.add_port(
            name='anchorport4',
            midpoint=[self.xsize - dx_anchor, 0],
            width=width_anchor,
            orientation=-90
        )

        # add achors
        anchor = Taper(
            self.layer,
            'anchor',
            length_anchor,
            width_anchor,
            widthmax_anchor
        )

        # Addd to new device
        if 1 in which:
            anchor_1 = self << anchor
            anchor_1.connect(
                port='tpport1',
                destination=self.ports['anchorport1']
            )
            self.add_port(
                port=anchor_1.ports['tpport2'],
                name='extanchorport1'
            )
        if 2 in which:
            anchor_2 = self << anchor
            anchor_2.connect(
                port='tpport1',
                destination=self.ports['anchorport2']
            )
            self.add_port(
                port=anchor_2.ports['tpport2'],
                name='extanchorport2'
            )
        if 3 in which:
            anchor_3 = self << anchor
            anchor_3.connect(
                port='tpport1',
                destination=self.ports['anchorport3']
            )
            self.add_port(
                port=anchor_3.ports['tpport2'],
                name='extanchorport3'
            )
        if 4 in which:
            anchor_4 = self << anchor
            anchor_4.connect(
                port='tpport1',
                destination=self.ports['anchorport4']
            )
            self.add_port(
                port=anchor_4.ports['tpport2'],
                name='extanchorport4'
            )


class SlabStyleWaveguide(SividdleDevice):
    """Device describing a slab-style rectangular waveguide.

    This device will hace two ports associated with the axis defined
    by the 'height' dimension, which are named ''wgport1' and 'wgport2'.

    Parameters
    ----------
    layer: int
        Layer of waveguide.
    length: float
        length of waveguide.
    height: float
        Height of waveguide.
    window_width: float
        Defines width of windows around waveguide
    photonic_crystal: boolean
        True if there are holes in the waveguide
    photonic_crystal_params
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='ss_waveguide')

        # retrieve parameters
        layer = params['layer']
        length = params['length']
        height = params['height']
        phc = params['photonic_crystal']
        slab_style = params['slab_style']

        window_width = params['window_width']
        self.add_polygon(
            [(0, 0), (length, 0), (length, window_width),
                (0, window_width)], layer=layer
        )

        self.add_polygon(
            [(0, -height), (length, -height),
                (length, -height - window_width),
                (0, -height - window_width)], layer=layer
        )

        if phc:

            holes = AdiabaticTaperedEllipseArray(
                params['holes_layer'],
                params['hx_init'],
                params['hy_init'],
                params['hy_final'],
                params['a_const'],
                int(params['num_taper']),
                int(params['num_cells']),
                both=params['both'],
                flip=True
            )

            # wg
            dy_holes = self.ysize * 0.5
            if slab_style:
                dy_holes = -height / 2
            self << holes.move(
                (
                    self.xsize * 0.5 + params['dx_holes'],
                    dy_holes  # self.ysize * 0.5
                )
            )

        # Store Layer
        self.layer = layer

    def add_anchors(self, dx_anchor, width_anchor,
                    widthmax_anchor, length_anchor, which):
        """Add anchors to the waveguide.

        Parameters
        ----------
        dx_anchor: float
            Distance of anchoring point to end of waveguide.
        width_anchor: float
            Width of anchor at anchoring point.
        widthmax_anchor: float
            Width of anchor opposite of anchoring point.
        length_anchor: float
            Length of anchor.
        which: list
            Which anchor to add, labelled:
            1 ..... 2
            3 ..... 4
            Default: All anchors.
        """
        if which is None:
            which = [1, 2, 3, 4]

        # add ports
        self.add_port(
            name='anchorport1',
            midpoint=[dx_anchor, self.ysize],
            width=width_anchor,
            orientation=90
        )

        self.add_port(
            name='anchorport2',
            midpoint=[self.xsize - dx_anchor, self.ysize],
            width=width_anchor,
            orientation=90
        )

        self.add_port(
            name='anchorport3',
            midpoint=[dx_anchor, 0],
            width=width_anchor,
            orientation=-90
        )

        self.add_port(
            name='anchorport4',
            midpoint=[self.xsize - dx_anchor, 0],
            width=width_anchor,
            orientation=-90
        )

        # add anchors
        anchor = Taper(
            self.layer,
            'anchor',
            length_anchor,
            width_anchor,
            widthmax_anchor
        )

        # Add to new device
        if 1 in which:
            anchor_1 = self << anchor
            anchor_1.connect(
                port='tpport1',
                destination=self.ports['anchorport1']
            )
            self.add_port(
                port=anchor_1.ports['tpport2'],
                name='extanchorport1'
            )
        if 2 in which:
            anchor_2 = self << anchor
            anchor_2.connect(
                port='tpport1',
                destination=self.ports['anchorport2']
            )
            self.add_port(
                port=anchor_2.ports['tpport2'],
                name='extanchorport2'
            )
        if 3 in which:
            anchor_3 = self << anchor
            anchor_3.connect(
                port='tpport1',
                destination=self.ports['anchorport3']
            )
            self.add_port(
                port=anchor_3.ports['tpport2'],
                name='extanchorport3'
            )
        if 4 in which:
            anchor_4 = self << anchor
            anchor_4.connect(
                port='tpport1',
                destination=self.ports['anchorport4']
            )
            self.add_port(
                port=anchor_4.ports['tpport2'],
                name='extanchorport4'
            )


class TaperedWaveGuide(SividdleDevice):
    """Device describing tapered waveguide.

    Parameters
    ----------
    params: dict
        Dictionary containing the parameters of the tapered waveguide
    params['len_wg']: float
        Length of waveguide section.
    params['height_wg']: float
        Height of waveguide section.
    params['len_tp_left']: float
        Lenght of left tapered section.
    params['len_tp_right']: float
        Lenght of right tapered section.
    params['width_tp']: float
        Final width of tapered section.
    params['which_anchors']: list
        Which anchor to add, labelled:
        1 ..... 2
        3 ..... 4
    params['dx_anchor']: float
        Distance of anchoring point to end of waveguide.
    params['width_anchor']: float
        Width of anchor at anchoring point.
    params['widthmax_anchor']: float
        Width of anchor opposite of anchoring point.
    params['length_anchor']: float
        Length of anchor.
    params['invert']: boolean
        If True, invert layer by bounding box.
        This is useful if we're working with positive
        photoresists.
    params['invert_layer']: float
        Target layer for inverted design.
    params['taper_gap_left']: float
        Gap between left taper and etch box.
    params['taper_gap_right']: float
        Gap between right taper and etch box.
    params['photonic_crystal_params']: dict
        Dictionary containing the parameters for the photonic crystal
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='Tapered_Waveguide')

        self.params = params

        # Generate waveguide.
        waveguide = WaveGuide(
            params['layer'],
            params['len_wg'],
            params['height_wg'],
            params['photonic_crystal_params']
        )

        # Add anchors.
        waveguide.add_anchors(
            params["dx_anchor"],
            params["width_anchor"],
            params["widthmax_anchor"],
            params["length_anchor"],
            params['which_anchors']
        )

        # Generate tapered sections.
        taper_left, taper_right = [
            Taper(
                layer=params['layer'],
                name='taper_left',
                dy_min=params['width_tp'],
                dy_max=params['height_wg'],
                length=length
            ) for length in [params['len_tp_left'], params['len_tp_right']]
        ]

        # Add waveguide and taper to device.
        waveguide_on_device = self << waveguide
        taper_left_on_device = self << taper_left
        taper_right_on_device = self << taper_right

        # connect the ports
        waveguide_on_device.connect(
            port='wgport1',
            destination=taper_left_on_device.ports['tpport2']
        )

        taper_right_on_device.connect(
            port='tpport2',
            destination=waveguide_on_device.ports['wgport2']
        )

        # add ports to new device
        self.add_port(
            port=taper_left_on_device.ports['tpport1'],
            name='tpwgport1'
        )

        self.add_port(
            port=taper_right_on_device.ports['tpport1'],
            name='tpwgport2'
        )

        # Add the inverse of the design to another layer
        # (for positive photoresist usage)
        if params['invert']:
            self << self.invert(
                params['invert_layer'],
                [
                    params['taper_gap_left'],
                    params['taper_gap_right'],
                    0,
                    0
                ]
            )

        # The below doesn't work unless you have defined this parameter,
        # which one does not generically do.
        # if params['add_taper_marker']:
        #     self.add_arrow_markers()
        #     # Compensate for the arrows in assignment of center.
        #     self.center = [0, -params['width_tp'] * 0.5]
        # else:
        #     self.center = [0, 0]

    def add_arrow_markers(self):
        """Add arrow markers.

        This will generate three arrow markers,
        two indicating where the tapered regions start and one
        indicating which side has the left taper.

        """
        params = self.params

        up_arrow_params = {
            'fontsize'   : 12,
            'name'       : 'taperlabel',
            'style'     : 'normal',
            'layer'     : 3,
            'text'      : '↑',
        }

        # Render arrow, and move to beginning of left taper.

        init_ysize = self.ysize  # Will get modified by adding of arrow.
        up_arrow_1 = RenderedText(up_arrow_params)
        up_arrow_2 = RenderedText(up_arrow_params)
        up_arrow_3 = RenderedText(up_arrow_params)

        # Arrow pointing at end of taper section, left.
        self << up_arrow_1.move(
            [
                params['len_tp_left'],
                -((init_ysize + up_arrow_1.ysize) * 0.5
                  + 0.15 * init_ysize + 5)
            ]
        )

        # Arrow pointing at end of taper section, right.
        self << up_arrow_2.move(
            [
                params['len_tp_left'] + params['len_wg'],
                -((init_ysize + up_arrow_2.ysize) * 0.5
                  + 0.15 * init_ysize + 5)
            ]
        )

        # Arrow pointing at end of taper section, right.
        self << up_arrow_3.move(
            [
                params['len_tp_left']
                + params['len_tp_right'] + params['len_wg'],
                -((init_ysize + up_arrow_3.ysize) * 0.5
                  + 0.15 * init_ysize + 5)
            ]
        )

        # Render arrow, and move to beginning of right taper.
        right_arrow_params = up_arrow_params
        right_arrow_params['text'] = '→'
        right_arrow = RenderedText(right_arrow_params)
        self << right_arrow.move(
            [
                right_arrow.xsize * 0.5,
                -((init_ysize + right_arrow.ysize) * 0.5
                  + 0.15 * init_ysize + 5)
            ]
        )


class EllipseArray(SividdleDevice):
    """Device containing an Array of ellipses.

    These arrays can be used to construct photonic crystal cavities.
    Nonclature carried over from https://arxiv.org/pdf/1907.13200.pdf

    Parameters
    ----------
    layer: int
        Target layer.
    hx_array: numpy.array
        Width of ellipses.
    hy_array: numpy.array
        Height of ellipses.
    a_consts_array: numpy.array:
        Lattice constants.
    """

    def __init__(self, layer, hx_array, hy_array, a_consts_array):

        # Check correct array lengths
        assert len(hy_array) is len(hx_array), \
            "Please provide arrays of equal length"

        assert len(a_consts_array) is (len(hx_array) - 1), \
            "Please provide a lattice constant array with one fewer \
             element than the hole arrays."

        SividdleDevice.__init__(self, name='EllipseArray')

        # Prepend 0 to lattice constant array for convenience
        a_consts_array = np.insert(a_consts_array, 0, 0)
        shift = 0

        # Generate ellipses
        for (hx, hy, a) in zip(hx_array, hy_array, a_consts_array):
            shift += a
            ellipse = gdspy.Round(
                (0 + shift, 0),
                [hx * 0.5, hy * 0.5],
                tolerance=1e-4,
                layer=layer
            )
            self.add(ellipse)


class AdiabaticTaperedEllipseArray(SividdleDevice):
    """Device containing an Array of ellipses tapering down.

    The ellipses are tapered down from hy_init to hy_final while maintaining
    the aspect ratio. The lattice constant is assumed to be constant

    These arrays can be used to construct photonic crystal cavities.
    Nonclature carried over from https://arxiv.org/pdf/1907.13200.pdf.

    Parameters
    ----------
    hx_init: float
        Initial value of hx.
    hy_init: float
        Initial value of hx.
    hy_final: float
        Final value of hy.
    a_const: float:
        Lattice constant.
    num_taper: int
        Number of tapered cells.
    num_cells: int
        Number of non-tapered cells.
    flip: boolean
        If True, flip array. Tapered section
        now is left.
    both: boolean:
        If True, add taper on both sides.
    """

    def __init__(self, layer, hx_init, hy_init,
                 hy_final, a_const, num_taper, num_cells, flip, both):

        SividdleDevice.__init__(self, name='AdiabaticTaperedEllipseArray')

        # Arrays of untapered devices
        hx_array_init = np.ones(num_cells) * hx_init
        hy_array_init = np.ones(num_cells) * hy_init

        # Construct Hx array.
        hy_array_tapered = np.linspace(
            hy_init,
            hy_final,
            num_taper
        )

        # Construct Hx array.
        hx_array_tapered = np.linspace(
            hx_init,
            hy_final,
            num_taper
        )

        # Appending tapered and nontapered arrays
        hx_array = np.append(hx_array_init, hx_array_tapered)
        hy_array = np.append(hy_array_init, hy_array_tapered)

        if both:
            hx_array = np.append(np.flip(hx_array_tapered), hx_array)
            hy_array = np.append(np.flip(hy_array_tapered), hy_array)
            num_taper += num_taper

        if flip:
            hx_array = np.flip(hx_array)
            hy_array = np.flip(hy_array)

        # Construct array of constant lattice consants.
        a_consts_array = np.ones(num_taper + num_cells - 1) * a_const

        self << EllipseArray(layer, hx_array, hy_array, a_consts_array)

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class DoubleTaperedDevice(SividdleDevice):
    """Makes a double-tapered device.

    A device comprised of waveguide tapers on both sides,
    two quadratic tapered supports on either side of the cavity,
    and a photonic crystal cavity (PCC) zone at the center. That is,
    the PCC is defined by the following parameters and is allowed
    to be asymmetric both in hole dimensions and hole spacing.

    Parameters
    ----------
    params: dict
        Dictionary containing the following parameters:

    params['cavity_length']: float
        length of the cavity region
    params['layer_wg']: int
        Layer of waveguide
    params['name']: string
        Name of GDS cell.
    params['tapered_coupler_length']: float
        length of the tapers on either side of the device
    params['tapered_coupler_minWidth']: float
        minimum width of the tapers on either side of the device
    params['tapered_support_length']: float
        length of the taperedSupports
    params['tapered_support_width']: float
        tapered support width. Typically, this is 1.4*width
    params['waveguide_spacer_length']: float
        length of the spacers that separate the tapered coupler
        from the tapered support and that space out the two
        tapered supports on either side.
    params['width']: float
        cavity and waveguide width
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        self.layer_wg = params['layer_wg']
        self.tapered_coupler_length = params['tapered_coupler_length']
        self.waveguide_spacer_length = params['waveguide_spacer_length']
        self.cavity_length = params['cavity_length']
        self.tapered_support_length = params['tapered_support_length']
        self.tapered_coupler_minWidth = params['tapered_coupler_minWidth']
        self.width = params['width']
        self.tapered_support_width = params['tapered_support_width']

        self.tapered_support_params = {
            'name'           : "tapered_support",
            'layer'          : self.layer_wg,
            'taper_length_1' : self.tapered_support_length / 2,
            'straight_length_center': 0,
            'taper_length_2' : self.tapered_support_length / 2,
            'width_1'        : self.width,
            'width_center'   : self.tapered_support_width,
            'width_2'        : self.width
        }

        # Create paths to define segments of device
        spacer_path = pp.straight(length=self.waveguide_spacer_length)
        cavity_path = pp.straight(length=self.cavity_length)

        # Create blank CrossSection objects to be used for each path
        spacer_xs = CrossSection()
        cavity_xs = CrossSection()

        # Add a section to each CrossSection to define the width for extrusion
        spacer_xs.add(width=self.width, offset=0, layer=self.layer_wg,
                      ports=('in_spacer', 'out_spacer'))
        cavity_xs.add(width=self.width, offset=0, layer=self.layer_wg,
                      ports=('in_cavity', 'out_cavity'))

        # Combine the Paths and the CrossSections to make Devices
        spacer_dev = spacer_path.extrude(spacer_xs)
        cavity_dev = cavity_path.extrude(cavity_xs)

        # Create TaperedSupport Device
        tapered_support = TaperedSupport(self.tapered_support_params)

        # Create Taper Device for tapered couplers
        tapered_coupler = Taper(layer=self.layer_wg, name="TaperedCoupler",
                                length=self.tapered_coupler_length,
                                dy_min=self.tapered_coupler_minWidth,
                                dy_max=self.width)

        # Create Device_References of the Devices
        spacer1_ref = self.add_ref(spacer_dev)
        spacer2_ref = self.add_ref(spacer_dev)
        spacer3_ref = self.add_ref(spacer_dev)
        spacer4_ref = self.add_ref(spacer_dev)
        cavity_ref = self.add_ref(cavity_dev)
        tapered_support1_ref = self.add_ref(tapered_support)
        tapered_support2_ref = self.add_ref(tapered_support)
        tapered_support3_ref = self.add_ref(tapered_support)
        tapered_support4_ref = self.add_ref(tapered_support)
        tapered_coupler1_ref = self.add_ref(tapered_coupler)
        tapered_coupler2_ref = self.add_ref(tapered_coupler)

        # Connect up the references to build the geometry
        spacer1_ref.connect('in_spacer',
                            tapered_coupler1_ref.ports['tpport2'])
        tapered_support1_ref.connect('tpport1',
                                     spacer1_ref.ports['out_spacer'])
        spacer2_ref.connect('in_spacer', tapered_support1_ref.ports['tpport2'])
        tapered_support2_ref.connect('tpport1',
                                     spacer2_ref.ports['out_spacer'])
        cavity_ref.connect('in_cavity', tapered_support2_ref.ports['tpport2'])
        tapered_support3_ref.connect('tpport1', cavity_ref.ports['out_cavity'])
        spacer3_ref.connect('in_spacer', tapered_support3_ref.ports['tpport2'])
        tapered_support4_ref.connect('tpport1',
                                     spacer3_ref.ports['out_spacer'])
        spacer4_ref.connect('in_spacer', tapered_support4_ref.ports['tpport2'])
        tapered_coupler2_ref.connect('tpport2',
                                     spacer4_ref.ports['out_spacer'])
        # Shift center of bounding box to origin.
        self.center = [0, 0]


class DoubleTaperedDeviceWithCouplerSupports(SividdleDevice):
    """A double tapered device with coupler supports.

    A device comprised of waveguide tapers on both sides,
    two quadratic tapered supports on either side of the cavity,
    and a photonic crystal cavity (PCC) zone at the center. The PCC
    holes are defined in a separate class. The tapered waveguide
    couplers are supported at their tips with larger, non-freestanding
    structures.

    Parameters
    ----------
    params: dict
        Dictionary containing the following parameters:

    params['cavity_length']: float
        length of the cavity region
    params['layer_wg']: int
        Layer of waveguide
    params['name']: string
        Name of GDS cell.
    params['tapered_coupler_length']: float
        length of the tapers on either side of the device
    params['tapered_coupler_minWidth']: float
        minimum width of the tapers on either side of the device
    params['tapered_coupler_support_length']: float
        length of tapered coupler support
    params['tapered_coupler_support_width']: float
        max width of tapered coupler support
    params['tapered_support_length']: float
        length of the taperedSupports
    params['tapered_support_width']: float or array of floats
        tapered support width. Typically, this is 1.4*width.
        If array, then first element is width for left and second
        is for width on right.
    params['waveguide_spacer_length']: float
        length of the spacers that separate the tapered coupler
        from the tapered support and that space out the two
        tapered supports on either side.
    params['width']: float
        cavity and waveguide width
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        self.layer_wg = params['layer_wg']
        self.tapered_coupler_length = params['tapered_coupler_length']
        self.waveguide_spacer_length = params['waveguide_spacer_length']
        self.cavity_length = params['cavity_length']
        self.tapered_support_length = params['tapered_support_length']
        self.tapered_coupler_minWidth = params['tapered_coupler_minWidth']
        self.tapered_coupler_support_beam_length = params[
            'tapered_coupler_support_beam_length']
        self.tapered_coupler_support_beam_width = params[
            'tapered_coupler_support_beam_width']
        self.tapered_coupler_support_roof_height = params[
            'tapered_coupler_support_roof_height']
        self.tapered_coupler_support_house_width = params[
            'tapered_coupler_support_house_width']
        self.tapered_coupler_support_house_length = params[
            'tapered_coupler_support_house_length']
        self.width = params['width']
        self.tapered_support_width = params['tapered_support_width']
        if isinstance(self.tapered_support_width, Iterable):
            self.tapered_support_width_L = self.tapered_support_width[0]
            self.tapered_support_width_R = self.tapered_support_width[1]
        else:
            self.tapered_support_width_L = self.tapered_support_width
            self.tapered_support_width_R = self.tapered_support_width

        self.tapered_support_L_params = {
            'name'           : "tapered_support",
            'layer'          : self.layer_wg,
            'taper_length_1' : self.tapered_support_length / 2,
            'straight_length_center': 0,
            'taper_length_2' : self.tapered_support_length / 2,
            'width_1'        : self.width,
            'width_center'   : self.tapered_support_width_L,
            'width_2'        : self.width
        }

        self.tapered_support_R_params = {
            'name'           : "tapered_support",
            'layer'          : self.layer_wg,
            'taper_length_1' : self.tapered_support_length / 2,
            'straight_length_center': 0,
            'taper_length_2' : self.tapered_support_length / 2,
            'width_1'        : self.width,
            'width_center'   : self.tapered_support_width_R,
            'width_2'        : self.width
        }

        # Create paths to define segments of device
        spacer_path = pp.straight(length=self.waveguide_spacer_length)
        cavity_path = pp.straight(length=self.cavity_length)

        # Create blank CrossSection objects to be used for each path
        spacer_xs = CrossSection()
        cavity_xs = CrossSection()

        # Add a section to each CrossSection to define the width for extrusion
        spacer_xs.add(width=self.width, offset=0, layer=self.layer_wg,
                      ports=('in_spacer', 'out_spacer'))
        cavity_xs.add(width=self.width, offset=0, layer=self.layer_wg,
                      ports=('in_cavity', 'out_cavity'))

        # Combine the Paths and the CrossSections to make Devices
        spacer_dev = spacer_path.extrude(spacer_xs)
        cavity_dev = cavity_path.extrude(cavity_xs)

        # Create TaperedSupport Device
        tapered_support_L = TaperedSupport(self.tapered_support_L_params)
        tapered_support_R = TaperedSupport(self.tapered_support_R_params)

        # Create Taper Device for tapered couplers

        tapered_coupler = TaperedCouplerWSupport(
            layer=self.layer_wg,
            name="TaperedCoupler",
            length=self.tapered_coupler_length,
            dy_min=self.tapered_coupler_minWidth, dy_max=self.width,
            beam_length=self.tapered_coupler_support_beam_length,
            beam_width=self.tapered_coupler_support_beam_width,
            roof_height=self.tapered_coupler_support_roof_height,
            house_width=self.tapered_coupler_support_house_width,
            house_length=self.tapered_coupler_support_house_length
        )
        # # and for coupler supports
        # couplerSupport = Taper(layer = self.layer_wg,
        #         name = "TaperedCouplerSupport",
        #         length = self.tapered_coupler_support_length,
        #         dy_min = self.tapered_coupler_support_width,
        #         dy_max=3*self.tapered_coupler_support_width)

        # Create Device_References of the Devices
        spacer1_ref = self.add_ref(spacer_dev)
        spacer2_ref = self.add_ref(spacer_dev)
        spacer3_ref = self.add_ref(spacer_dev)
        spacer4_ref = self.add_ref(spacer_dev)
        cavity_ref = self.add_ref(cavity_dev)
        tapered_support1_ref = self.add_ref(tapered_support_L)
        tapered_support2_ref = self.add_ref(tapered_support_L)
        tapered_support3_ref = self.add_ref(tapered_support_R)
        tapered_support4_ref = self.add_ref(tapered_support_R)
        tapered_coupler1_ref = self.add_ref(tapered_coupler)
        # tapered_couplerSupport1_ref = self.add_ref(couplerSupport)
        tapered_coupler2_ref = self.add_ref(tapered_coupler)
        # tapered_couplerSupport2_ref = self.add_ref(couplerSupport)

        # Connect up the references to build the geometry
        # tapered_coupler1_ref.connect('tpport1',tapered_couplerSupport1_ref.ports['tpport1'])
        spacer1_ref.connect('in_spacer', tapered_coupler1_ref.ports['tpport2'])
        tapered_support1_ref.connect('tpport1',
                                     spacer1_ref.ports['out_spacer'])
        spacer2_ref.connect('in_spacer', tapered_support1_ref.ports['tpport2'])
        tapered_support2_ref.connect('tpport1',
                                     spacer2_ref.ports['out_spacer'])
        cavity_ref.connect('in_cavity', tapered_support2_ref.ports['tpport2'])
        tapered_support3_ref.connect('tpport1', cavity_ref.ports['out_cavity'])
        spacer3_ref.connect('in_spacer', tapered_support3_ref.ports['tpport2'])
        tapered_support4_ref.connect('tpport1',
                                     spacer3_ref.ports['out_spacer'])
        spacer4_ref.connect('in_spacer', tapered_support4_ref.ports['tpport2'])
        tapered_coupler2_ref.connect('tpport2',
                                     spacer4_ref.ports['out_spacer'])
        # tapered_couplerSupport2_ref.connect('tpport1',tapered_coupler2_ref.ports['tpport1'])
        # Shift center of bounding box to origin.
        self.center = [0, 0]


class AirholePCC_PeriodOnly(SividdleDevice):
    """Airholes for PCC which only modulate period.
    These are similar to V0p4p2 Devices, but include the right most
    airhole of the wvg-mirr region which was dropped for some reason 
    in previous designs.

    Parameters
    ----------
    params: dict
        Dictionary containing the parameters of the photonic crystal cavity
    params['layer']: int
        Target layer for photonic crystal cavity holes
    params['aL']: float
        Lattice constant for the left hand side mirror
    params['aR']: float
        Lattice constant for the right hand side mirror
    params['hxL']: float
        Hole dimension along the waveguide long axis
        for the left hand side mirror.
    params['hyL']: float
        Hole dimension orthogonal to the waveguide long axis
        for the left hand side mirror.
    params['hxR']: float
        Hole dimension along the waveguide long axis
        for the right hand side mirror.
    params['hyR']: float
        Hole dimension orthogonal to the waveguide long axis
        for the right hand side mirror.
    params['maxDef']: float
        cavity defect parameter
    params['nholesLMirror']: int
        Number of holes in the bragg mirror zone on the LHS
    params['nholesRMirror']: int
        Number of holes in the bragg mirror zone on the RHS
    params['nholes_wvg-mirr_trans_L']: int
        Number of holes transitioning from mirror dimensions
        to unperforated waveguide on the LHS.
    params['nholes_wvg-mirr_trans_R']: int
        Number of holes transitioning from mirror dimensions
        to unperforated waveguide on the RHS.
    params['nholes_defect']: int
        Number of holes in each half of the defect region.
        That is, the total number of holes in the cavity zone
        will be nholes_defect*2.
    params['min_hole_dim']: float
        minimum fabable hole dimension
    params['effective_index']: float
        Effective refractive index of the waveguide as determined
        by a mode source in lumerical. This is used to calculate
        the standing wave "lattice constant" that the mirror-to-
        waveguide lattice constant transitions to.
    params['resonance_wavelength']: float
        resonance wavelength of the cavity in um. This is used to calculate
        the standing wave "lattice constant" that the mirror-to-
        waveguide lattice constant transitions to.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='PCC')
        # Breakout params to help code read cleaner
        self.hxL = params['hxL']
        self.hyL = params['hyL']
        self.hxR = params['hxR']
        self.hyR = params['hyR']
        self.nholesLMirror = params['nholesLMirror']
        self.nholes_wvgmirr_trans_L = params['nholes_wvg-mirr_trans_L']
        self.nholes_wvgmirr_trans_R = params['nholes_wvg-mirr_trans_R']
        self.nholesRMirror = params['nholesRMirror']
        self.min_hole_dim = params['min_hole_dim']

        self.aL = params['aL']
        self.aR = params['aR']

        self.nholes_defect = params['nholes_defect']
        self.maxDef = params['maxDef']

        self.neff = params['effective_index']
        self.res_wavelen = params['resonance_wavelength']

        self.layer = params['layer']

        self.params = params

        self.all_hx, self.all_hy, self.all_a = self.compute_geometry()
        print(self.all_hx, self.all_hy, self.all_a)
        print(self.all_hx.shape, self.all_hy.shape, self.all_a.shape)
        self << EllipseArray(self.layer, self.all_hx, self.all_hy, self.all_a)
        # Shift center of bounding box to origin.
        self.center = [0, 0]

    def _defect(self,a,i,n,maxDef):
        return (a*(1-maxDef*(2*((i-1)/n)**3 - 3*((i-1)/n)**2 + 1)))

    def compute_geometry(self):
        # waveguide-mirror hole size transition from smallest fab-able to
        # mirror hole size while maintaining ellipse aspect ratio.
        leftest_hole_hx = ((self.hxL / self.hyL) * (self.hxL > self.hyL)
                           + (self.hxL <= self.hyL)) * self.min_hole_dim
        lhs_wvg_mirror_hx = np.linspace(leftest_hole_hx, self.hxL,
                                        self.params['nholes_wvg-mirr_trans_L'])
        leftest_hole_hy = ((self.hyL / self.hxL) * (self.hxL < self.hyL)
                           + (self.hxL >= self.hyL)) * self.min_hole_dim
        lhs_wvg_mirror_hy = np.linspace(leftest_hole_hy, self.hyL,
                                        self.params['nholes_wvg-mirr_trans_L'])

        # LHS Photonic crystal mirror holes
        lhs_mirror_hx = self.hxL * np.ones(self.nholesLMirror)
        lhs_mirror_hy = self.hyL * np.ones(self.nholesLMirror)

        # defect holes
        # In general, the mirror hole dimensions may be different on each side,
        # and hole dims have even been varied as a cavity defect. However, it
        # is practically difficult to fabricate holes of varying sizes. This is
        # a constant hole size design, so no transition between hole sizes is
        # necessary
        lhs_cavity_hx = np.ones(self.nholes_defect) * self.hxL
        lhs_cavity_hy = np.ones(self.nholes_defect) * self.hyL
        rhs_cavity_hx = np.ones(self.nholes_defect) * self.hxR
        rhs_cavity_hy = np.ones(self.nholes_defect) * self.hyR

        rhs_mirror_hx = self.hxR * np.ones(self.nholesRMirror)
        rhs_mirror_hy = self.hyR * np.ones(self.nholesRMirror)

        # waveguide-mirror hole size transition from mirror hole size to
        # smallest fab-able dimension while maintaining ellipse aspect ratio.
        rightest_hole_hx = ((self.hxR / self.hyR) * (self.hxR > self.hyR)
                            + (self.hxR <= self.hyR)) * self.min_hole_dim
        rhs_wvg_mirror_hx = np.linspace(self.hxR, rightest_hole_hx,
                                        self.params['nholes_wvg-mirr_trans_R'])
        rightest_hole_hy = ((self.hyR / self.hxR) * (self.hxR < self.hyR)
                            + (self.hxR >= self.hyR)) * self.min_hole_dim
        rhs_wvg_mirror_hy = np.linspace(self.hyR, rightest_hole_hy,
                                        self.params['nholes_wvg-mirr_trans_R'])

        # append all holes together
        all_hx = np.append(lhs_wvg_mirror_hx,
                           np.append(lhs_mirror_hx, lhs_cavity_hx))
        all_hx = np.append(
            all_hx,
            np.append(rhs_cavity_hx,
                      np.append(rhs_mirror_hx, rhs_wvg_mirror_hx))
        )

        all_hy = np.append(lhs_wvg_mirror_hy,
                           np.append(lhs_mirror_hy, lhs_cavity_hy))
        all_hy = np.append(
            all_hy,
            np.append(rhs_cavity_hy,
                      np.append(rhs_mirror_hy, rhs_wvg_mirror_hy))
        )

        # lattice constant calculations

        # standing wave lattice constant based on effective index and
        # resonance wavelength
        standing_wave_a = 0.5 * self.res_wavelen / self.neff
        leftest_a = standing_wave_a
        # linearly transition the lattice constant from the waveguide's
        # standing wave spacing to the mirror spacing. The cartoon heuristic
        # is that this should help keep the field mostly in the dielectric,
        # and minimize scattering from surfaces.
        lhs_wvg_mirror_a = np.linspace(leftest_a, self.aL,
                                       self.nholes_wvgmirr_trans_L, endpoint=False)

        # LHS Photonic crystal mirror holes
        lhs_mirror_a = self.aL * np.ones(self.nholesLMirror)

        # defect holes
        # linearly transition the lattice constant from aL to aR
        a = np.linspace(self.aL,self.aR,2*self.nholes_defect-1)

        lhs = np.linspace(self.nholes_defect,1,self.nholes_defect)
        rhs = np.linspace(2,self.nholes_defect,self.nholes_defect-1)
        lhs_cavity_a = self._defect(a[:self.nholes_defect],lhs,self.nholes_defect,self.maxDef)
        rhs_cavity_a = self._defect(a[self.nholes_defect:],rhs,self.nholes_defect,self.maxDef)

        # RHS Photonic crystal mirror holes
        rhs_mirror_a = self.aR * np.ones(self.nholesRMirror)

        # RHS mirror-to-waveguide transition
        rightest_a = standing_wave_a
        # have to make rhs_wvg_mirror_a in a kind of run around way
        # to get the endpoints right
        rhs_wvg_mirror_a = np.linspace(rightest_a, self.aR,
                                       self.nholes_wvgmirr_trans_R, endpoint=False)[::-1]

        # append all holes together
        all_a = np.append(lhs_wvg_mirror_a,
                          np.append(lhs_mirror_a, lhs_cavity_a))
        all_a = np.append(all_a,
                          np.append(rhs_cavity_a,
                                    np.append(rhs_mirror_a, rhs_wvg_mirror_a)))

        return all_hx, all_hy, all_a




class OvercoupledPCCv0p4p2(SividdleDevice):
    """Airholes for V0p4p2 Devices.

    Parameters
    ----------
    params: dict
        Dictionary containing the parameters of the photonic crystal cavity
    params['layer']: int
        Target layer for photonic crystal cavity holes
    params['aL']: float
        Lattice constant for the left hand side mirror
    params['aR']: float
        Lattice constant for the right hand side mirror
    params['hxL']: float
        Hole dimension along the waveguide long axis
        for the left hand side mirror.
    params['hyL']: float
        Hole dimension orthogonal to the waveguide long axis
        for the left hand side mirror.
    params['hxR']: float
        Hole dimension along the waveguide long axis
        for the right hand side mirror.
    params['hyR']: float
        Hole dimension orthogonal to the waveguide long axis
        for the right hand side mirror.
    params['maxDef']: float
        cavity defect parameter
    params['nholesLMirror']: int
        Number of holes in the bragg mirror zone on the LHS
    params['nholesRMirror']: int
        Number of holes in the bragg mirror zone on the RHS
    params['nholes_wvg-mirr_trans_L']: int
        Number of holes transitioning from mirror dimensions
        to unperforated waveguide on the LHS.
    params['nholes_wvg-mirr_trans_R']: int
        Number of holes transitioning from mirror dimensions
        to unperforated waveguide on the RHS.
    params['nholes_defect']: int
        Number of holes in each half of the defect region.
        That is, the total number of holes in the cavity zone
        will be nholes_defect*2.
    params['min_hole_dim']: float
        minimum fabable hole dimension
    params['effective_index']: float
        Effective refractive index of the waveguide as determined
        by a mode source in lumerical. This is used to calculate
        the standing wave "lattice constant" that the mirror-to-
        waveguide lattice constant transitions to.
    params['resonance_wavelength']: float
        resonance wavelength of the cavity in um. This is used to calculate
        the standing wave "lattice constant" that the mirror-to-
        waveguide lattice constant transitions to.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='PCC')
        # Breakout params to help code read cleaner
        self.hxL = params['hxL']
        self.hyL = params['hyL']
        self.hxR = params['hxR']
        self.hyR = params['hyR']
        self.nholesLMirror = params['nholesLMirror']
        self.nholes_wvgmirr_trans_L = params['nholes_wvg-mirr_trans_L']
        self.nholes_wvgmirr_trans_R = params['nholes_wvg-mirr_trans_R']
        self.nholesRMirror = params['nholesRMirror']
        self.min_hole_dim = params['min_hole_dim']

        self.aL = params['aL']
        self.aR = params['aR']

        self.nholes_defect = params['nholes_defect']
        self.maxDef = params['maxDef']

        self.neff = params['effective_index']
        self.res_wavelen = params['resonance_wavelength']

        self.layer = params['layer']

        self.params = params

        self.all_hx, self.all_hy, self.all_a = self.compute_geometry()
        print(self.all_hx, self.all_hy, self.all_a)
        print(self.all_hx.shape, self.all_hy.shape, self.all_a.shape)
        self << EllipseArray(self.layer, self.all_hx, self.all_hy, self.all_a)
        # Shift center of bounding box to origin.
        self.center = [0, 0]

    def compute_geometry(self):
        # waveguide-mirror hole size transition from smallest fab-able to
        # mirror hole size while maintaining ellipse aspect ratio.
        leftest_hole_hx = ((self.hxL / self.hyL) * (self.hxL > self.hyL)
                           + (self.hxL <= self.hyL)) * self.min_hole_dim
        lhs_wvg_mirror_hx = np.linspace(leftest_hole_hx, self.hxL,
                                        self.params['nholes_wvg-mirr_trans_L'])
        leftest_hole_hy = ((self.hyL / self.hxL) * (self.hxL < self.hyL)
                           + (self.hxL >= self.hyL)) * self.min_hole_dim
        lhs_wvg_mirror_hy = np.linspace(leftest_hole_hy, self.hyL,
                                        self.params['nholes_wvg-mirr_trans_L'])

        # LHS Photonic crystal mirror holes
        lhs_mirror_hx = self.hxL * np.ones(self.nholesLMirror)
        lhs_mirror_hy = self.hyL * np.ones(self.nholesLMirror)

        # defect holes
        # In general, the mirror hole dimensions may be different on each side,
        # and hole dims have even been varied as a cavity defect. However, it
        # is practically difficult to fabricate holes of varying sizes. This is
        # a constant hole size design, so no transition between hole sizes is
        # necessary
        lhs_cavity_hx = np.ones(self.nholes_defect) * self.hxL
        lhs_cavity_hy = np.ones(self.nholes_defect) * self.hyL
        rhs_cavity_hx = np.ones(self.nholes_defect) * self.hxR
        rhs_cavity_hy = np.ones(self.nholes_defect) * self.hyR

        rhs_mirror_hx = self.hxR * np.ones(self.nholesRMirror)
        rhs_mirror_hy = self.hyR * np.ones(self.nholesRMirror)

        # waveguide-mirror hole size transition from mirror hole size to
        # smallest fab-able dimension while maintaining ellipse aspect ratio.
        rightest_hole_hx = ((self.hxR / self.hyR) * (self.hxR > self.hyR)
                            + (self.hxR <= self.hyR)) * self.min_hole_dim
        rhs_wvg_mirror_hx = np.linspace(self.hxR, rightest_hole_hx,
                                        self.params['nholes_wvg-mirr_trans_R'])
        rightest_hole_hy = ((self.hyR / self.hxR) * (self.hxR < self.hyR)
                            + (self.hxR >= self.hyR)) * self.min_hole_dim
        rhs_wvg_mirror_hy = np.linspace(self.hyR, rightest_hole_hy,
                                        self.params['nholes_wvg-mirr_trans_R'])

        # append all holes together
        all_hx = np.append(lhs_wvg_mirror_hx,
                           np.append(lhs_mirror_hx, lhs_cavity_hx))
        all_hx = np.append(
            all_hx,
            np.append(rhs_cavity_hx,
                      np.append(rhs_mirror_hx, rhs_wvg_mirror_hx))
        )

        all_hy = np.append(lhs_wvg_mirror_hy,
                           np.append(lhs_mirror_hy, lhs_cavity_hy))
        all_hy = np.append(
            all_hy,
            np.append(rhs_cavity_hy,
                      np.append(rhs_mirror_hy, rhs_wvg_mirror_hy))
        )

        # lattice constant calculations

        # standing wave lattice constant based on effective index and
        # resonance wavelength
        standing_wave_a = 0.5 * self.res_wavelen / self.neff
        leftest_a = standing_wave_a
        # linearly transition the lattice constant from the waveguide's
        # standing wave spacing to the mirror spacing. The cartoon heuristic
        # is that this should help keep the field mostly in the dielectric,
        # and minimize scattering from surfaces.
        lhs_wvg_mirror_a = np.linspace(leftest_a, self.aL,
                                       self.nholes_wvgmirr_trans_L + 1)

        # LHS Photonic crystal mirror holes
        lhs_mirror_a = self.aL * np.ones(self.nholesLMirror - 1)

        # defect holes
        # linearly transition the lattice constant from aL to aR
        a_trans_slope = (self.aL - self.aR) / (2 * self.nholes_defect - 1)
        lhs_cavity_a = np.linspace(1, self.nholes_defect, self.nholes_defect)
        #
        a_l_trans = self.aL + (
            lhs_cavity_a - self.nholes_defect) * a_trans_slope
        lhs_cavity_a[0] = a_l_trans[0] * (1 - self.maxDef)
        lhs_cavity_a[1:] = a_l_trans[1:] * (
            1 - self.maxDef * (2 * (
                (lhs_cavity_a[1:] - 1) / self.nholes_defect)**3 - 3 * (
                (lhs_cavity_a[1:] - 1) / self.nholes_defect)**2 + 1)
        )
        lhs_cavity_a = lhs_cavity_a[::-1]
        # blindly copying the form of the defect from the matlab codebase

        rhs_cavity_a = np.linspace(1, self.nholes_defect, self.nholes_defect)
        # This is almost certainly wrong
        a_r_trans = self.aR - (
            rhs_cavity_a - self.nholes_defect) * a_trans_slope
        rhs_cavity_a[0] = a_r_trans[0] * (1 - self.maxDef)
        rhs_cavity_a[1:] = a_r_trans[1:] * (1 - self.maxDef * (
            2 * ((rhs_cavity_a[1:] - 1) / self.nholes_defect)**3 - 3 * (
                (rhs_cavity_a[1:] - 1) / self.nholes_defect)**2 + 1
        ))
        # Turns out the simulated cavity actually didn't use this smallest
        # 'a' at rhs_cavity_a[0]. So we remove it here.
        rhs_cavity_a = rhs_cavity_a[1:]

        # RHS Photonic crystal mirror holes
        rhs_mirror_a = self.aR * np.ones(self.nholesRMirror - 1)

        # RHS mirror-to-waveguide transition
        rightest_a = standing_wave_a
        rhs_wvg_mirror_a = np.linspace(self.aR, rightest_a,
                                       self.nholes_wvgmirr_trans_R + 1)

        # append all holes together
        all_a = np.append(lhs_wvg_mirror_a,
                          np.append(lhs_mirror_a, lhs_cavity_a))
        all_a = np.append(all_a,
                          np.append(rhs_cavity_a,
                                    np.append(rhs_mirror_a, rhs_wvg_mirror_a)))

        return all_hx, all_hy, all_a

class AirholeDevice(SividdleDevice):
    def __init__(self, pcc_params, dt_params, scaling, coupler_supports=False):

        SividdleDevice.__init__(self, name='FreeGeom_Device')

        self.pcc_params = pcc_params.copy()
        self.dt_params = dt_params.copy()
        self.scaling = scaling
        print("scaling = {},hxL = {}".format(scaling, self.pcc_params['hxL']))
        self.pcc_params['aL'] *= self.scaling
        self.pcc_params['aR'] *= self.scaling
        self.pcc_params['hxL'] *= self.scaling
        self.pcc_params['hyL'] *= self.scaling
        self.pcc_params['hxR'] *= self.scaling
        self.pcc_params['hyR'] *= self.scaling

        self.dt_params['tapered_support_width'] *= self.scaling
        self.dt_params['width'] *= self.scaling

        device_holes = AirholePCC_PeriodOnly(self.pcc_params)
        if(coupler_supports):
            device_waveguide = DoubleTaperedDeviceWithCouplerSupports(
                               self.dt_params)
        else:
            device_waveguide = DoubleTaperedDevice(self.dt_params)

        

        self << device_holes
        self << device_waveguide
        # Shift center of bounding box to origin.
        self.center = [0, 0]


class OvercoupledAirholeDevicev0p4p2(SividdleDevice):
    def __init__(self, pcc_params, dt_params, scaling):

        SividdleDevice.__init__(self, name='FreeGeom_Device')

        self.pcc_params = pcc_params.copy()
        self.dt_params = dt_params.copy()
        self.scaling = scaling
        print("scaling = {},hxL = {}".format(scaling, self.pcc_params['hxL']))
        self.pcc_params['aL'] *= self.scaling
        self.pcc_params['aR'] *= self.scaling
        self.pcc_params['hxL'] *= self.scaling
        self.pcc_params['hyL'] *= self.scaling
        self.pcc_params['hxR'] *= self.scaling
        self.pcc_params['hyR'] *= self.scaling

        self.dt_params['tapered_support_width'] *= self.scaling
        self.dt_params['width'] *= self.scaling

        device_holes = OvercoupledPCCv0p4p2(self.pcc_params)
        device_waveguide = DoubleTaperedDevice(self.dt_params)

        self << device_holes
        self << device_waveguide
        # Shift center of bounding box to origin.
        self.center = [0, 0]


class OvercoupledAirholeDeviceWSupportv0p4p2(SividdleDevice):
    def __init__(self, pcc_params, dt_params, scaling):

        SividdleDevice.__init__(self, name='FreeGeom_Device')

        self.pcc_params = pcc_params.copy()
        self.dt_params = dt_params.copy()
        self.scaling = scaling
        print("scaling = {},hxL = {}".format(scaling, self.pcc_params['hxL']))
        self.pcc_params['aL'] *= self.scaling
        self.pcc_params['aR'] *= self.scaling
        self.pcc_params['hxL'] *= self.scaling
        self.pcc_params['hyL'] *= self.scaling
        self.pcc_params['hxR'] *= self.scaling
        self.pcc_params['hyR'] *= self.scaling

        self.dt_params['tapered_support_width'] *= self.scaling
        self.dt_params['width'] *= self.scaling

        device_holes = OvercoupledPCCv0p4p2(self.pcc_params)
        device_waveguide = DoubleTaperedDeviceWithCouplerSupports(
            self.dt_params)

        self << device_holes
        self << device_waveguide
        # Shift center of bounding box to origin.
        self.center = [0, 0]


class ImplantationWindow(SividdleDevice):
    """This is pretty much just a wrapper for gdspy's round."""

    def __init__(self, width, height, layer):
        SividdleDevice.__init__(self, name='Implant_Window')
        self.width = width
        self.height = height
        self.layer = layer
        implant_window = pg.rectangle(
            size=(width, height),
            layer=layer
        )
    #    implantWindow = gdspy.Round(
    #             (0, 0),
    #             [width * 0.5, height * 0.5],
    #             tolerance=1e-4,
    #             layer=layer
    #         )
        self.add(implant_window)
        self.center = [0, 0]


class RetroReflector(SividdleDevice):
    """Device containing retroreflector.

    Parameters
    ----------
    params: dict
        Dictionary containing the parameters of the retroreflector.
    params['len_wg']: float
        Length of waveguide section.
    params['height_wg']: float
        Height of waveguide section.
    params['len_tp_left']: float
        Lenght of left tapered section.
    params['dx_anchor']: float
        Distance of anchoring point to end of waveguide.
    params['width_anchor']: float
        Width of anchor at anchoring point.
    params['widthmax_anchor']: float
        Width of anchor opposite of anchoring point.
    params['length_anchor']: float
        Length of anchor.
    params['invert']: boolean
        If True, invert layer by bounding box.
        This is useful if we're working with positive
        photoresists.
    params['invert_layer']: float
        Target layer for inverted design.
        hx_init: float
        Initial value of hx.
    params['hy_init']: float
        Initial value of hx.
    params['hy_final']: float
        Final value of hy.
    params['a_const']: float:
        Lattice constant.
    params['num_taper']: int
        Number of tapered cells.
    params['num_cells']: int
        Number of non-tapered cells.
    params['dx_holes']: float
        Distance from holes to start of
        waveguide.
    params['which_anchors']: list
        Which anchor to add, labelled:
        1 ..... 2
        3 ..... 4
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='Retroreflector')

        params['add_taper_marker'] = False
        params['len_tp_right'] = 0

        waveguide = TaperedWaveGuide(params)

        self << waveguide

        # Carry over ports
        self.add_port(
            port=waveguide.ports['tpwgport1'],
            name='retroport1'
        )

        self.add_port(
            port=waveguide.ports['tpwgport2'],
            name='retroport2'
        )


class InterPoserRetroReflector(SividdleDevice):
    """Tapered waveguide - retroreflector combination.

    Parameters
    ----------
    params: dict
        Device parameter.

    For keys, see docstring of RetroReflector
    and TaperedWaveGuide.

    """

    def __init__(self, params):

        SividdleDevice.__init__(
            self,
            name='tapered_waveguide_with_retroreflector'
        )

        # Re-wrap photonic crystal parameters
        photonic_crystal_params = {
            'hx_init'                           : params['hx_init'],
            'hy_init'                           : params['hy_init'],
            'hy_final'                          : params['hy_final'],
            'a_const'                           : params['a_const'],
            'num_taper'                         : params['num_taper'],
            'num_cells'                         : params['num_cells'],
            'dx_holes'                          : params['dx_holes'],
            'holes_layer'                       : params['holes_layer'],
        }

        # Adjust params array for waveguide
        params['taper_gap_right'] = 0
        params['which_anchors'] = params['which_anchors_waveguide']
        params['len_tp_left'] = params['len_tp_left_waveguide']
        params['len_tp_right'] = params['len_tp_right_waveguide']
        params['photonic_crystal_params'] = None
        params['dx_anchor'] = params['dx_anchor_waveguide']
        params['len_wg'] = params['len_wg_waveguide']
        params['width_anchor'] = params['width_anchor_waveguide']
        params['widthmax_anchor'] = params['widthmax_anchor_waveguide']
        params['height_wg'] = params['height_waveguide']

        # Make waveguide
        waveguide = TaperedWaveGuide(params)

        # Adjust params array for retroreflector
        params['taper_gap_right'] = 0
        params['taper_gap_left'] = 0
        params['which_anchors'] = params['which_anchors_retroreflector']
        params['len_tp_left'] = params['len_tp_left_retroreflector']
        params['dx_anchor'] = params['dx_anchor_retroreflector']
        params['photonic_crystal_params'] = photonic_crystal_params
        params['len_wg'] = params['len_wg_retroreflector']
        params['width_anchor'] = params['width_anchor_retroreflector']
        params['widthmax_anchor'] = params['widthmax_anchor_retroreflector']
        params['height_wg'] = params['height_retroreflector']

        retro = RetroReflector(params)

        waveguide_on_device = self << waveguide
        retro_on_device = self << retro.movex(
            -(waveguide.xsize - retro.xsize) * 0.5
            + waveguide.xsize + params['dist_conn']
        )

        # connect the ports
        waveguide_on_device.connect(
            port='tpwgport2',
            destination=retro_on_device.ports['retroport1']
        )

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class RectangularSweep(SividdleDevice):
    """Lays out equidistant grid of devices generated using different parameters.

    Parameters
    ----------
    sweep_params: dict
        Dictionary containing the settings of sweep,
        with the following keys.
    sweep_params['params']: dict
        Setting dictionary for device class to be used.
    sweep_params['device_name']: string
        Name used for gds cell.
    sweep_params['device_class']: SividdleDevice
        Sividdle device to be replicated
        and for which the params dictionary is furnished.
    sweep_params['varsx']: array
        Array of values to be changed in x-direction,
        will replace params at key 'keyx'.
    sweep_params['keyx']: string
        Dictionary key of entry in params for which
        varsx furnishes the new variables.
    sweep_params['varxz'] /  sweep_params['keyy']: array / string
        Same for y-direction.
    sweep_params['pitchx'] / sweep_params['pitchy']: float
        Separation of bounding boxes of different
        devices in x-direction / y-direction.
    sweep_params['equidistant_grid']: boolean
       If true, adjust pitch in order to take device dimensions into account.
    sweep_params['grid_label']: boolean
        If True, label grid by assigning each
        device a coordinate 'A0', 'A1', etc.
    sweep_params['grid_label_params']: dict
        Parameters of grid label.
    sweep_params['grid_label_params']['label_layer']: int
        Layer of grid label to be exposed.
    sweep_params['grid_label_params']['textsize']: int
        Textsize of label.
    sweep_params['grid_label_params']['font']: string
        Font style,
        from here:
        https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/fonts_demo.html
    sweep_params['grid_label_params']['label_dist']: float
        Distance between label and device.
    sweep_params['grid_label_params']['revert_numbers']: boolean
        Revert ordering of numbers
        (sometimes usefull if same sweep
        is replicated in mirrored way).
    sweep_params['grid_label_params']['revert_letters']: boolean
        Revert odering of letters.
    sweep_params['staggered']: boolean
        If True, stagger columns in y direction.
    sweep_params['staggerd_y_pitch']: float
        Pitch in y-direction applied it staggerd_y_pitch
        is set to True.
    """

    def __init__(self, sweep_params):

        # Shorten
        sp = copy.deepcopy(sweep_params)

        SividdleDevice.__init__(self, name=sp['sweep_name'])

        num_iter_x = len(sp['varsx'])
        num_iter_y = len(sp['varsy'])

        # keep track of necessary paddings
        padding_x = np.zeros([num_iter_x, num_iter_y])
        padding_y = np.zeros_like(padding_x)

        # Stores dimensions of devices
        device_dimensions = np.zeros([num_iter_x, num_iter_y, 2])

        # Generate labels
        letter_label = list(string.ascii_uppercase)[0:num_iter_x]
        number_label = [str(i) for i in range(num_iter_y)]

        if sp['grid_label_params']['revert_numbers']:
            number_label = number_label[::-1]

        if sp['grid_label_params']['revert_letters']:
            letter_label = letter_label[::-1]

        # make devices
        # TODO: Adjust for non quadratic grid

        for i, j in it.product(range(num_iter_x), range(num_iter_y)):
            sp['device_params'][sp['keyx']] = sp['varsx'][i]
            sp['device_params'][sp['keyy']] = sp['varsy'][j]
            sp['device_params']['id_string'] = '{}{}'.format(
                letter_label[i],
                number_label[j]
            )

            device = sp['device_class'](sp['device_params'])
            device_dimensions[i, j, 0] = device.xsize
            device_dimensions[i, j, 1] = device.ysize

        for i, j in it.product(range(num_iter_x), range(num_iter_y)):

            # Refresh copy in case key entries have been altered
            sp = copy.deepcopy(sweep_params)

            sp['device_params'][sp['keyx']] = sp['varsx'][i]
            sp['device_params'][sp['keyy']] = sp['varsy'][j]
            sp['device_params']['id_string'] = '{}{}'.format(
                letter_label[i],
                number_label[j]
            )

            new_device = sp['device_class'](sp['device_params'])

            if i == 0 and sp['staggered']:
                padding_y[i, j] = sp['staggerd_y_pitch'] * j

            if j < num_iter_y - 1:
                # One to the right, takes into account xsize.
                current_xsize = device_dimensions[i, j, 0]
                right_xsize = device_dimensions[i, j + 1, 0]
                # Only adjust padding if equidistant_grid = True

                padding_x[i, j + 1] = padding_x[i, j] + \
                    (current_xsize + right_xsize) * 0.5 * \
                    sp['equidistant_grid'] + sp['pitchx']

            if i < num_iter_x - 1:
                current_ysize = device_dimensions[i, j, 1]
                top_ysize = device_dimensions[i + 1, j, 1]
                # Only adjust padding if equidistant_grid = True
                padding_y[i + 1, j] = padding_y[i, j] + \
                    (current_ysize + top_ysize) * 0.5 * \
                    sp['equidistant_grid'] + sp['pitchy']

            # Add grid labels.
            if sp['grid_label']:

                if i == 0 or i == num_iter_x - 1:
                    sp['grid_label_params']['text'] = number_label[j]
                    if i == 0:
                        sp['grid_label_params']['orientation'] = 'b'
                    else:
                        sp['grid_label_params']['orientation'] = 't'

                    new_device.add_device_label(sp['grid_label_params'])
                if j == 0 or j == num_iter_y - 1:
                    sp['grid_label_params']['text'] = letter_label[i]
                    if j == 0:
                        sp['grid_label_params']['orientation'] = 'l'
                    else:
                        sp['grid_label_params']['orientation'] = 'r'

                    new_device.add_device_label(sp['grid_label_params'])

            self << new_device.move([padding_x[i, j], padding_y[i, j]])

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class ImageArray(SividdleDevice):
    """Transforms a color image to BW array and generates an according pattern.

    Parameters
    ----------
    params: dict
        Dictionary containing setting
    params['name']: string
        Name of GDS cell.
    params['image']: string
        Path to image
    params['threshold']: int
        Threshold from 0 - 255 separating black from white. Not used
        if picture dithering is enabled.
    params['pixel_size']: int
        Physical size of one pixel on the design in um.
    params['layer']: int
        Layer where picture will be displayed.
    params['image_device']: SividdleDevice
        Device which is used as pixel. If none, a rectangle is used.
    params['dither']: Boolean
        If true, picture dithering is used to approximate greyscale.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        # Generate binary bitmap out of image
        bitmap = image_to_binary_bitmap(
            params['image'],
            params['threshold'],
            params['dither']
        )
        x_image = bitmap.shape[0]
        y_image = bitmap.shape[1]

        if params['image_device'] is None:
            # Define pixel polygon
            pixel = pg.rectangle(
                size=(
                    params['pixel_size'],
                    params['pixel_size']
                ),
                layer=params['layer']
            )
        else:
            pixel = params['image_device']

        for x in range(x_image):
            for y in range(y_image):
                if bitmap[x, y] == 1:
                    reference = self.add_ref(pixel)
                    reference.move(
                        [
                            x * pixel.xsize,
                            y * pixel.ysize
                        ]
                    )

        # Shift center of bounding box to origin.
        self.center = [0, 0]


class RenderedText(SividdleDevice):
    """Device containing rendered text.

    Parameters
    ----------
    params: dict
        Dictionary containing all settings
    params['name']: string
        Name of gds cell.
    params['text']: string
        Text to render.
    params['style']: string
        Font style,
        from here:
        https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/fonts_demo.html
    params['fontsize']: float
        Font size to render.
    params['layer']:
        Layer of rendered text.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        font_prop = FontProperties(style=params['style'])

        text = gdspy.PolygonSet(
            render_text(
                params['text'],
                params['fontsize'],
                font_prop=font_prop
            ),
            layer=params['layer']
        )

        self.add(text)

        # Shift center of bounding box to origin.
        self.center = [0, 0]
