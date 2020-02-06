
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

import gdspy
from matplotlib.font_manager import FontProperties
import numpy as np
from phidl import Device
import phidl.geometry as pg
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

        # perform substraction for positive e-beam resist
        inverse = pg.boolean(
            A=bounding_box,
            B=self,
            operation='A-B',
            layer=layer
        )

        return inverse

    def add_label(self, params):
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


class CrossAligmentMark(SividdleDevice):
    """Write alignment marker.

    Parameters
    ----------
    positive (Boolean): Positive tone resist if true
    layer: int
        Layer of aligment mark.
    d_small: int
        Width of small rectangle.
    d_large:  int
        Width of large rectangle.
    sep: int
        Gap between rectangles.
    """

    def __init__(self, layer, d_small=1.75, d_large=1.975, sep=0.275):
        SividdleDevice.__init__(self, name='aligment_mark')
        self << pg.rectangle(size=(d_large, d_large), layer=layer)
        self << pg.rectangle(size=(d_small, d_small), layer=layer)\
            .movex(d_large + sep)
        self << pg.rectangle(size=(d_small, d_small), layer=layer)\
            .movey(d_large + sep)
        self << pg.rectangle(size=(d_large, d_large), layer=layer)\
            .move([d_small + sep, d_small + sep])


class WriteFieldCrossAligmentMark(SividdleDevice):
    """Writefiled with four cross-type alignment markers.

    Parameters
    ----------
    params: dict
        Contains the writefield parameters:
    params['bounding_box_size']: float
        Dimension of write field.
    params['positive']: boolean
        If True, pattern for positive tone resist is created.
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
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name='writefield')

        # make bounding box
        bounding_box = BoundingBox(
            params['bounding_box_layer'], params['bounding_box_size']
        )

        self << bounding_box

        # make alignment marks
        alignment_mark = CrossAligmentMark(params['alignment_layer'])

        if params['positive']:
            alignment_mark = alignment_mark.invert(params['alignment_layer'])

        # Position alignment marks on writefield
        delta_x = params['alignment_offset_dx']
        delta_y = params['alignment_offset_dy']

        self << pg.copy(alignment_mark).move(
            (
                bounding_box.xsize * 0.5 - delta_x,
                bounding_box.ysize * 0.5 - delta_y
            )
        )
        self << pg.copy(alignment_mark).move(
            (
                bounding_box.xsize * 0.5 - delta_x,
                -(bounding_box.ysize * 0.5 - delta_y)
            )
        )
        self << pg.copy(alignment_mark).move(
            (
                -(bounding_box.xsize * 0.5 - delta_x),
                bounding_box.ysize * 0.5 - delta_y
            )
        )
        self << pg.copy(alignment_mark).move(
            (
                -(bounding_box.xsize * 0.5 - delta_x),
                -(bounding_box.ysize * 0.5 - delta_y)
            )
        )

        # Shift center of bounding box to origin
        self.center = [0, 0]


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
        self.label(
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


class Taper(SividdleDevice):
    """Device describing a tapering section of a waveguide.

    This device will hace two ports associated left and right ends
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


class WaveGuide(SividdleDevice):
    """Device describing a rectangular waveguide.

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
    params['len_tp_rightt']: float
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

        if params['add_taper_marker']:
            self.add_arrow_markers()
            # Compensate for the arrows in assignment of center.
            self.center = [0, -params['width_tp'] * 0.5]
        else:
            self.center = [0, 0]

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

        assert len(a_consts_array) is len(hx_array) - 1, \
            "Please provide arrays of equal length"

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
    flip: boolen
        If True, flip array. Tapered section
        now is left.
    """

    def __init__(self, layer, hx_init, hy_init,
                 hy_final, a_const, num_taper, num_cells, flip):

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

        if flip:
            hx_array = np.flip(hx_array)
            hy_array = np.flip(hy_array)

        # Construct array of constant lattice consants.
        a_consts_array = np.ones(num_taper + num_cells - 1) * a_const

        self << EllipseArray(layer, hx_array, hy_array, a_consts_array)

        # Shift center of bounding box to origin.
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

                    new_device.add_label(sp['grid_label_params'])
                if j == 0 or j == num_iter_y - 1:
                    sp['grid_label_params']['text'] = letter_label[i]
                    if j == 0:
                        sp['grid_label_params']['orientation'] = 'l'
                    else:
                        sp['grid_label_params']['orientation'] = 'r'

                    new_device.add_label(sp['grid_label_params'])

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
        Threshold from 0 - 255 separating black from white
    params['pixel_size']: int
        Physical size of one pixel on the design in um.
    params['layer']: int
        Layer where picture will be displayed.
    """

    def __init__(self, params):

        SividdleDevice.__init__(self, name=params['name'])

        # Generate binary bitmap out of image
        bitmap = image_to_binary_bitmap(params['image'], params['threshold'])
        x_image = bitmap.shape[0]
        y_image = bitmap.shape[1]

        # Define pixel polygon
        pixel = pg.rectangle(
            size=(
                params['pixel_size'],
                params['pixel_size']
            ),
            layer=params['layer']
        )

        for x in range(x_image):
            for y in range(y_image):
                if bitmap[x, y] == 1:
                    self << pg.copy(pixel).move(
                        (
                            x * params['pixel_size'],
                            y * params['pixel_size']
                        )
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
