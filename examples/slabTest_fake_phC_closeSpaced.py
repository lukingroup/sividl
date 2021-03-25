# Bring your packages onto the path
# import os
import sys
import numpy as np

sys.path.append('/Users/elizacornell/Documents/GitHub/sividl')

import sividl.sividl_devices as sivp

# ==============================================================================
# Demonstrates use of sividl.devices
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometers.
# ==============================================================================

alignment_mark_params = {
    'layer'               : 1,
    'exposure_box'     : True,
    'exposure_box_dx'     : 7,
    'exposure_box_layer'  : 12,
    'invert'              : True,
    'd_small'             : 1.75,
    'd_large'             : 1.975,
    'sep'                 : 0.275,
    'make_dot'            : True,
    'dot_layer'           : 13,
    'dot_size'            : 0.010
}

# Setup writefield.
writefield_parameters = {
    'bounding_box_size'     : 500,
    'bounding_box_layer'    : 255,
    'positive'              : True,
    'alignment_layer'       : 1,
    'alignment_offset_dx'   : 235,
    'alignment_offset_dy'   : 235,
    'add_text_label'        : True,
    'text_label_layer'      : 260,
    'alignment_mark_params' : alignment_mark_params
}

# Generate Writefield.
write_field = sivp.WriteFieldCrossAligmentMark(writefield_parameters)

# Initial slab parameters.
slab_cluster_parameters = {
    'expose_layer'  : 1,
    'num_slabs'     : 3,
    'lengths_slab'  : [20,20,20],
    'widths_slit'   : [5,5,5],
    'widths_slab'   : [2,2,2],
    'slab_dx'       : [2,4],
    'label_layer'   : 254
}

# Slab style defines whether the waveguide is defined by two windows around a slab, or by a single rectangle
phC_parameters = {
    'layer'         : 5,
    'length'        : 20,
    'height'        : 1,
    'slab_style'    : True,
    'window_width'  : 3,
    'photonic_crystal':True,
    'holes_layer'   : 5,
    'hx_init'       : 0.4,
    'hy_init'       : 0.4,
    'hy_final'      : 0.4,
    'a_const'       : 0.5,
    'num_taper'     : 2,
    'num_cells'     : 15,
    'dx_holes'      : 0,
    'both'          : True
#         'hx_array'      : 0.5*np.ones(5),
#         'hy_array'      : 0.5*np.ones(5), 
#         'a_consts_array': np.ones(5)
    # in this case width_slab = width_phC
#         'expose_layer'  : 5,
#         'length_slab'   : 20, 
#         'width_slit'    : 3,
#         'width_slab'    : 1,
#         'label_layer'   : 253,
}

tapered_waveguide_params = {
    'layer'                     : 6,
    'len_wg'                    : 20,
    'height_wg'                 : 1,
    'len_tp_left'               : 15,
    'len_tp_right'              : 15,
    'width_tp'                  : 0.2,
    'which_anchors'             : None,
    'dx_anchor'                 : 10,
    'width_anchor'              : 1.5,
    'widthmax_anchor'           : 2,
    'length_anchor'             : 3,
    'invert'                    : True,
    'invert_layer'              : 5,
    'taper_gap_left'            : 2,
    'taper_gap_right'           : 2,
    'photonic_crystal_params'   : phC_parameters
}

# Setup slab sweep parameters.
slabs_per_cluster = 3
num_devices_x = 3

width_slab_min = 2
width_slab_max = 2
num_iter_slab_widths = 3

width_slit_min = 5.5
width_slit_max = 8
num_iter_slit_widths = 6

# Setup phC sweep parameters.
#     num_holes_min = 6
#     num_holes_max = 11
#     num_iter_num_holes = 6

a_min = 0.5
a_max = 1
num_iter_a = 6

window_width_min = 1
window_width_max = 3.5
num_iter_window_widths = 3



# Arrays containing different parameters used for 2D sweep.
# slab_widths = np.linspace(
#     width_slab_min,
#     width_slab_max,
#     num_iter_slab_widths
# )

slit_width_coeffs = np.linspace(
    width_slit_min,
    width_slit_max,
    num_iter_slit_widths
)

slab_width_arrays = width_slab_min * np.ones((slabs_per_cluster,num_devices_x))
slit_width_arrays = np.outer(slit_width_coeffs, np.ones(slabs_per_cluster))


#     num_holes_vals = np.linspace(
#         num_holes_min,
#         num_holes_max,
#         num_iter_num_holes
#     )


a_vals = np.linspace(
    a_min,
    a_max,
    num_iter_a
)


window_widths = np.linspace(
    window_width_min,
    window_width_max,
    num_iter_window_widths
)

# We will add a label to the devices, with the following parameters:
grid_label_params = {
    'fontsize'      : 5,
    'style'         : 'normal',
    'layer'         : 1,
    'distance'      : 15,
    'revert_numbers': False,
    'revert_letters': False
}

# Sweep parameters fully specify the labelled 2D sweep.
slab_sweep_params = {
    'device_params'     : slab_cluster_parameters,
    'sweep_name'        : 'horizontal_slab_sweep',
    'device_class'      : sivp.EtchSlabCluster,
    'varsx'             : slit_width_arrays,
    'varsy'             : slab_width_arrays,
    'keyx'              : 'widths_slit',
    'keyy'              : 'widths_slab',
    'pitchx'            : 75,
    'pitchy'            : 30,
    'grid_label'        : False,
    'grid_label_params' : grid_label_params,
    'equidistant_grid'  : False,
    'staggered'         : False
}

waveguide_sweep_params = {
    'device_params'     : tapered_waveguide_params,
    'sweep_name'        : 'horizontal_wg_sweep',
    'device_class'      : sivp.TaperedWaveGuide,
    'varsx'             : a_vals,
    'varsy'             : window_widths,
    'keyx'              : 'a_const',
    'keyy'              : 'window_width',
    'pitchx'            : 75,
    'pitchy'            : 30,
    'grid_label'        : False,
    'grid_label_params' : grid_label_params,
    'equidistant_grid'  : False,
    'staggered'         : False
}

# Generate sweep in horizontal direction.
sweep_slab_horizontal = sivp.RectangularSweep(slab_sweep_params)
sweep_wg_horizontal = sivp.RectangularSweep(waveguide_sweep_params)

#     # Mirror device and sweep in vertical direction.

#     # Mirror sweep by reversing the order of sweep parameters.
slab_sweep_params['sweep_name'] = 'vertical_slab_sweep'
slab_sweep_params['varsx'] = slab_sweep_params['varsx'][::-1]
slab_sweep_params['varsy'] = slab_sweep_params['varsy'][::-1]

waveguide_sweep_params['sweep_name'] = 'vertical_phC_sweep'
waveguide_sweep_params['varsx'] = waveguide_sweep_params['varsx'][::-1]
waveguide_sweep_params['varsy'] = waveguide_sweep_params['varsy'][::-1]


#     # Reverse lettering order to match labelling of horizontal array.
grid_label_params['revert_numbers'] = True
grid_label_params['revert_letters'] = True
slab_sweep_params['grid_label_params'] = grid_label_params

#     # Generate sweep in horizontal direction.
sweep_slab_vertical = sivp.RectangularSweep(slab_sweep_params)
sweep_wg_vertical = sivp.RectangularSweep(waveguide_sweep_params)

#     # Add two sweeps to write field and move accordingly.
write_field << sweep_slab_horizontal.move([-95, +120])
write_field << sweep_wg_horizontal.move([-95, +120])
write_field << sweep_slab_vertical.rotate(45).move([+70, -70])
write_field << sweep_wg_vertical.rotate(45).move([+70, -70])

# Add Arrow pointing to top right alignment marker
arrow_params = {
    'name'          : 'arrow',
    'text'          : '→',
    'fontsize'      : 45,
    'style'         : 'normal',
    'layer'         : 4,
}

arrowH1 = sivp.RenderedText(arrow_params)
arrowH2 = sivp.RenderedText(arrow_params)
arrow_params['fontsize'] = 200
arrowV1 = sivp.RenderedText(arrow_params)
arrowV2 = sivp.RenderedText(arrow_params)



tapered_waveguide = sivp.TaperedWaveGuide(tapered_waveguide_params)


write_field << arrowH1.move([200, 235])
write_field << arrowH2.rotate(np.pi).move([-200, -235])
write_field << arrowV1.rotate(np.pi/2).move([200, 140])
write_field << arrowV2.rotate(-np.pi/2).move([-200, -140])


#     image = sivp.ImageArray(picture_parameters).move([-160, 160]).rotate(-90)
#     write_field << image

# Export as GDS file.
write_field.write_gds('slabTest_fake_phC_closeSpaced.gds')
