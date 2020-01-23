import numpy as np
import sividl_devices as sivp

# ==============================================================================
# Demonstrates use of sividl.devices
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometer.
# ==============================================================================


# Setup writefield.
writefield_parameters = {
    'bounding_box_size':   500,
    'bounding_box_layer':  255,
    'positive':            True,
    'alignment_layer':     1,
    'alignment_offset_dx': 25,
    'alignment_offset_dy': 25,
}

# Generate Writefield.
write_field = sivp.WriteFieldCrossAligmentMark(writefield_parameters)

# Initial slab parameters.
slab_parameters = {
    'expose_layer': 1,
    'length_slab':  20,
    'width_slit':   1,
    'width_slab':   2,
    'label_layer':  254,
}

# Setup sweep parameters.
width_slab_min = 0.5
width_slab_max = 7
num_iter_slab_widths = 7

width_slit_min = 0.5
width_slit_max = 7
num_iter_slit_widths = 7

# Arrays containing different parameters used for 2D sweep.
slab_widths = np.linspace(width_slab_min, width_slab_max, num_iter_slab_widths)
slit_widths = np.linspace(width_slit_min, width_slit_max, num_iter_slit_widths)

# We will add a label to the devices, with the following parameters:
grid_label_params = {
    'label_layer':      1,
    'textsize':         5,
    'font':             'normal',
    'label_dist':       15,
    'revert_numbers':   False,
    'revert_letters':   False
}

# Sweep parameters fully specify the labelled 2D sweep.
sweep_params = {
    'device_params':         slab_parameters,
    'sweep_name':           'horizontal_sweep',
    'device_class':         sivp.EtchSlap,
    'varsx':                slab_widths,
    'varsy':                slit_widths,
    'keyx':                 'width_slab',
    'keyy':                 'width_slit',
    'pitchx':               13,
    'pitchy':               13,
    'grid_label':           True,
    'grid_label_params':    grid_label_params
}

# Generate sweep in horizontal direction.
sweep_box_horizontal = sivp.EquidistantRectangularSweep(sweep_params)

# Mirror device and sweep in vertical direction.

# Mirror sweep by reversing the order of sweep parameters.
sweep_params['sweep_name'] = 'vertical_sweep'
sweep_params['varsx'] = sweep_params['varsx'][::-1]
sweep_params['varsy'] = sweep_params['varsy'][::-1]

# Reverse lettering order to match labelling of horizontal array.
grid_label_params['revert_numbers'] = True
grid_label_params['revert_letters'] = True
sweep_params['grid_label_params'] = grid_label_params

# Generate sweep in horizontal direction.
sweep_box_vertical = sivp.EquidistantRectangularSweep(sweep_params)

# Add two sweeps to write field and move accordingly.
write_field << sweep_box_horizontal.move([-100, +105])
write_field << sweep_box_vertical.rotate(45).move([+65, -80])

# Export as GDS file.
write_field.write_gds('example.gds')
