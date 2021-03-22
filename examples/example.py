# Bring your packages onto the path
import os
import sys

import numpy as np
sys.path.append(os.path.abspath(os.path.join('.')))
import sividl.sividl_devices as sivp  # noqa: E402

# ==============================================================================
# Demonstrates use of sividl.devices
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometer.
# ==============================================================================


def run_example():
    # Setup writefield.
    # Alignment mark parameters
    alignment_mark_params = {
        'layer'               : 1,
        'exposure_box'        : True,
        'exposure_box_dx'     : 8,
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
        'exposure_box'          : True,
        'exposure_box_dx'       : 8,
        'exposure_box_layer'    : 12,
        'add_text_label'        : True,
        'text_label_layer'      : 260,
        'alignment_mark_params' : alignment_mark_params
    }

    # Generate Writefield.
    write_field = sivp.WriteFieldCrossAligmentMark(writefield_parameters)

    # Initial slab parameters.
    slab_parameters = {
        'expose_layer'  : 1,
        'length_slab'   : 20,
        'width_slit'    : 1,
        'width_slab'    : 2,
        'label_layer'   : 254,
    }

    # Setup sweep parameters.
    width_slab_min = 0.5
    width_slab_max = 7
    num_iter_slab_widths = 7

    width_slit_min = 0.5
    width_slit_max = 7
    num_iter_slit_widths = 7

    # Arrays containing different parameters used for 2D sweep.
    slab_widths = np.linspace(
        width_slab_min,
        width_slab_max,
        num_iter_slab_widths
    )

    slit_widths = np.linspace(
        width_slit_min,
        width_slit_max,
        num_iter_slit_widths
    )

    # We will add a label to the devices, with the following parameters:
    grid_label_params = {
        'fontsize'      : 5,
        'style'          : 'normal',
        'layer'         : 1,
        'distance'      : 15,
        'revert_numbers': False,
        'revert_letters': False
    }

    # Sweep parameters fully specify the labelled 2D sweep.
    sweep_params = {
        'device_params'     : slab_parameters,
        'sweep_name'        : 'horizontal_sweep',
        'device_class'      : sivp.EtchSlap,
        'varsx'             : slab_widths,
        'varsy'             : slit_widths,
        'keyx'              : 'width_slab',
        'keyy'              : 'width_slit',
        'pitchx'            : 30,
        'pitchy'            : 13,
        'grid_label'        : True,
        'grid_label_params' : grid_label_params,
        'equidistant_grid'  : True,
        'staggered'         : False
    }

    # Generate sweep in horizontal direction.
    sweep_box_horizontal = sivp.RectangularSweep(sweep_params)

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
    sweep_box_vertical = sivp.RectangularSweep(sweep_params)

    # Add two sweeps to write field and move accordingly.
    write_field << sweep_box_horizontal.move([-100, +105])
    write_field << sweep_box_vertical.rotate(45).move([+65, -80])

    # Add Arrow pointing to top right alignment marker
    arrow_params = {
        'name'          : 'arrow',
        'text'          : '→',
        'fontsize'      : 45,
        'style'         : 'normal',
        'layer'         : 4,
    }

    arrow = sivp.RenderedText(arrow_params)
    write_field << arrow.move([200, 227])

    # For fun, add a picture of the Harvard Logo to the empty space
    picture_parameters = {
        'name'          : 'image',
        'image'         : 'staticfiles/harvard_logo.jpeg',
        'threshold'     : 140,
        'dither'        : False,
        'image_device'  : None,
        'pixel_size'    : 0.5,
        'layer'         : 1
    }

    image = sivp.ImageArray(picture_parameters).move([-160, 160]).rotate(-90)
    write_field << image

    # Export as GDS file.
    write_field.write_gds('example.gds')     


if __name__ == "__main__":
    run_example()


def test_run_example():
    """Pytest launcher for example test.

    This pytest will run the example and throw will fail if example
    code is thothrowing an error.

    TODO: Extend these tests.
    """
    run_example()
