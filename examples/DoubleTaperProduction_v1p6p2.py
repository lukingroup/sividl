import os
import sys

import numpy as np
sys.path.append(os.path.abspath(os.path.join('.')))
import phidl.geometry as pg
import sividl.sividl.sividl_devices as sivp

# ==============================================================================
# Adapted from example.py for production design of OvercoupledAirhole v0p4p2
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometers.
# ==============================================================================

FILENAME = "20231023_Ovrcpld_v1p6p2"

def run_example():
    waveguide_width = 0.48015718
    hsq_scaling_factor = 1.126
    waveguide_width_hsq = hsq_scaling_factor * waveguide_width

    # Scale the cavity design to hit target wavelegnth
    nominal_resonance = 737.69  # Simulated wavelength of cavity
    target_resonance = 737.134    # Target wavelength of cavity
    resonance_scaling = target_resonance / nominal_resonance

    aperture_size = 0.065
    aperture_offsets = np.array([0.7252, 0.6004, 0.7252, 0.6004,
                                 0.7252, 0.6004, 0.7252, 0.6004,
                                 0.7252, 0.6004, 0.7252, 0.6004,
                                 0.7252, 0.6004, 0.7252])
    aperture_ncols = 3
    aperture_nrows = 3
    aperture_aX    = 0.500
    aperture_aY    = 0.500

    pcc_params_7_5 = {
        'layer'                  : 2,
        'aL'                     : 0.27184455,
        'aR'                     : 0.24995467,
        'hxL'                    : 0.10824489,
        'hyL'                    : 0.16361557,
        'hxR'                    : 0.10824489,
        'hyR'                    : 0.16361557,
        'maxDef'                 : 0.14187774,
        'nholesLMirror'          : 7,
        'nholesRMirror'          : 5,
        'nholes_wvg-mirr_trans_L': 5,
        'nholes_wvg-mirr_trans_R': 5,
        'nholes_defect'          : 5,
        'min_hole_dim'           : 0.05,
        'effective_index'        : 1.6,
        'resonance_wavelength'   : 0.737134 # used to calculate standing wave "lattice constant"
        # that the mirror-to-waveguide lattice constant transitions to.
    }

    pcc_params_7_4 = {
        'layer'                  : 2,
        'aL'                     : 0.27184455,
        'aR'                     : 0.24995467,
        'hxL'                    : 0.10824489,
        'hyL'                    : 0.16361557,
        'hxR'                    : 0.10824489,
        'hyR'                    : 0.16361557,
        'maxDef'                 : 0.14187774,
        'nholesLMirror'          : 7,
        'nholesRMirror'          : 4,
        'nholes_wvg-mirr_trans_L': 5,
        'nholes_wvg-mirr_trans_R': 5,
        'nholes_defect'          : 5,
        'min_hole_dim'           : 0.05,
        'effective_index'        : 1.6,
        'resonance_wavelength'   : 0.737134 # used to calculate standing wave "lattice constant"
        # that the mirror-to-waveguide lattice constant transitions to.
    }

    double_taper_w_support_thermal_params = {
        'cavity_length'                       : 11.0,
        'layer_wg'                            : 1,
        'name'                                : "DT_wSupport",
        'tapered_coupler_length'              : 17,
        'tapered_support_length'              : 10,
        'tapered_coupler_minWidth'            : 0.11,
        'tapered_coupler_support_beam_length' : 6.0,
        'tapered_coupler_support_beam_width'  : 0.33,
        'tapered_coupler_support_roof_height' : 1.5,
        'tapered_coupler_support_house_width' : 1.1,
        'tapered_coupler_support_house_length': 1.5,
        'tapered_support_width' : np.array([1.8 * waveguide_width_hsq,
                                            1.4 * waveguide_width_hsq]),
        'waveguide_spacer_length'             : 6,
        'width'                               : waveguide_width_hsq
    }

    device_length = (
          2 * double_taper_w_support_thermal_params['tapered_coupler_length']
        + 2 * double_taper_w_support_thermal_params['tapered_coupler_support_beam_length']
        + 4 * double_taper_w_support_thermal_params['tapered_support_length']
        + 4 * double_taper_w_support_thermal_params['waveguide_spacer_length']
        + double_taper_w_support_thermal_params['cavity_length'])
    print(device_length)

    # Setup writefield.
    # Alignment mark parameters
    alignment_mark_params = {
        'layer'               : 1,
        'exposure_box'        : True,
        'exposure_box_dx'     : 8,
        'exposure_box_layer'  : 10,
        'invert'              : False,
        'd_small'             : 1.75,
        'd_large'             : 1.975,
        'sep'                 : 0.275,
        'make_dot'            : True,
        'dot_layer'           : 12,
        'dot_size'            : 0.010
    }

    # Setup writefield.
    writefield_parameters = {
        'bounding_box_size'     : 500,
        'bounding_box_layer'    : 25,
        'positive'              : True,
        'alignment_layer'       : 1,
        'alignment_offset_dx'   : 235,
        'alignment_offset_dy'   : 235,
        'exposure_box'          : True,
        'exposure_box_dx'       : 8,
        'exposure_box_layer'    : 10,
        'add_text_label'        : True,
        'text_label_layer'      : 199,
        'alignment_mark_params' : alignment_mark_params
    }

    # Generate Writefield.
    write_field = sivp.WriteFieldCrossAligmentMark(writefield_parameters)

    num_cols = 3
    num_rows = 15
    wf_width = writefield_parameters['bounding_box_size']

    margin_text = 30.0
    margin_y = 40.0
    margin_large = device_length / 2 + 10.0
    offset = (wf_width - 2 * margin_y) / num_rows / num_cols
    print("offset = ", offset)

    photonic_scaling = np.linspace(1.0, 0.965, num_cols)
    fab_scaling = resonance_scaling * photonic_scaling
    print(fab_scaling)

    col_coords = np.linspace(
        -wf_width / 2 + margin_large,
         wf_width / 2 - margin_large,
         num_cols)

    # Each column is a different photonic scaling to account for fab variance
    for i, x in enumerate(col_coords):
        text_params = {
            'name'     : f'label_{photonic_scaling[i]:1.5}',
            'text'     : f'{photonic_scaling[i]:1.5}',
            'style'    : 'normal',
            'fontsize' : 5,
            'layer'    : 5
        }

        text_label = sivp.RenderedText(text_params)
        label_top_ref = write_field.add_ref(text_label)
        label_top_ref.move([x, wf_width / 2 - margin_text / 2  - offset * i])
        label_bot_ref = write_field.add_ref(text_label)
        label_bot_ref.move([x, -wf_width / 2 + margin_text / 2 + offset * (num_cols - i - 1)])

        # Define the column of devices! They are defined here as they need to
        # have a differnt fab scaling per column.
        v0p4p2_7_5_dev1 = sivp.OvercoupledAirholeDeviceWSupportv0p4p2(
            pcc_params_7_5, double_taper_w_support_thermal_params,
            fab_scaling[i])
        v0p4p2_7_4_dev1 = sivp.OvercoupledAirholeDeviceWSupportv0p4p2(
            pcc_params_7_4, double_taper_w_support_thermal_params,
            fab_scaling[i])

        column_layout = [v0p4p2_7_4_dev1, v0p4p2_7_5_dev1, v0p4p2_7_4_dev1, v0p4p2_7_5_dev1,
                         v0p4p2_7_4_dev1, v0p4p2_7_5_dev1, v0p4p2_7_4_dev1, v0p4p2_7_5_dev1,
                         v0p4p2_7_4_dev1, v0p4p2_7_5_dev1, v0p4p2_7_4_dev1, v0p4p2_7_5_dev1,
                         v0p4p2_7_4_dev1, v0p4p2_7_5_dev1, v0p4p2_7_4_dev1]

        # TODO: why not fab_scaling?
        aperture_dev = sivp.ImplantationWindow(
            aperture_size * photonic_scaling[i],
            aperture_size * photonic_scaling[i], layer=11)
        aperture_offsets_scaled = aperture_offsets * photonic_scaling[i]
        aperture_aX_scaled = aperture_aX * photonic_scaling[i]
        aperture_aY_scaled = aperture_aY * photonic_scaling[i]

        row_coords = np.linspace(
             wf_width / 2 - margin_y - offset * i,
            -wf_width / 2 + margin_y + offset * (num_cols - i - 1),
            num_rows)

        # For each row, we just write the appropriate device and aperture
        # based on the column layout
        for j, y in enumerate(row_coords):
            dev_ref = write_field.add_ref(column_layout[j])
            dev_ref.move([x, y])

            ap_ref = write_field.add_array(aperture_dev,
                        columns=aperture_ncols, rows=aperture_nrows,
                        spacing=[aperture_aX_scaled, aperture_aY_scaled])
            ap_ref.move([x + aperture_offsets_scaled[j] - aperture_aX_scaled,
                         y - aperture_aY_scaled])

    # Add Arrow pointing to top right alignment marker
    rarrow_params = {
        'name'          : 'arrow',
        'text'          : '>',
        'fontsize'      : 40,
        'style'         : 'normal',
        'layer'         : 4,
    }

    larrow_params = {
        'name'          : 'arrow',
        'text'          : '<',
        'fontsize'      : 40,
        'style'         : 'normal',
        'layer'         : 4,
    }

    # Coarse alignment spots in apex of arrow
    dot_size = 0.01
    dot = pg.rectangle(size=(dot_size, dot_size), layer=10)

    rarrow = sivp.RenderedText(rarrow_params)
    rarrow_top_ref = write_field.add_ref(rarrow)
    rarrow_top_ref.move([200, 235])
    spot_top_ref = write_field.add_ref(dot)

    larrow = sivp.RenderedText(larrow_params)
    larrow_bot_ref = write_field.add_ref(larrow)
    larrow_bot_ref.move([-200, -235])
    spot_bot_ref = write_field.add_ref(dot)
    spot_bot_ref.move([-206.306 - dot_size/2, -235.022 - dot_size/2])

    # Alignment crosses at writefield corners
    alignment_cross = pg.cross(length=12.0, width=2.0, layer=1)

    cross_topleft_ref = write_field.add_ref(alignment_cross)
    cross_topleft_ref.move([-210, 235])
    cross_topright_ref = write_field.add_ref(alignment_cross)
    cross_topright_ref.move([210, 235])
    cross_botleft_ref = write_field.add_ref(alignment_cross)
    cross_botleft_ref.move([-210, -235])
    cross_botright_ref = write_field.add_ref(alignment_cross)
    cross_botright_ref.move([210, -235])
    write_field.move([10000, 10000])

    write_field.write_gds(f'{FILENAME}.gds', precision=1e-10)

if __name__ == "__main__":
    run_example()

def test_run_example():
    """Pytest launcher for example test. This pytest will run the example and
    will fail if example code is throwing an error.
    """
    run_example()
