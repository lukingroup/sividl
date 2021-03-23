# Bring your packages onto the path
import os
import sys

import numpy as np
sys.path.append(os.path.abspath(os.path.join('.')))
import sividl.sividl_devices as sivp  # noqa: E402

# ==============================================================================
# Adapted from example.py for production design of OvercoupledAirhole v0p4p1
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometers.
# ==============================================================================


def run_example():
    # Setup writefield.
    # Alignment mark parameters
    alignment_mark_params = {
        'layer'               : 1,
        'exposure_box'        : True,
        'exposure_box_dx'     : 8,
        'exposure_box_layer'  : 12,
        'invert'              : False,
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
        'bounding_box_layer'    : 25,
        'positive'              : True,
        'alignment_layer'       : 1,
        'alignment_offset_dx'   : 235,
        'alignment_offset_dy'   : 235,
        'exposure_box'          : True,
        'exposure_box_dx'       : 8,
        'exposure_box_layer'    : 12,
        'add_text_label'        : True,
        'text_label_layer'      : 200,
        'alignment_mark_params' : alignment_mark_params
    }

    # Generate Writefield.
    write_field = sivp.WriteFieldCrossAligmentMark(writefield_parameters)

    waveguide_width = 0.493
    #scale the design width for fab
    waveguide_width_HSQ = 1.1*waveguide_width
    
    nominalResonance = 0.7377
    targetResonance = 0.737
    resonance_scaling = targetResonance/nominalResonance


    PCC_params = {
        'layer'               : 2,
        'aL'                  : 0.2717,
        'aR'                  : 0.2502,
        'hxL'                 : 0.1135849,
        'hyL'                 : 0.1605274,
        'hxR'                 : 0.1135849,
        'hyR'                 : 0.1605274,
        'maxDef'              : 0.1392,
        'nholesLMirror'       : 7,
        'nholesRMirror'       : 3,
        'nholes_wvg-mirr_trans_L': 5,
        'nholes_wvg-mirr_trans_R': 5,
        'nholes_defect'       : 5,
        'min_hole_dim'        : 0.05,
        'effective_index'     : 1.6,
        'resonance_wavelength': 0.737
    }

    DT_params = {
        'cavity_length'         : 15.0,
        'layer_wg'              : 1,
        'name'                  : "DT",
        'tapered_coupler_length': 17.5,
        'tapered_support_length': 10,
        'tapered_coupler_minWidth':0.11,
        'tapered_support_width' : 1.4*waveguide_width_HSQ,
        'waveguide_spacer_length': 6,
        'width'               : waveguide_width_HSQ
    }

    DT_wSupport_params = {
        'cavity_length'         : 15.0,
        'layer_wg'              : 1,
        'name'                  : "DT_wSupport",
        'tapered_coupler_length': 17.5,
        'tapered_support_length': 10,
        'tapered_coupler_minWidth': 0.11,
        'tapered_coupler_support_beam_length' : 6.0,
        'tapered_coupler_support_beam_width'  : 0.33,
        'tapered_coupler_support_roof_height' : 1.5,
        'tapered_coupler_support_house_width': 1.1,
        'tapered_coupler_support_house_length': 1.5,
        'tapered_support_width' : 1.4*waveguide_width_HSQ,
        'waveguide_spacer_length': 6,
        'width'               : waveguide_width_HSQ
    }


    num_cols = 5
    num_rows = 15
    wf_width = writefield_parameters['bounding_box_size']
    deviceLength = 2*DT_params['tapered_coupler_length']+\
                    4*DT_params['tapered_support_length']+\
                        4*DT_params['waveguide_spacer_length']+\
                            DT_params['cavity_length']
    print(deviceLength)
    margin_y = 30.0
    margin_large = deviceLength/2+10.0    
    offset = 6.0
    photonic_scaling = np.linspace(0.965,1.035,num_cols)
    fab_scaling = resonance_scaling*photonic_scaling
    print(fab_scaling)
    for i,x in enumerate(np.linspace(-wf_width/2+margin_large,wf_width/2-margin_large,num_cols)):
        text_params = {
            'name'                  : 'label_{:1.5}'.format(photonic_scaling[i]),
            'text'                  : '{:1.5}'.format(photonic_scaling[i]),
            'style'                 : 'normal',
            'fontsize'              : 5,
            'layer'                 : 1
        }
        
        text_label = sivp.RenderedText(text_params)
        label_top_ref = write_field.add_ref(text_label)
        label_top_ref.move([x,wf_width/2-margin_y/2])
        label_bot_ref = write_field.add_ref(text_label)
        label_bot_ref.move([x,-wf_width/2+margin_y/2])
    
        # x += (i>num_cols/2)*margin_small
        v0p4p2_dev0 = sivp.OvercoupledAirholeDevice_v0p4p2(PCC_params, DT_params, fab_scaling[i])
        v0p4p2_dev1 = sivp.OvercoupledAirholeDevice_wSupport_v0p4p2(PCC_params, DT_wSupport_params, fab_scaling[i])
        # if(i==int(num_cols/2)):
        #     pass
        # else:
        #     for y in np.linspace(wf_width/2-margin_small-offset*i,-wf_width/2+margin_small+offset*(num_cols-i),num_rows):
        #         dev_ref = write_field.add_ref(v0p4p2_dev)
        #         dev_ref.move([x,y])
        for j,y in enumerate(np.linspace(wf_width/2-margin_y-offset*i,-wf_width/2+margin_y+offset*(num_cols-i-1),num_rows)):
            if(j%2==0):
                dev_ref = write_field.add_ref(v0p4p2_dev0)
            else:
                dev_ref = write_field.add_ref(v0p4p2_dev1)
            dev_ref.move([x,y])

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

    rarrow = sivp.RenderedText(rarrow_params)
    rarrow_top_ref = write_field.add_ref(rarrow)
    rarrow_bot_ref = write_field.add_ref(rarrow)
    rarrow_top_ref.move([200, 235])
    rarrow_bot_ref.move([200, -235])
    larrow = sivp.RenderedText(larrow_params)
    larrow_top_ref = write_field.add_ref(larrow)
    larrow_bot_ref = write_field.add_ref(larrow)
    larrow_top_ref.move([-200, 235])
    larrow_bot_ref.move([-200, -235])

    
    # # Export as GDS file.
    write_field.write_gds('Ovrcpld_v0p4p2.gds')


if __name__ == "__main__":
    run_example()


def test_run_example():
    """Pytest launcher for example test.

    This pytest will run the example and will fail if example
    code is throwing an error.

    TODO: Extend these tests.
    """
    run_example()
