# Bring your packages onto the path
import os
import sys

import numpy as np
from math import floor
import phidl.geometry as pg

sys.path.append(os.path.abspath(os.path.join('.')))
import sividl.sividl_devices as sivp  # noqa: E402

# ==============================================================================
# Adapted from example.py for aperture array for NV implantation on Bart's failed COVID chips
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometers.
# ==============================================================================

def make_half_cavity(num_mirrs):
    adef = [252.49, 254.01, 257.94, 263.32, 269.17, 274.55, 278.48]
    amir = [280.0]
    atap = [274.98, 263.00, 248.70, 236.72, 231.70]
    abarrier = []
    for i in range(num_mirrs):
        abarrier += amir
    a = np.array(adef + abarrier + atap)
    a *= 1e-3
    a *= 0.965
    a[-3:-2] *= 0.93
    return np.copy(a)

def add_lhs(write_field, implant_window, lhs, x0, y0):
    xl = x0
    for ai in lhs:
        aperture_ref = write_field.add_ref(implant_window)
        aperture_ref.move([xl, y0])
        xl -= ai
    aperture_ref = write_field.add_ref(implant_window)
    aperture_ref.move([xl, y0])

def add_rhs(write_field, implant_window, rhs, x0, y0):
    xr = x0
    for ai in rhs:
        aperture_ref = write_field.add_ref(implant_window)
        aperture_ref.move([xr, y0])
        xr += ai
    aperture_ref = write_field.add_ref(implant_window)
    aperture_ref.move([xr, y0])

def run_example():

    write_field = pg.import_gds(filename='Lukinsiv91119.gds',
                                cellname='lukinsiv91119', flatten=False)
    # write_field.write_gds('BartProd_NVapertures.gds')
    # write_field = pg.import_gds(filename='BartProd_NVapertures.gds',
    #                             cellname='toplevel', flatten=False)
    
    implant_window = sivp.ImplantationWindow(0.062, 0.075, 11)

    half_cavity_7 = make_half_cavity(7)
    half_cavity_6 = make_half_cavity(6)
    half_cavity_5 = make_half_cavity(5)
    half_cavity_4 = make_half_cavity(4)
    half_cavity_3 = np.array([0.245, 0.245, 0.253, 0.253, 0.259, 0.2685, 0.2697, 0.274, 0.262, 0.270, 0.260, 0.2506, 0.230, 0.2317])

    mirr_num = np.array([[3, 6], [4, 6], [4, 6], [5, 7], [4, 4]])

    optical_scaling = [1.01, 1, 0.99, 0.98, 0.97, 0.96]

    x0 = 139
    y0 = 176

    y = y0

    for j in range(60):
        if (j > 19 and j < 40):
            x_shift = -15
        else:
            x_shift = 0
        # L and R Taper side of Cavity
        x0 = -165 + x_shift
        for x in np.linspace(x0 , x0 + 17.5, 50):
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([x, y])
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([-x, y])
        # L and R notch side of Cavity
        x0 = -129.5 + x_shift
        for x in np.linspace(x0 , x0 + 7, 20):
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([x, y])
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([-x, y])
        # CL and CR taper side of cavity
        x0 = -63 + x_shift
        for x in np.linspace(x0 , x0 + 10, 35):
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([x, y])
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([-x, y])
        # CL and CR notch side of Cavity
        x0 = -37 + x_shift
        for x in np.linspace(x0 , x0 + 7, 20):
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([x, y])
            aperture_ref = write_field.add_ref(implant_window)
            aperture_ref.move([-x, y])
        # L and R Cavities
        x0 = -139 + x_shift
        scaling = optical_scaling[floor(j/2 / 5)]
        add_lhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][0]), x0, y)
        add_rhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][1]), x0, y)
        # R
        add_rhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][0]), -x0, y)
        add_lhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][1]), -x0, y)
        # CL Cavities and CR
        x0 = -45 + x_shift
        scaling = optical_scaling[floor(j/2 / 5)]
        add_lhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][0]), x0, y)
        add_rhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][1]), x0, y)
        # CR
        add_rhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][0]), -x0, y)
        add_lhs(write_field, implant_window,
                scaling * make_half_cavity(mirr_num[int(j/2) % 5][1]), -x0, y)
        
        # y -= 13
        y-=6.5

    # add marker dots for fine alignment
    marker_dot = sivp.ImplantationWindow(0.02, 0.02, 11)
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([-155, -217.5])
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([155, -217.5])
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([-155, 207.5])
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([155, 207.5])

    # add marker dots for coarse alignment
    marker_dot = sivp.ImplantationWindow(0.02, 0.02, 10)
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([-155, -232.5])
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([155, -232.5])

    # add marker clear circles
    ellipse_params = {
        'rx'     : 8,
        'ry'     : 9,
        'layer'  : 10,
    }
    marker_clear = sivp.Ellipse(ellipse_params)

    marker_clear_ref = write_field.add_ref(marker_clear)
    marker_clear_ref.move([-155, -217.5])
    marker_clear_ref = write_field.add_ref(marker_clear)
    marker_clear_ref.move([155, -217.5])
    marker_clear_ref = write_field.add_ref(marker_clear)
    marker_clear_ref.move([-155, 207.5])
    marker_clear_ref = write_field.add_ref(marker_clear)
    marker_clear_ref.move([155, 207.5])

    write_field.write_gds('BartProd_NVapertures.gds')


if __name__ == "__main__":
    run_example()


def test_run_example():
    """Pytest launcher for example test.

    This pytest will run the example and will fail if example
    code is throwing an error.

    TODO: Extend these tests.
    """
    run_example()
