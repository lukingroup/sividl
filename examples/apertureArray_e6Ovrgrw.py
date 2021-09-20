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

    aperture_sizes = np.linspace(0.02, 0.07, 3)
    print(aperture_sizes)
    implant_windows = [sivp.ImplantationWindow(ap_dim, ap_dim, 11) for ap_dim in aperture_sizes]
    device_length = 80 # um
    x0 = 139
    y0 = 176

    y = y0 - 6.5
    num_apps = 0
    # in each row, len(aperture_sizes) groups of apertures. device_length/num_groups apertures per group
    # j is row index
    for j in range(30):
        if (j > 9 and j < 20):
            x_shift = -15
        else:
            x_shift = 0
        # L and R Cavities
        x0 = -184 + x_shift
        for group in range(len(aperture_sizes)):
            group_shift = (device_length / len(aperture_sizes)) * (group - len(aperture_sizes)/2)
            for x in np.linspace(x0 - group_shift, x0 - group_shift + (device_length / len(aperture_sizes)), int(device_length / len(aperture_sizes))):
                aperture_ref = write_field.add_ref(implant_windows[group])
                aperture_ref.move([-x, y])
                aperture_ref = write_field.add_ref(implant_windows[group])
                aperture_ref.move([x, y])
                num_apps += 2
        # CL Cavities and CR
        x0 = -90 + x_shift
        for group in range(len(aperture_sizes)):
            group_shift = (device_length / len(aperture_sizes)) * (group - len(aperture_sizes)/2)
            for x in np.linspace(x0 - group_shift, x0 - group_shift + (device_length / len(aperture_sizes)), int(device_length / len(aperture_sizes))):
                aperture_ref = write_field.add_ref(implant_windows[group])
                aperture_ref.move([-x, y])
                aperture_ref = write_field.add_ref(implant_windows[group])
                aperture_ref.move([x, y])
                num_apps += 2
        y -= 13
        # y-=6.5
    print("num apps = ", num_apps)

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

    write_field.write_gds('BartProd_e6ovrgrth.gds')


if __name__ == "__main__":
    run_example()


def test_run_example():
    """Pytest launcher for example test.

    This pytest will run the example and will fail if example
    code is throwing an error.

    TODO: Extend these tests.
    """
    run_example()
