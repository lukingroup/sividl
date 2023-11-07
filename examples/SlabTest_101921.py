# Bring your packages onto the path
import os
import sys

import numpy as np
import phidl.geometry as pg

sys.path.append(os.path.abspath(os.path.join('.')))
import sividl.sividl_devices as sivp  # noqa: E402

# ==============================================================================
# Adapted from aperture array we used for kai-me but now 
# for e6 overgrowth project.
# ==============================================================================
# This code will generate a writefield and an array of labelled devices.
# Devices are two slits of varying width separated by varying distance.
#
# Note that all dimensions are in micrometers.
# ==============================================================================


def run_example():
    write_field = pg.import_gds(filename='slabTest_closeSpaced_staggered_smooth.gds',
                                cellname='toplevel', flatten=False)

    

    # aperture_sizes = np.linspace(5.5, 5.5, 1)
    aperture_sizes = [5.5]
    print(aperture_sizes)
    implant_windows = [sivp.ImplantationWindow(ap_dim, ap_dim, 5) 
                       for ap_dim in aperture_sizes]
    x0 = -122
    y0 = 180
    
    dx = 79
    dy = 33
    Dy = 1

    x = x0
    y = y0

    nrows = 12
    ncols = 4
    counter = 1
    for i in range(ncols):
        for j in range(nrows):
            if(counter % 5 == 0):
                aperture_ref = write_field.add_ref(implant_windows[0])
                aperture_ref.move([x, y])
                counter = 0
            counter +=  1
            y -= dy
        x += dx
        y = y0
        y += Dy * (i + 1)

    x = -102
    y = y0
    counter = 1
    for i in range(ncols):
        for j in range(nrows):
            if(counter % 5 == 0):
                aperture_ref = write_field.add_ref(implant_windows[0])
                aperture_ref.move([x, y])
                counter = 0
            counter +=  1
            y -= dy
        x += dx
        y = y0
        y += Dy * (i + 1)


    # add marker dots for coarse alignment
    marker_dot = sivp.ImplantationWindow(0.02, 0.02, 12)
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([-210.315, -235])
    marker_dot_ref = write_field.add_ref(marker_dot)
    marker_dot_ref.move([210.315, 235])

    write_field.move([250, 250])

    write_field.write_gds('SlabChop.gds')


if __name__ == "__main__":
    run_example()


def test_run_example():
    """Pytest launcher for example test.

    This pytest will run the example and will fail if example
    code is throwing an error.

    TODO: Extend these tests.
    """
    run_example()
