import gdspy
from matplotlib.textpath import TextPath
import numpy as np
from PIL import Image


# ==============================================================================
# sividl.utils -
# ==============================================================================
# This module helper function for sividl.devices
# ==============================================================================


# From https://gdspy.readthedocs.io/en/stable/gettingstarted.html
def render_text(text, size=None, position=(0, 0),
                font_prop=None, tolerance=0.1):
    path = TextPath(position, text, size=size, prop=font_prop)
    polys = []
    xmax = position[0]
    for points, code in path.iter_segments():
        if code == path.MOVETO:
            c = gdspy.Curve(*points, tolerance=tolerance)
        elif code == path.LINETO:
            c.L(*points)
        elif code == path.CURVE3:
            c.Q(*points)
        elif code == path.CURVE4:
            c.C(*points)
        elif code == path.CLOSEPOLY:
            poly = c.get_points()
            if poly.size > 0:
                if poly[:, 0].min() < xmax:
                    i = len(polys) - 1
                    while i >= 0:
                        if gdspy.inside(poly[:1], [polys[i]],
                                        precision=0.1 * tolerance)[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean([p], [poly], 'xor',
                                                 precision=0.1 * tolerance,
                                                 max_points=0).polygons[0]
                            break
                        elif gdspy.inside(polys[i][:1], [poly],
                                          precision=0.1 * tolerance)[0]:
                            p = polys.pop(i)
                            poly = gdspy.boolean([p], [poly], 'xor',
                                                 precision=0.1 * tolerance,
                                                 max_points=0).polygons[0]
                        i -= 1
                xmax = max(xmax, poly[:, 0].max())
                polys.append(poly)
    return polys


def image_to_binary_bitmap(filename, threshhold):
    """Convert image into binary array.

    From:
    https://stackoverflow.com/questions/50494541/bw-image-to-binary-array

    Parameters
    ----------
    filename: string
        Filename of the image.
    threshhold: int
        Value from 0-255 separating the white
        from the dark pixels.

    Return
    ----------
    bitmap: numpy.array
        Binary bitmap of the array.
    """
    img = Image.open(filename).convert('L')
    np_img = np.array(img)
    np_img = ~np_img  # invert B&W
    np_img[np_img <= threshhold] = 0
    np_img[np_img > threshhold] = 1

    return np_img
