from distutils.core import setup

setup(
  name = "sividl",
  packages = ["sividl"],
  py_modules=["sividl.sividl_devices", "sividl.sividl_utils"],
  version = "0.2",
  author = "Can Knaut, Erik Knall",
  author_email = "cknaut@g.harvard.edu",
  description = "1D Cavity Waveguide and Unit Cell gds generation",
)