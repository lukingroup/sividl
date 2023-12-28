from distutils.core import setup 
# from setuptools import setup, find_packages #EK suggests to change to setuptools in future versions 23/12/27

setup(
  name = "sividl",
  packages = ["sividl"],
  py_modules=["sividl.sividl_devices", "sividl.sividl_utils"],
  version = "0.2",
  author = "Can Knaut, Erik Knall",
  author_email = "cknaut@g.harvard.edu",
  description = "1D Cavity Waveguide and Unit Cell gds generation",
)
