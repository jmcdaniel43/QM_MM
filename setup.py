"""
setup.py: Builds repository and distribution information.
"""
import setuptools

__author__ = "Jesse McDaniel and John Pederson"
__version__ = "1.1.0"

with open("README.rst", "r") as fh:
    description = fh.read()

setuptools.setup(name="QM_MM",
                 version="1.1.0",
                 description="QM/MM with PME for long-range electrostatics",
                 author="Jesse McDaniel and John Pederson",
                 author_email="jpederson6@gatech.edu",
                 packages=["qm_mm"],
                 python_requires=">3.0")
