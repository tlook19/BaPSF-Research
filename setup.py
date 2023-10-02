import os
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

setup(
    name="bapsfda",
    version="0.1.3a4",
    description="An data analysis package for BaPSF/LAPD",
    author="Thomas R. Look",
    author_email="tlook@physics.ucla.edu",
    packages=["bapsfda", "bapsfda.analysis", "bapsfda.visualization"],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "astropy",
        "pycine",
        "scaleogram",
    ],
)
