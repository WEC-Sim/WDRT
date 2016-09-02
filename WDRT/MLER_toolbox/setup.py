from setuptools import setup, find_packages
 
setup(
    name = "WDRT-MLER",
    version = "0.1",
    packages = find_packages(),

    # metadata for upload to PyPI
    author = "Eliot Quon",
    author_email = "eliot.quon@nrel.gov",
    description = "WEC Design Response Toolbox (WDRT) implementation of the Most Likely Extreme Response (MLER) method",
    license = "Apache License, Version 2.0",
    url = "http://wec-sim.github.io/WDRT/index.html",
)
