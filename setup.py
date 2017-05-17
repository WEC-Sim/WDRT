from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


setup(
    name = "WDRT",
    version = "1.0.0",
    url = "https://github.com/WEC-Sim/WDRT",
    packages=['WDRT', 'examples'],
)