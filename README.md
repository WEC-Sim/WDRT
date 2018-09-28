# WDRT documentation
## Accessing documentation
This directory contains the source for the WDRT's documentation, located at http://wec-sim.github.io/WDRT/.

## Updating documentation
The following procedure is similar to that used for the WEC-Sim users guide.

### Required packages
1. Python
1. Sphinx
  1. Install using ``pip install sphinx``
  1. Install the bibtex extension for sphinx using ``pip install sphinxcontrib-bibtex``
  1. Install the rtd theme for sphinx using ``pip install sphinx_rtd_theme``. You might have to manually move it to the ``sphinx/themes/`` directory.

### Edit and update html users guide
The users guide is developed using [Sphinx](http://sphinx-doc.org/) and rendered in html. To edit or add to the users guide, modify the source files located in the ``doc`` folder using syntax and methods described in the [Sphinx Documentation](http://sphinx-doc.org/contents.html).
Once you are done editing, move to the ``doc`` folder type ``make html`` from the command line to build the documentation.
This builds a html formatted copy of the documentation in the ``doc/_build/html`` folder.
After building the HTML users guide, you can view the local copy of the documentation by opening the ``doc/_build/html/index.html`` file in a web browser


### Update the documentation on the http://wec-sim.github.io/wdrt website
The github.io website renders the documentation in the ``gh-pages`` branch as a website located at http://wec-sim.github.io/wdrt.
The easiest way to update the website is to make the ``doc/_build/html`` folder a clone of ``gh-pages`` branch.
The user can then push changes in the html documentation directly to the ``gh-pages`` branch.
Here are the steps to do this in a Linux/Mac Terminal, note that windows instructions are very similar:

  ```Shell
  # Move to the _build directory
  cd doc/_build

  # Remove the html folder in the _build directory
  rm -rf html

  # Clone the gh-pages branch into a folder named html
  git clone --depth 1 -b gh-pages https://github.com/WEC-Sim/wdrt.git html

  # Move back to the users guide directory
  cd doc

  # Build the html documentation
  make html

  # Move doc/_build/html directory
  cd doc/_build/html

  # Use git to check the status of the gh-pages branch, then commit and push changes. Once this step is performed the WEC-Sim website should be updated with any changes that were made to the source code.
  git status
  git add -A
  git commit -m 'update to documentation'
  git push
  ```
