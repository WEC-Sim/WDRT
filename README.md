# WDRT documentation
## Accessing documentation
This directory contains the source for the WDRT's documentation, located at http://wec-sim.github.io/WDRT/.

## Updating documentation
The following procedure is similar to that used for the WEC-Sim users guide.

### Required packages
1. Python
1. Sphinx
  1. Install using ``pip install sphinx``
  1. Install sphinx extensions 
      1. ``pip install pip install sphinx_rtd_theme sphinxcontrib-bibtex''


### Edit and update html users guide
The users guide is developed using [Sphinx](http://sphinx-doc.org/) and rendered in html. To edit or add to the users guide, modify the source files located in the ``src`` folder using syntax and methods described in the [Sphinx Documentation](http://sphinx-doc.org/contents.html).
Once you are done editing, move to the base folder (this folder) and type ``make html`` from the command line to build the documentation.
This builds a html formatted copy of the documentation.
After building the HTML users guide, you can view the local copy of the documentation by opening the ``index.html`` file in a web browser


### Update the documentation on the http://wec-sim.github.io/wdrt website
The github.io website renders the documentation in the ``gh-pages`` branch as a website located at http://wec-sim.github.io/wdrt.
Here are the steps to do this in a terminal:

  ```Shell
  # Create a gh-pages branch in your repository
  git checkout --orphan gh-pages
  
  # Remove all files from the repository
  git rm -rf .

  # pull the gh-pages branch from master
  git pull origin gh-pages
  
  # push to the gh-pages branch on your fork
  git push -u remoteName gh-pages

  # Make your changes

  # Build the html documentation
  make html

  # Use git to check the status of the gh-pages branch, then commit and push changes. Once this step is performed the WEC-Sim website should be updated with any changes that were made to the source code.
  git status
  git add -A
  git commit -m 'update to documentation'
  git push

  # Move back to the master branch
  git checkout master
  ```
