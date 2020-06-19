# -*- coding: utf-8 -*-
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# make sure sphinx always uses the current branch
sys.path.insert(0, os.path.abspath('../..'))

# https://sphinx-rtd-theme.readthedocs.io/en/stable/
import sphinx_rtd_theme


# -- Project information -----------------------------------------------------

project = 'PROPKA 3'
copyright = '2020, Jan H. Jensen, Chresten R. Søndergaard, Mats H. M. Olsson, Michał Rostkowski, Nathan A. Baker, Matvey Adzhigirey, Oliver Beckstein, Jimmy Charnley Kromann, Mike Beachy, Toni G, Thomas Holder'
author = 'Jan H. Jensen, Chresten R. Søndergaard, Mats H. M. Olsson, Michał Rostkowski, Nathan A. Baker, Matvey Adzhigirey, Oliver Beckstein, Jimmy Charnley Kromann, Mike Beachy, Toni G, Thomas Holder'

# The full version, including alpha/beta/rc tags
release = '3.2.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.intersphinx',
              'sphinx.ext.mathjax', 'sphinx.ext.viewcode',
              'sphinx.ext.napoleon', 'sphinx.ext.todo',
              'sphinx.ext.autosummary',
##              'sphinx_sitemap',
              'sphinx_rtd_theme']

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

# for sitemap with https://github.com/jdillard/sphinx-sitemap
# change if we put it under a custom domain; for right now, assume RTD
##site_url = "https://propka.readthedocs.io"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# https://www.sphinx-doc.org/en/master/usage/extensions/autosummary.html
# We use a customized _templates/autosummary/module.rst to document members, too.
autosummary_generate = True
autosummary_imported_members = False
autosummary_generate_overwrite = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
