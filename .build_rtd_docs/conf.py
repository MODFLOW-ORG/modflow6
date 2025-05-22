# ruff: noqa: E402
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import shutil

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys
from subprocess import PIPE, Popen

sys.path.insert(0, os.path.abspath(os.path.join("..", "doc")))
sys.path.insert(0, os.path.abspath(os.path.join("..", "distribution")))

# -- determine if running on readthedocs ------------------------------------
on_rtd = os.environ.get("READTHEDOCS") == "True"

# -- print current directory
print(f"Current Directory...'{os.path.abspath(os.getcwd())}'")

# -- clean up doxygen files -------------------------------------------------
dox_pths = ("_mf6io",)
for dox_pth in dox_pths:
    print(f"cleaning....{dox_pth}")
    for root, dirs, files in os.walk(dox_pth):
        for name in files:
            fpth = os.path.join(root, name)
            os.remove(fpth)

# -- Update the modflow 6 version -------------------------------------------
print("Update the modflow6 version")
from update_version import update_version

update_version()

# -- import version from doc/version.py -------------------------------------
from version import __version__

dstdir = "_mf6run"
if os.path.isdir(dstdir):
    shutil.rmtree(dstdir)
os.makedirs(dstdir)

print(f"Copy run-time comparison table to {dstdir}")
fpth = "run-time-comparison.md"
src = os.path.join("..", "distribution", fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

dstdir = "_dev"
if os.path.isdir(dstdir):
    shutil.rmtree(dstdir)
os.makedirs(dstdir)

print(f"Copy developer docs to {dstdir}")
fpth = "DEVELOPER.md"
src = os.path.join("..", fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

fpth = "CONTRIBUTING.md"
src = os.path.join("..", fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

fpth = "styleguide.md"
src = os.path.join(fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

fpth = "readme.md"
src = os.path.join("..", "doc", "mf6io", "mf6ivar", fpth)
dst = os.path.join(dstdir, "dfn.md")
shutil.copy(src, dst)

fpth = "EXTENDED.md"
src = os.path.join("..", fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

fpth = "CODE_OF_CONDUCT.md"
src = os.path.join("..", fpth)
dst = os.path.join(dstdir, fpth)
shutil.copy(src, dst)

dstdir = "_migration"
if os.path.isdir(dstdir):
    shutil.rmtree(dstdir)
os.makedirs(dstdir)

print(f"Copy migration guides to {dstdir}")
fpth = "mf6_6_0_prt_migration_guide.md"
src = os.path.join(fpth)
dst = os.path.join(dstdir, fpth)
shutil.move(src, dst)

# -- build the deprecations table --------------------------------------------
print("Build the deprecations markdown table")
pth = os.path.join("..", "doc", "mf6io", "mf6ivar")
args = (sys.executable, "deprecations.py")
# run the command
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:")
    print(stderr.decode("utf-8"))

# -- copy deprecations markdown ---------------------------------------------
print("Copy the deprecations table")
dstdir = "_mf6run"
fpth = "deprecations.md"
src = os.path.join("..", "doc", "mf6io", "mf6ivar", "md", fpth)
dst = os.path.join(dstdir, fpth)
# copy the file
shutil.copy(src, dst)

# -- build the mf6io markdown files -----------------------------------------
print("Build the mf6io markdown files")
pth = os.path.join("..", "doc", "mf6io", "mf6ivar")
args = (sys.executable, "mf6ivar.py")
# run the command
proc = Popen(args, stdout=PIPE, stderr=PIPE, cwd=pth)
stdout, stderr = proc.communicate()
if stdout:
    print(stdout.decode("utf-8"))
if stderr:
    print("Errors:")
    print(stderr.decode("utf-8"))

# -- update the doxygen version number ---------------------------------------
print("Update the Doxyfile with the latest version number")
with open("Doxyfile", "r") as fp:
    lines = fp.readlines()

tag = "PROJECT_NUMBER"
with open("Doxyfile", "w") as fp:
    for line in lines:
        if tag in line:
            line = f'{tag}         = "version {__version__}"\n'
        fp.write(line)

# -- Project information -----------------------------------------------------

project = "MODFLOW 6"
copyright = "2024, MODFLOW Development Team"
author = "MODFLOW Development Team"

# -- Project version ---------------------------------------------------------
version = __version__
release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "IPython.sphinxext.ipython_console_highlighting",  # lowercase didn't work
    "sphinx.ext.autosectionlabel",
    "nbsphinx",
    "nbsphinx_link",
    "myst_parser",
    "sphinx_markdown_tables",
    "sphinxcontrib.mermaid",
]

myst_fence_as_directive = ["mermaid"]

# # Tell sphinx what the pygments highlight language should be.
# highlight_language = 'fortran'
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = [
    "_static/theme_overrides.css",  # override wide tables in RTD theme
]

# html_theme_options = {
#     "github_url": "https://github.com/MODFLOW-ORG/modflow6",
#     "use_edit_page_button": False
# }

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True
