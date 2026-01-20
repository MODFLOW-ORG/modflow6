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
import shutil

# -- determine if running on readthedocs ------------------------------------
on_rtd = os.environ.get("READTHEDOCS") == "True"

# -- Copy files from .build_rtd_docs if not on readthedocs ------------------
if not on_rtd:
    src_pth = os.path.join("..", ".build_rtd_docs")

    # copy directories
    dirs = (
        "_mf6io",
        "_static",
    )
    for on_dir in dirs:
        src = os.path.join(src_pth, on_dir)
        dst = os.path.join(".", on_dir)
        if os.path.exists(dst):
            print(f"deleting...{dst}")
            shutil.rmtree(dst)
        print(f"copying {src} -> {dst}")
        shutil.copytree(src, dst)

    # copy files
    files = ("mf6io.rst", "index.rst")
    for file_name in files:
        src = os.path.join(src_pth, file_name)
        dst = os.path.join(".", file_name)
        if os.path.exists(dst):
            print(f"deleting...{dst}")
            os.remove(dst)
        print(f"copying {src} -> {dst}")
        shutil.copy(src, dst)


# -- Project information -----------------------------------------------------

project = "MODFLOW 6"
copyright = "2024, MODFLOW Development Team"
author = "MODFLOW Development Team"

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
    "myst_parser",
    "sphinx_markdown_tables",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Settings for GitHub actions integration
if on_rtd:
    extensions.append("rtds_action")
    rtds_action_github_repo = "MODFLOW-ORG/modflow6"
    rtds_action_path = "."
    rtds_action_artifact_prefix = "rtd-files-for-"
    rtds_action_github_token = os.environ.get("GITHUB_TOKEN", None)

# set master doc for readthedoce
master_doc = "index"

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
html_context = {
    "github_repo": "modflow6",
    "doc_path": ".doc",
}
html_css_files = [
    "_static/theme_overrides.css",  # override wide tables in RTD theme
]

html_theme_options = {}

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
