# -*- coding: utf-8 -*-

project = 'FEniCS hands-on'
author = 'Jan Blechta, Roland Herzog, Jaroslav Hron, Gerd Wachsmuth'
copyright = '2014, 2015, 2018 ' + author

# The short X.Y version
version = '2018.1'
# The full version, including alpha/beta/rc tags
release = '2018.1.0.dev0'

extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinxcontrib.contentui',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
html_theme = 'sphinx_rtd_theme'
html_theme_options = {}
html_static_path = ['_static']

latex_elements = {
    'papersize': 'a4paper',
    'fontsize': '10pt',
    #'preamble': '',
    #'figure_align': 'htbp',
}
latex_documents = [
    (master_doc, 'FEniCS-hands-on.tex', 'FEniCS hands-on',
     'Jan Blechta', 'howto'),
]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# Enable ".. only:: solution" sections
tags.add('solution')

# Add scroll bar to long equations
def setup(app):
    app.add_stylesheet('math-scroll.css')
