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
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.rst']
pygments_style = 'sphinx'
highlight_language = 'python3'
html_theme = 'sphinx_rtd_theme'
html_theme_options = {}
html_static_path = ['_static']

latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    #'preamble': '',
    #'figure_align': 'htbp',
}
latex_documents = [
    (master_doc, 'FEniCS-hands-on.tex', 'FEniCS hands-on', author, 'howto'),
]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

## Enable ".. only:: solution" sections
#tags.add('solution')

# Fixes for typesetting math
def setup(app):
    app.add_stylesheet('math.css')
    app.add_javascript('copybutton.js')
