FEniCS hands-on
===============

This is the source code for building FEniCS hands-on tutorial.


Build on Read the Docs
----------------------

.. image:: https://readthedocs.org/projects/fenics-handson/badge/?version=latest
    :target: https://fenics-handson.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


Building locally
----------------

1. Install requirements::

       python3 -mpip install [--user|--prefix=...] -r requirements.txt

2. Build by one of::

       make
       make html
       make html O='-t solution'
       make html-solution  # equivalent to previous one
       make latexpdf
       make latexpdf O='-t solution'


Instant refresh in browser
--------------------------

To have automatic refresh in browser when developing
locally:

1. Install some auto reload plugin.
   For example: `Live Reload for Firefox
   <https://addons.mozilla.org/en-US/firefox/addon/live-reload>`_,
   set reload rule::

       Host URL: http://127.0.0.1:8000/*
       Source file URLs: http://127.0.0.1:8000/*

2. Run http server by::

       make servehtml
