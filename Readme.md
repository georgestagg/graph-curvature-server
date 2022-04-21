# Graph Curvature Calculator

This is the back end server for the [Graph Curvature Calculator](https://www.mas.ncl.ac.uk/graph-curvature/), a mathematical tool for calculating various notions of discrete graph curvature. Discrete curvature is an exciting new research area with possible applications to Bayesian networks, natural language processing, quantum gravity, and more.

If you have found this software useful, please cite the following article:
  * The Graph Curvature Calculator and the curvatures of cubic graphs, Experimental Mathematics, 2019 (arXiv:1712.03033 [math.CO])

## Requirements ##

 * Python 3+
 * numpy
 * scipy
 * sympy
 * web.py

## Installation Locally ##
To install the graph calculator locally, run the following commands in order.
 * `git clone git@github.com:georgestagg/graph-curvature-webapp.git`
 * `git clone git@github.com:georgestagg/graph-curvature-server.git`
 * `cd graph-curvature-server`
 * `python3 -m venv graph-curvature-venv`
 * `source graph-curvature-venv/bin/activate`
 * `pip install -r requirements.txt`
 * Start the server: `python graph.py 8090&`
 * `cd ../graph-curvature-webapp`
 * Start the webapp server: `python -m SimpleHTTPServer&`

Then visit `http://localhost:8000` in your browser.
