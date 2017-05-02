### Requirements ###
* Python 2.7
* virtualenv for Python

### Installation Locally ###
To install the graph calculator locally, run the following commands in order.
* `git clone https://@mas-gitlab.ncl.ac.uk/graph-curvature/graph-curvature-webapp.git`
* `git clone https://@mas-gitlab.ncl.ac.uk/graph-curvature/graph-curvature-server.git`
* `cd graph-curvature-server`
* `virtualenv -p /usr/bin/python2 graph-curvature-venv`
* `source graph-curvature-venv/bin/activate`
* `pip install -r requirements.txt`
* Start the server: `python graph.py 8090&`
* `cd ../graph-curvature-webapp`
* Start the webapp server: `python -m SimpleHTTPServer&`

Then visit `http://localhost:8000` in your browser.

### Installation Online ###
To install the graph calculator online, run the following commands:
* `git clone https://@mas-gitlab.ncl.ac.uk/graph-curvature/graph-curvature-server.git`
* `cd graph-curvature-server`
* `virtualenv -p /usr/bin/python2 graph-curvature-venv`
* `source graph-curvature-venv/bin/activate`
* `pip install -r requirements.txt`
* Start the server: `python graph.py 8090&`
* Install the [graph-curvature-webapp](https://mas-gitlab.ncl.ac.uk/graph-curvature/graph-curvature-webapp) and follow the Installation Online instructions.

