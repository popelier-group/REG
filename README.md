# REG.py

This repository is a continuation of the development of the Relative Energy Gradient (REG) method [[1]](#1) implemented in Python3.

The manual.pdf contains the information regarding the theory and the code implementation plus a full walk-through tutorial on a simple system.

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk

Please, cite [[1]](#1) if using this code for your studies/research. The code is under MIT license.
If using the Ramer-Douglas-Peucker algorithm on a potential energy surface from electronic structure calculations, please cite [[2]](#2)[[3]](#3).

# Installation

REG.py is installable locally using pip inside the cloned REG.py directory.
`python -m pip install -e .`

Note that _auto_reg.py_ is not part of the library but it is an external script used for automation.

# Usage

**_Manual has been update with a full tutorial. Functions section needs to be updated_**

To run a REG analysis (with IQA and DFT-D3):

- Install REG in your machine (preferably in a conda environment).
- Copy the `auto_reg.py` script in the folder where each REG step has been saved. Note: Each REG step should be saved as a numbered folder and contain the gaussian single point energy output, gaussian wavefunction output and the atomic-files folder obtained with AIMAll [[4]](#4).
- Enter the script with a text editor, change all the possible options as explained in the script.
- Run the command :
  `python3 auto_reg.py > reg.log &`

# Setup Ramer-Douglas-Peucker

The `setup_rdp.py` is a standalone script that can be run with specific options to find the suitable number of points on which to run the following REG-IQA analysis.
It uses the Ramer-Douglas-Peucker algorithm [[2]](#2)[[3]](#3) on the PES obtained from electronic structure calculations.

Run the command
`python3 setup_rdp.py -h (or --help)`
for the usage.

NOTE: if REG is installed with pip , this script is automatically added as an executable command named `rdp`.
Minimum requirements:

- A txt file containing the energies at each step of an IRC or PES scan (easily obtained in GaussView).

# Results

All results of the REG analysis will be saved in a folder called "SYSTEM_results" with SYSTEM being the chosen system name.
Inside the folder different files for each type of energy will be found together with png images of the REG tables and various csv files for later data analysis.
The main file is called "REG.xlsx" which contains all the results together in one place.
The "REG.xlsx" file can also be input into PyMol for better visualisation through a custom [PyMol-REG](https://github.com/popelier-group/PyMol-REG) plugin.

## References

<a id="1">[1]</a>
Thacker, Joseph CR, and Paul LA Popelier. "The ANANKE relative energy gradient (REG) method to automate IQA analysis over configurational change." Theoretical chemistry accounts 136.7 (2017): 1-13.

<a id="2">[2]</a>
Fabio Falcioni and Paul L. A. Popelier
Journal of Chemical Information and Modeling 2023 63 (14), 4312-4327
DOI: 10.1021/acs.jcim.3c00404

<a id="3">[3]</a>
Douglas, D.H. and T.K. Peucker, Algorithms for the reduction of the number of points required to represent a digitized line or its caricature. Cartographica: the international journal for geographic information and geovisualization, 1973. 10(2): p. 112-122.

<a id="4">[4]</a>
AIMAll (Version 19.10.12), Todd A. Keith, TK Gristmill Software, Overland Park KS, USA, 2019 (aim.tkgristmill.com)
