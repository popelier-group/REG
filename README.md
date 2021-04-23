# REG.py
This repository is a continuation of the development of the Relative Energy Gradient (REG) method implemented in python3.

The manual.pdf contains the information regarding the theory and the code implementation. 

Please, report bugs and issues to fabio.falcioni@manchester.ac.uk

# Usage
***Manual is not up to date but still useful for tips and tricks***

To run a REG analysis (with IQA and DFT-D3):
- Save this repository in your machine.
- Copy the `auto_reg.py` script in the folder where each REG step has been saved. Note: Each REG step should be saved as a numbered folder and contain the gaussian single point energy output, gaussian wavefunction output and the atomic-files folder obtained with AIMAll. 
- Enter the script, change the REG installation path and all the possible options as explained in the script.
- Run the command : 
`python3 auto_reg.py > reg.log &`
  
# Results
All results of the REG analysis will be saved in a folder called "SYSTEM_results" with SYSTEM being the chosen system name. 
Inside the folder different files for each type of energy will be found together with png images of the REG tables and various csv files for later data analysis.
The main file is called "REG.xlsx" which contains all the results together in one place.
