# paleointensity.org

# Introduction

Paleointensity.org is an online, open source, application to analyze paleointensity data produced by the most common paleointensity techniques.

# What can be found where?

Here we present the raw Python codes used for all calculations of the selection criteria and the graphing of the different methods. These codes can be tested with the test.py script. More about the python codes can be found in the “How to run the scripts” at the bottom of this readme document. 
The Python scripts are compiled to Javascript to run on paleointensity.org. 


# Authors

Annemarieke Béguin -- Utrecht University

Greig A. Paterson -- University of Liverpool

Andrew J. Biggin -- University of Liverpool

Lennart V. de Groot -- Utrecht University



# Brief Version History:

Beta version -- The beta version of Paleointensity.org is released alongside the submitted manuscript 


# How to run the scripts:
All the calculations are done in the Python scripts. There are multiple scripts and within a script multiple functions are defined.  

All functions are called from the “sc_main.py” script. Within “sc_main.py” four functions are defined, one for each of the methods: thermal_thellier, microwave_thellier, pseudo_thellier and MSP. The input for these four functions is the **_sample data_** information and the **_selection_** of this data used for the interpretation of the data. The output is organized in a so called “suitcase”. The structure of the suitcases is defined in “sc_config.py”. Each suitcase consists of multiple levels of dictionaries. There are three suitcases: the Thellier_suitcase (used for both Thermal and Microwave Thellier), the pseudo_Thellier_suitcase, and the MSP_suitcase. 

Within the four main functions a pipeline is defined, this is the order in which the different functions will be called. All scripts that are part of the suitcase pipeline are named “sc_XXX.py”. All functions within these scripts have the suitcase as input and use the information stored in the suitcase for the calculations, each function then returns the findings of the calculations in the correct places in the suitcase as defined in the “sc_config.py” script.

Additionally to the suitcase scripts there are scripts for the visualization of the data, graphing.py and graphing_pTh.py. Calculations that are used throughout the scripts are placed in helpers.py and specific helper functions for the Arai-plot curvature are in helpers_curve.py. The calculation of the MSP parameters is done outside the suitcase, in MSP_param.py, as these parameters do not depend on the data selection and therefore only have to be calculated once. The bootstrapped results of the MSP module are also calculated outside the suitcase and are done in MSP_bootsrap.py. 

For testing purposes, the “test.py” and corresponding “demo_data.py” files can be used. This demo data shows examples of samples and their selection for the interpretation, it follows the data structure of the website. There is a demo_data file for each method.


**Summary:** The data is interpreted using predefined suitcases that travel the different functions within the pipelines. These functions use the information in the suitcases for their calculations and return the findings to the same suitcase. 
