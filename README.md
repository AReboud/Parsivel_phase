# Parsivel_phase
Determine the phase of the hydrometeors from OTT Parsivel spectrum.

# Dependencies
The code is developed for python 3.9+ and should run on any recent Windows system (and most likely also Linux, but not tested).

The following python packages are required:

    numpy
    matplotlib
    seaborn
    h5py (to read the HDF files)
    pandas
    scipy
    os
    sklearn
    pydsd (available here https://github.com/josephhardinee/PyDSD)
    pytmatrix (available here https://github.com/jleinonen/pytmatrix)
    
Other useful functions are available in the func.py file.

# What does this code do?
This code processes the raw files (require 2 files .raw of OTT Parsivel2 records) to determine the phase of the precipitation at each time step.
Two methods are proposed:
 - Calculation of the precipitation phase according to the Synop Code, output from the Parsivel.
 - Calculation based from the particles size and velocity classes recorded by the Parsivel. The phase is infered through an experimental matrix of classification based on the schematic classification from Löffler-Mang (2000).

The major phase observed according those both methods are compared through a confusion matrix (with associated metrics) to determine the optimal parameters (snowflakes density, particles mass threshold, etc.).
Plots illustrating the reasonning are also provided.

# Data
Dataset samples are provided in the data folder (Parsivel file and spectrum). The dataset was recorded in Chamrousse (France) at 1750 m.a.s.l on December 2021.
The drop size distribution can be computed through the reader provided in Parsivel_ReaderCampbell.py (created by F. Cazenave, IRD) using the pydsd packages.

# Questions
In case of any questions, please don't hesitate to contact Arnaud Reboud: arnaud [dot] reboud [at] univ [dash] grenoble [dash] alpes [dot] fr
