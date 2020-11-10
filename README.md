# THE COSMIC TRIGGER #

This project is part of a Bachelor Thesis at the Mu3e experiment @ Physikalisches Insitut Heidelberg.

*Author:* Konstantin Neureither

*Supervisor:* Prof. Dr. Andre Schöning

### What is this repository for? ###

* Analysis of Cosmic trigger data from the mu3e simulation package
* Streamline simulation data in order to improve evaluation efficiency
* Implemented helix fitting for cosmic muons
* Create superpixel and superpixel template mapping over detector area
* Evaluate and optimize template mapping
* Pattern Recognition of simulation data

### How do I get set up? ###

Requirements: 

* root package installed
* f2c libraries (fortran to c)

Components:

* The software contains several subcomponents, which are
	* ```CTPreSlimData``` -- slims the cosmic muon TriRec output from the Mu3e simulation into a CosmicData file, which contains the connected and fitted cosmic muon tracks. It can append several MC simulation runs to one Cosmic dataset file.
	* ```CTCoreModules``` -- this folder contains the ```Pattern Engine``` class that creates and handles super pixel mappings as well as the ```Template Bank```, which simulates an associative memory.
	* ```CTCosPatTrain``` -- scripts to build a database from a cosmic dataset. 
	* ```CTBgAna``` -- scripts to evaluate the false-positive trigger rate for non-cosmic hit frames. These scripts use a Mu3e simulation MC file as input and compute the full combinatorics in each frame, create Templates from it and test, if they occur in a given template bank.
	* ```CTPlottingScripts``` -- several scripts to create all different curves, such as ROC, SPC-T-Count fit, etc.
	
* Additional Folders
	* ```data``` a data directory is created where the template database files are stored after they are created.
	* ```util``` contains several functions that are used within several modules and scripts.
	* ```karimaki``` an implementation of a three dimensional helix fit, used in ```CTPreSlimData```.
	* ```root``` classes for ROOT file reading and writing.
	* ```Mu3eRecAcc``` produces analysis and control plots for a Mu3e simulation in cosmic mode.

Get it running:

* Running the CMake config should be straight forward, maybe the f2clib path needs to be changed in the top dir CMake file.
* The ```Configuration.h``` file in ```CTCoreModules``` is used to set the parameters. In a future upgrade, this will be replaced by actual configuration files - it is still a bit messy.

### Who do I talk to? ###

The software will be slightly improved and refractored in the next weeks, before it can be reused for further projects. If you need help, you can still contact me [neureither@physi.uni-heidelberg.de](mailto:neureither@physi.uni-heidelberg.de)