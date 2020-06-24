# THE COSMIC TRIGGER #

This project is part of a Bachelor Thesis at the Mu3e experiment @ Physikalisches Insitut Heidelberg.
Author: Konstantin Neureither
Supervisor: Prof. Dr. Andre Schöning

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

* Up till now, the project contains two executables
	* ```Mu3eRecAcc``` produces analysis and control plots
	* ```Mu3eCosPat``` runs a prototype of SuperPixel computation on root files of the simulation

Get it running:

* Running the CMake config should be straight forward, maybe the f2clib path needs to be changed in the top dir CMake file.

### Who do I talk to? ###

Contact: [neureither@physi.uni-heidelberg.de](mailto:neureither@physi.uni-heidelberg.de)