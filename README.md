# The Cosmic Trigger

This simulation software project was developed within a Bachelorthesis carried out at the *Physikalisches Institut Uni Heidelberg* by Konstantin Neureither under the supervision of Prof. Dr. André Schöning in 2020. 

Its main goal is the implementation of a *Pattern Recognition Hardware Trigger* designed to detect cosmic muon tracks in the detector data at the Mu3e experiment. As an input the *Cosmic Trigger* uses Monte-Carlo Simulation data produced by the Mu3e Simualation package (including data from the Cosmic Simulation, the Cosmic Reconstruction and the Beam Simulation). 


### Module overview

The evaluation is split into two main analysis chains, the **Database Training** (or Database Building) aswell as the **Background Evaluation** where the false-positive rate of the trigger is examined using "normal" beam decay hits.

The Associative Memory, that is used as a hardware database in the future implementation as well as the Pattern Management is implemented in the Core Modules, the **Pattern Engine** and the **Template Bank**. In order to create the patterns from data of the pixel detector, pixels are grouped to *Super Pixels*. For a cosmic muon track, the super pixels that were hit by a cosmic are combined to a *Template*. To detect cosmics in the detector data, the occurence of super pixel hits according to these templates is looked up in that data.


![Module Overview](img/Software_Module_Overview.png)


### Project File Structure

* The software contains several subcomponents, which are
    * ```CTPreSlimData``` -- slims the cosmic muon TriRec output from the Mu3e simulation into a CosmicData file, which contains the connected and fitted cosmic muon tracks. It can append several MC simulation runs to one Cosmic dataset file.
    * ```CTCoreModules``` -- this folder contains the ```Pattern Engine``` class that creates and handles super pixel mappings as well as the ```Template Bank```, which simulates an associative memory.
    * ```CTCosPatTrain``` -- scripts to build a database from a cosmic dataset.
    * ```CTBgAna``` -- scripts to evaluate the false-positive trigger rate for non-cosmic hit frames. These scripts use a Mu3e simulation MC file as input and compute the full combinatorics in each frame, create Templates from it and test, if they occur in a given template bank.
    * ```CTPlottingScripts``` -- several scripts to create all different curves, such as ROC, SPC-T-Count fit, etc.

* Additional Folders
    * ```data``` a data directory is created. Here the template database and CosmicSIDs files are stored after they are created. Also the Slimmed Files are stored here and the MC simulation files can be provided in this directory as well.
      * Data provided by Mu3e Simulation:
        * ```BackgroundData``` Mu3e Simulation .root files containing beam simulations (![Bkg-MCData][Module Overview])
        * ```SimulationData``` Mu3e Simulation TriRec .root files containing reconstructed cosmic muon tracks (![Cosmic-MCData][Module Overview])
      * Data produced by scripts of Cosmic Trigger:
        * ```SlimmedData``` Slimmed versions of the Cosmic Simulation Data (![CosmicData][Module Overview])
        * ```TemplateData``` The Cosmic Template Database files (![CosmicTDB][Module Overview])
        * ```CosmicSIDtrackData``` Cosmic Track SIDs used for external signal efficiency benchmarking
    * ```util``` contains several functions that are used within several modules and scripts.
    * ```karimaki``` an implementation of a three dimensional helix fit, used in ```CTPreSlimData```.
    * ```root``` classes for ROOT file reading and writing.
    * ```Mu3eRecAcc``` Mu3e Reconstruction Accuracy produces analysis and control plots for a Mu3e simulation TriRec output in cosmic mode.
    * ```output``` this folder is created to store all PDF and .root files created by the Cosmic Trigger (and are not stored in /data such as the template databases)
      * ```0_RecAcc``` Reconstruction Accuracy output
      * ```1_SlimSegs``` Slim Segs Script Control Plots
      * ```2_DBTraining``` Training Control Plots, Training overview plots, training control data
      * ```3_DBEvaluation``` Bkg Evaluation output data, Bkg Eval Plots 
      * ```4_CosmicSIDs``` 
      * ```5_TemplateFilterEff``` Histograms showing the cosmic efficiency with filters on
      * ```6_Tests``` Plots of test scripts
    

### How to set it up?

The project can be entirely build with cmake. As requirements, one needs to install

* root package (v6.20/04) 
* f2c libraries (fortran to c) *needed for karimaki helix fit*.

The f2c libraries are horrible to install, but there is a nice bash script available online, that does the job.
However, if not explicitly needed, the helix fit can just be excluded from the code and the project can then 
be compiled without.

To compile, follow these commands in the terminal:

```
mkdir build
cd build/
cmake ..
make
```

Once it is build, the executable should be executed from the projects top-level directory, e.g.

```
build/CTCosPatTrain/BuildDBSingleConfig 13 200 1 0.5 10
```

Yes, this is not very convenient, but I would just recommend to use an IDE, where you can choose your woking directory yourself and just hit the run button.

Also add the following directories to your top level:

```
mkdir data
mkdir output
```


### How to run it?

#### Data generation with the Mu3e Simulation

The Cosmic Trigger is a stand-alone program, that uses generated data from the Mu3e simulation package. 
This package can be installed from the bitbucket repository, if you need help, talk to Nik Berger.

Two different kinds of files are produced with the Mu3e sim, that are relevant for this study. 
First, it is necessary to simulate cosmic muon tracks. The simulation has a cosmic mode. To simulate the data just set 
the settings in the *mu3e/run/digi.json* file in the mu3e repo to the values shown in *Mu3eSimConfigFiles/digi_cosmic.json*.

Additionally, for the background evalutation, background (beam) data is needed. This can be obtained by using the normal 
mode of the mu3e simulation. The settings are provided in *Mu3eSimConfigFiles/digi_bg.json*. 

Once these steps are done and you can produce datat with the Mu3e Simulation Package, the files must be provided to the
cosmic trigger, and the Cosmic Trigger can be used.

#### Analysis with the Cosmic Trigger

Usually the Cosmic Trigger Code is used in several phases, which are:

1. Preparation of the data with ``PreSlimCosData``. This will result in a dataset file stored in data/
2. Database Training with the Executables in CTCosPatTrain, namely e.g. ```BuildDBMultiConfigs```. 
   The Parameters for super pixel mapping and the simulation in general can be specified in the file 
   *CTCoreModules/Configuration.h*. It centralises all parameter choices. For a future use of the Cosmic Trigger,
   this parameter selection should be moved to external config files to avoid recompiling.
3. Database Evaluation with Background data. The code is provided in the folder *CTBkgAna/* and Executables are e.g. ```EvalBkgDBMulti```
4. Finally, once the data is simulated and stored in *data/* and *output/*, the Plotting scripts can be used. Those can 
   found in *CTPlottingScripts/*. 
   
More information on how to run the Cosmic Trigger and how to choose the parameters is contained in *CTCoreModules/Configuration.h*
This directory also contains the Code for the ```TemplateBank``` and ```PatternEngine```.

### Get the data used in the first study

The data is provided on Tachyon at
```tachyon.physi.uni-heidelberg:/data/user_data/kneureither/mu3e/run/data```. This folder contains the raw data from 
the simulation and the corresponding *.json* config files. 

To reproduce the cosmic *dataset 13* that was used for most of the analysis within the thesis it is necessary to 
combine the *[...]trirec_cosmic.root* files of the following runs into one dataset by using ```PreSlimCosData```:

```
14,16,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37, 1000 to 1036, 1038 to 1065, 1100 to 1108 and 1110 to 1120.
```

Some of the run scripts are also included in this repository, in the folder *Mu3eSimConfigFiles/*.


### Who can I talk to?
If you need help, feel free to contact me via email neureither@physi.uni-heidelberg.de or dev@kneureither.de.

### One more thing
I do not really know how to use cmake, so if you do, feel free to clean up that project, unify commonly used libraries,
such as ```patterns``` or ```rootTreeFiles``` into the top level CMakeLists.txt or however this is done properly ;-)