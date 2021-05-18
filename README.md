# This is a repo holding my personal work for the DESY summer camp.

[ More details concerning the algorithm developed or the results obtained can be found here https://www.desy.de/f/students/2019/reports/Edoardo.Centamori.pdf ]

During my stay in Hamburg I have helped a research team at DESY to implement a Ring Polymer Molecular Dynamics symulation. I've worked there for roughly 3 months at the end of my bachelor and implemented in their code a few algorithms that I have personally created in order to impose certain constraint to the simulation and also optimized their code in a few ways. In this repo is contained not their code (which is private) but what I used to test out my algorithm, a toy model usefull for performing generic Path Integral Molecular Dynamics symulation (In practice is a complete molecular dynamics symulation but without a tuning of the parameters to represent a realistic physical model).

There are three 'front end' files (not really) which are used respectively to setup the system, run it and visualize the results:
1. setting.py : Here all the relevant properties of the system can be setted
2. main.py : This file is used to run the symulation and store the output data in a byte file
3. Visualizer.py : This script is used to visualize the results contained in a byte file, it allows to play a video of the particles moving and was used mainly for debugging purposes.

There are instead several files which are responsible for the various computation and algorithm implemented.


1. verletM.py : Responsible for the main computation, all the function that evolves the system during the simulations are contained here 
2. normal.py : Contains several utilities that are useful in one of the verletM.py algorithms concerning the propagation in terms of normal modes.
3. rattle.py : Responsible for the insertion in the algorithm of constraints that can be imposed during the evolution of the system.
4. thermostat.py : Sumperimpose to the system a random process that allows for thermalization. The constraints in rattle have been devices in such a way that they're not affected by the termalization process at all.
5. estimatorsM.py : All the relevant observables of the system are defined here.
6. There are other script and files which are less relevant and were used mainly for debugging and testing purpuses.


