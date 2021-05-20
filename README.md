# Ring Polymer Molecular Dynamics simulation, DESY summer camp.

A more detailed explanation concerning the algorithm developed or the results obtained can be found here https://www.desy.de/f/students/2019/reports/Edoardo.Centamori.pdf.

At the end of my bachelor degree I've spent roughly 3 months at the German synchrotron in DESY (Hamburg), I have helped a research team implementing a Ring Polymer Molecular Dynamics symulation (RPMD).

# What's RPMD?

In a few words, RPMD is an algorithm that allows you to calculate numerically the properties of out of equilibrium quantum system by performing symulations of equivalent systems that are classical in nature. In practice a single quantum particle can be simulated through the classical evolution of a necklace of beads connected by springs. 

# What's this project?

During my permanence in DESY I designed and developed an algorithm that allows you to impose certain constraints on the evolution of these classical systems. These constraints are usefull because allows you to symulate the kind of restiction existing in real molecules.
While I was working on the code developed by the research team both impoving his efficiency and adding my personal contribution, I've also wrote my personal molecular dynamic symulation from scratch.
This repository contains my simulation while I can't obviously share the algorithm of the research team. The algorithm here should be though as a toy model that was very useful to test and debug what I was developing without having to deal with the complexity of a full RPMD algorithm. 

# What's going on in the algortigm?

1. PIMD (Path integral molecular dynamic): As better explained in the report linked at the beginning, in equilibrium quantum system there is a correspondence with classical system which is what we're going to simulate.
2. Stormer-Verlet algorithm: It's a famous second-order integrator to solve differential equation, I've used it to simulate the evolution of the classical system.
3. Langevin thermostat and stochastic differential equations: The sistem we're trying to simulate is at thermal equilibrium and thus you should make it interact with an infinitely large thermostat to behave correctly, this is unduable computationally. Those who know stochastic calculus will understand that this kind of behaviour can be obtained by superimposing to the classical evolution a small noise with the right properties, this is what is called a langevin thermostat.
4. Rattle algorithm: It's a well known algorithm that allows to set certain constraints to a classical system (e.g. evolve an object fixing his ortientation, or fix the distance between two objects). The algorithm is not thought to work symultaneously with the thermostat and the small noise can cause some troubles, but I have adapted it to work just fine in my algorithm, some more details can be obtained in the project report linked at the beginning.

# What are the various files?

There are three files which are used respectively to setup the system, run it and visualize the results:

1. setting.py : Here all the relevant properties of the system can be setted (Real simulation can take some times, I used to run this in a cluster at DESY, I'll leave parameter that allows for a small simulation).
2. main.py : This file is used to run the symulation and store the output data in a byte file.
3. Visualizer.py : This script is used to visualize the results contained in a byte file, it allows to play a video of the particles moving and was used mainly for debugging purposes (I will leave in the repo also a few byte files obtained from previous run so that they can be visuliazed in Visualizer.py). 

There are also several files which are responsible for the various computation and algorithm implemented.

1. verletM.py : Responsible for the main computation, the various implementations of the Stormer-Verlet algorithm are contained in this file.
2. normal.py : Contains several utilities that are useful in one of the verletM.py algorithms concerning the propagation in terms of normal modes.
3. rattle.py : Responsible for the insertion in the algorithm of constraints that can be imposed during the evolution of the system.
4. thermostat.py : Sumperimpose to the system a random noise with the right properties that allows for canonical thermalization. The constraints in rattle have been devices in such a way that they're not affected by the termalization process at all.
5. estimatorsM.py : All the relevant observables of the system are defined here.
6. There are other script and files which are less relevant and were used mainly for debugging and testing purpuses.

# How to run the algorithm 

WIP
