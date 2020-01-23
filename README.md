# _k_-plex
This repository contains codes used in the computational study reported in the article for *"Approaches for finding cohesive subgroups in large‐scale social networks via maximum k‐plex detection"* that has been published in the _Nextworks_. If you wish to use or cite this code, please cite the paper:
```
@article{BBZMkplexGRASP2017,
Author = {Zhuqi Miao and Balabhaskar Balasundaram},
Journal = {Networks},
Month = {July},
Number = {4},
Pages = {388--407},
Title = {Approaches for finding cohesive subgroups in large-scale social networks via maximum $k$-plex detection},
Volume = {69},
Year = {2017}}
```

User instructions: Instances and packages required to run codes-- not provided here, other licenses required; compile and run instructions-- input parameters, etc.

* kplex.cpp: The code for solving the maximum _k_-plex problem.
* Makefile: Make file to compile kplex.cpp in a linux system.

# Compiling the code
1. Configure the Linux environment with Gurobi as described in the paper.
2. Copy kplex.cpp and Makefile into a same folder
3. In terminal, go to the folder. Then type "make" and hit enter to compile kplex.cpp

# Execution
1. Create an _input_ folder to place input instances, and create an _output_ folder for outputs
2. The code create an executale program named "prog". Use "qsub" command or double click to run it depending on the Linux environment
