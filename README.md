# The Right Way to code simulation studies in Stata
This repository accompanies my talk at the 25th London Stata Conference, 5–6 Sep 2019.
It contains three .do files which implement the same simulation study using different approaches in Stata.
## 1TheSimulateWay.ado
This file uses the `simulate` command to run the simulation study. The code is clear but has various shortcomings (as noted in my presentation), primarily that it produces messy data files, which may lead to errors, and that it cannot save random number generator states, which is often desirable.
## 2ThePostWay.ado
This file uses the `post` suite of commands to run the simulation study. The file is less clean in that code for simulating and analysing data are entangled with code for storing inputs and outputs. The approach produces a clean data file, PostEstimates.dta, and a file containing n_{sim}+1 random number states.
## 3TheRightWay.ado
The Right Way is a mashup of the previous two Ways: as with `simulate`, an rclass program is defined which returns outputs as scalars. The program is repeatedly called by a `forvalues` loop. Inputs and outputs are stored using the `post` suite. This has the advantages of clean code – as with The `simulate` Way – and clean data and states – as with The `post` Way.
