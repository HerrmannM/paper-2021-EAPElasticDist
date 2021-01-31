# Early Abandoned & Pruned for Elastic Distances

This code is designed to be used with the UCR archive using the ts format.
See http://timeseriesclassification.com/dataset.php
It implements DTW, CDTW, WDTW, ERP, MSM and TWE in four versions:
  *  base: a classic double buffered implementation
  *  base_ea: same as base + the classical way of early abandoning
  *  eap: double buffered implementation with pruning, using pruning to early abandon
  *  pru: pruning only. Same code as eap, but compute its own threshold based on a possible alignment.
     Hence, it cannot early abandon.
     
The code also implements SQED and LCSS, albeit without pruning.
Finally, some `_la` version of the code exists.
This is an attempt at speeding up further the computation by tightening the provided threshold using the
last alignement of the series. However, the results are hardware dependent.
In some cases, the cache faults are more expensive than the gain. Hence this option is not considered in our
experimentations.

## Building

All the C++ source code is available in the `code/src` directory.
Use cmake to build the project.
For example, in the `code` folder:
```
mkdir cmake-build-release
cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
./experiments
```

The last command shows a help message.

## Experimentation

The `code/launcher` contains a python script able to launch several commands in parallel,
along with a CSV containing parameters found by EE for the UCR Archive (in its 85 datasets version).
It also contains two folders.
One for an experimentation comparing run times of NN1 classifiers with elastic distances in several configurations.
One for an experimentation studying the impact of lower bounding.
Each folder contains:
  *  a python script generating a file containing commands, based on the EE parameters. This file can be used with the first script.
  *  a python script analysing the results and producing some graphs
  *  a results/json.zip file containing all our results
