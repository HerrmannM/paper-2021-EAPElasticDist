# Early Abandoned & Pruned for Elastic Distances

This code supports a paper currently submitted, and is not intended to be use as a library (although you can!).
If that's what you are looking for, please have a look at the [tempo library](https://github.com/MonashTS/tempo).

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

## Experimentations

The `nnclassification/launcher` contains a python script able to launch several commands in parallel,
along with a CSV containing parameters found by EE for the UCR Archive (in its 85 datasets version).
It also contains two folders.
One for an experimentation comparing run times of NN1 classifiers with elastic distances in several configurations.
One for an experimentation studying the impact of lower bounding.
Each folder contains:
  *  a python script generating a file containing commands, based on the EE parameters. This file can be used with the first script.
  *  a python script analysing the results and producing some graphs
  *  a results/json.zip file containing all our results
 
The third experiment presented in the paper on sub-sequene search is contained in the `subsequence` folder.
  
  
## Benchmarks
To produce the bechmarks output, unzip the json.zip files inside their results folder.
All the json files it contains should now reside under `code/launcher/exp_something/results/json`.
Then, launch the corresponding `analys_results.py` script to generate the figures.
