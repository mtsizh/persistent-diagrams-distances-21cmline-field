# persistent-diagrams-distances-21cmline-field
Code to reproduce results of the corresponding paper, studying distances between persistent diagrams of 21 cm line brightness field modeled with 21cmfast with varying cosmological parameters.

You will need the following packages installed: 
- 21cmfast
- gudhi
- scikit-bio

* lightcones_and_snapshots.py will produce a set of cosmological simulations via 21cmfast package. It will save the simulations in separate .h5 files.
* diag_dist.py will compute persistent diagrams of brightness fields taken from simulations produced with lightcones_and_snapshots.py, then it will compute Wasserstein distances for pairs of simulations with different cosmological parameters. After this, it will perform a PERMANOVA test with distance matrices to test the hypothesis that persistent diagrams of different cosmologies can be distinguished with Wasserstein distance. 
