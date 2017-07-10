# README #

Three states statistical modeling of animal movements (Moving-Resting-Hunting process).

### Installation ###

You need GNU Scientific Library (GSL) installed.
On ubuntu you can install GSL using `sudo aptitude install libgsl0-dev`

Then install the required R packages:

``` R
install.packages(c("argparser", "RcppGSL", "nloptr"))
```

To run the parallel simulation using the `--parallel` option,
ou also need an MPI library installed.
On ubuntu you can install openmpi using `sudo aptitude install libopenmpi-dev openmpi-bin`.
Then install the R packages:

``` R
install.packages(c("snow", "Rmpi"))
```

### Usage ###

The `simulths.R` script generates simulated data and fits the model.
You can run it with:

``` bash
# Explanation of the input arguments
Rscript simulths.R --help
```

``` bash
# Run simulation on HPC cluster using MPI and the SLURM scheduler
Rscript simulths.R --parallel
```

### Files Description ###

* `coga2dim.cpp` : Calculate coga2dim.

* `fitMovResHun.R` : estimation for thmam.

* `rMovResHun.R` : Generate random data for thmam.

* `likelihood.cpp` : Old version of likelihood.

* `simulths.R` : Simulation work for estimation.

* `thmam_hyper.cpp` : Current version of likelihood.

* `submit.slurm` : configure file for HPC.

* `thama_likelihood.cpp` : likelihood file contains serial and parallel versions.
