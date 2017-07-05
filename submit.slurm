#! /bin/bash
#SBATCH -n 3
#SBATCH -p debug

export OMPI_MCA_mpi_warn_on_fork=0

module load r/3.3.3 mpi/openmpi/1.10.1-gcc gsl/2.4
Rscript simulths.R parallel