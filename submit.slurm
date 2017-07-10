#! /bin/bash
#SBATCH -n 10
#SBATCH -p general
#SBATCH --mail-user=chaoran.hu@uconn.edu
#SBATCH --time=03-02:03:04

export OMPI_MCA_mpi_warn_on_fork=0

module load r/3.3.3 mpi/openmpi/1.10.1-gcc gsl/2.4
Rscript simulths.R --parallel
