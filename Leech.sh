#! /bin/bash
#
#SBATCH --job-name=M5LeechPlot
#SBATCH --output=M5log.txt
#SBATCH --error=M5err.txt
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --chdir=/home/ssavitz/QC
#SBATCH --open-mode=truncate

echo NEW RUN
date
echo NEW RUN

echo NEW RUN >&2
date >&2
echo NEW RUN >&2

./PlotMedium5

echo DONE
date
echo DONE

echo DONE >&2
date >&2
echo DONE >&2
