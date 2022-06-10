#! /bin/bash
#
#SBATCH --job-name=MLeechPlot
#SBATCH --output=logM.txt
#SBATCH --error=errM.txt
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --chdir=/home/ssavitz/QC
#SBATCH --open-mode=truncate

##SBATCH --signal=USR1@120 


#t=$((24*60*60 - 5))


echo NEW RUN
date
echo NEW RUN

echo NEW RUN >&2
date >&2
echo NEW RUN >&2

./PlotMedium

echo DONE
date
echo DONE

echo DONE >&2
date >&2
echo DONE >&2
