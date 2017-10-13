#!/bin/bash
#SBATCH --job-name=radiator0.0002_1.0
#SBATCH --output=/home/kaurov/scratch-midway/radiator/jobs/radiator0.0002_1.0.out
#SBATCH --error=/home/kaurov/scratch-midway/radiator/jobs/radiator0.0002_1.0.err
#SBATCH --time=23:00:00
#SBATCH --partition=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

cd /home/kaurov/scratch-midway/radiator


~/anaconda/bin/python primary_electron_de_more.py 1.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 2.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 3.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 4.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 6.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 10.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 16.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 25.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 40.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 63.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 100.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 158.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 251.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 398.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 631.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
~/anaconda/bin/python primary_electron_de_more.py 1000.0 1.0 0.0002 VFKY1996 RBEQ RBEQ & 
wait