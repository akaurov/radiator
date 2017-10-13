import numpy as np 

z_list = np.round(np.logspace(0, 3, 16))
d_list = [-0.5, 0.0, 0.5, 1.0, 10.0, 100.0, 1000.0, 10000.0]
x_ion = [2e-4, 1e-3, 1e-2, 1e-1, 0.5, 0.9, 0.99]

#python primary_electron_de.py 19.0 100.0 0.1 VFKY1996 RBEQ RBEQ


for i in range(len(x_ion)):
    for j in range(len(d_list)):
        f1=open('jobs/'+str(x_ion[i])+'_'+str(d_list[j])+'.sh', 'w')
        temp = '''#!/bin/bash
#SBATCH --job-name=radiator'''+str(x_ion[i])+'_'+str(d_list[j])+'''
#SBATCH --output=/home/kaurov/scratch-midway/radiator/jobs/radiator'''+str(x_ion[i])+'_'+str(d_list[j])+'''.out
#SBATCH --error=/home/kaurov/scratch-midway/radiator/jobs/radiator'''+str(x_ion[i])+'_'+str(d_list[j])+'''.err
#SBATCH --time=23:00:00
#SBATCH --partition=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

cd /home/kaurov/scratch-midway/radiator


'''
        for k in range(len(z_list)):
            temp+='~/anaconda/bin/python primary_electron_de_more.py '+str(z_list[k])+' '+str(d_list[j])+' '+str(x_ion[i])+' VFKY1996 RBEQ RBEQ & \n'
        temp+='wait'
        f1.write(temp)
