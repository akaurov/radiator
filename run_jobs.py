import glob
import os

for i in glob.glob('jobs/*'):
    os.system('qsub '+i)
