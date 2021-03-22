## CODE TO MOVE THE NEW FAA FILES INTO THE FOLDER
import os
import subprocess as sp

# read in genome ids with empty faa files from command line search
patric_table = pd.read_csv('combined_genomes.csv') 
genomeids = list(set(list(patric_table.genomeid)))
# iterate over genome ids
for gid in genomeids:
    cmd = 'cp ./prokkas/{n}/{n}.faa ./faafiles'.format(n = gid)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    print(out_cmd.decode())
