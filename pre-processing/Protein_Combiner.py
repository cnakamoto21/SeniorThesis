## CODE TO DOWNLOAD ALL OF THE FNA FILES   
import pandas as pd

# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
# open files to write all the proteins to
rpofiles = [open('rpo{}.faa'.format(i), 'w') for i in ['a', 'b', 'c', 'o']]
gyrfiles = [open('gyr{}.faa'.format(i), 'w') for i in ['a', 'b']]
errors = open('comberrors.txt', 'w')
# only probe each isolate once - purge redundant genome ids
gids = list(set(list(patric_table.genomeid)))
for gi in gids:
    # initialize global checking variables
    rpoidx = [[] for i in range(4)]
    gyridx = [[] for i in range(2)]
    # read in faa file contents
    with open ("./faafiles/{}.faa".format(gi), "r") as myfile:
        elements = myfile.read().split(">")
    # iterate over rows to find proteins, store indices
    for ele in elements:
        if 'DNA-directed RNA polymerase' in ele:
            if 'alpha' in ele:
                rpoidx[0].append(elements.index(ele))
            elif 'beta' in ele:
                rpoidx[1].append(elements.index(ele))
            elif 'omega' in ele:
                rpoidx[3].append(elements.index(ele))
            elif 'beta\'' in ele:
                rpoidx[2].append(elements.index(ele))
        if 'DNA gyrase subunit B' in ele:
            if 'alpha' in ele:
                gyridx[0].append(elements.index(ele))
            elif 'beta' in ele:
                gyridx[1].append(elements.index(ele))
        if len(rpoidx[0]) > 0 and len(rpoidx[1]) > 0 and len(rpoidx[2]) > 0 and len(rpoidx[3]) > 0 and len(gyridx[0]) > 0 and len(gyridx[1]) > 0: 
            break
    if [] in rpoidx or [] in gyridx: 
        errors.write("Error: rpo is {}, gyrb is {} for genome {}".format(rpoidx, len(gyrbidx), gi))
    # write corresponding aa sequences into files - only take first for each isolate
    for i in range(4):
        for rpo in rpoidx[i]: 
            rpofiles[i].write('>{}\n'.format(gi))
            rpofiles[i].writelines(elements[rpo].split('\n')[1:])
            rpofiles[i].write('\n')
            break
    for i in range(2):
        for gyr in gyridx[i]: 
            gyrfiles[i].write('>{}\n'.format(gi))
            gyrfiles[i].writelines(elements[gyr].split('\n')[1:])
            gyrfiles[i].write('\n')
            break
    # status check
    if idx % 100 == 0:
        print("completed", idx)
# close protein files
for file in rpofiles: 
    file.close()
for file in gyrfiles: 
    file.close()
errors.close()
