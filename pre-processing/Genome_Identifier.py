import os
import pandas as pd
import subprocess as sp
from pathlib import Path

## NOTE: THIS CODE CAN ONLY BE RUN ON PATRIC COMMAND LINE INTERFACE
## https://docs.patricbrc.org/cli_tutorial/cli_installation.html

# List of taxon ids: E. coli (562), TB (77643), Salmonella (28901), and staph (1280)
ids = [562, 77643, 28901, 1280]
for taxon in ids: 
    # I download the metadata from PATRIC
    try:
        os.makedirs("./metadata/sources/patric/{}".format(str(taxon)), exist_ok=True)
        os.makedirs("./resistance_data/sources/patric/{}".format(str(taxon)), exist_ok=True)
    except:
        raise(" * Creation of the directory ./metadata/sources/patric/ failed")
    # EDIT POINT! taxon_id
    cmd='p3-all-genomes --eq "taxon_id,{}" -a biosample_accession,geographic_location,isolation_country,collection_date,collection_year'.format(str(taxon))

    print(" * Downloading the IDs from PATRIC")
    print("    * executing: {}".format(cmd))
    p=sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd,err=p.communicate()
    print("... Finished!")
    with open("./metadata/sources/patric/{}/patric_ids_geo_isoltime_data.tsv".format(str(taxon)),"w") as outf:
        outf.write(out_cmd.decode())

    # I load the table of the IDs with pandas
    tab=pd.read_csv("./metadata/sources/patric/{}/patric_ids_geo_isoltime_data.tsv".format(str(taxon)), sep="\t", dtype={
    "genome.genome_id":object,
    "genome.biosample_accession":object,
    "genome.geographic_location":object,
    "genome.isolation_country":object,
    "genome.collection_date":object,
    "genome.collection_year":object
    })
    # I select only the strains that have a biosample
    tab_subset=tab.loc[tab['genome.biosample_accession'].notnull(),]
    list_ids=list(tab["genome.genome_id"])
    # split this list in chunks of 250 strains and request the data
    n=250
    chunks = [list_ids[i * n:(i + 1) * n] for i in range((len(list_ids) + n - 1) // n )]  
    # I request the data for each chunk
    counter=0
    for chunk in chunks:
        counter=counter+1
        if Path("./resistance_data/sources/patric/{}/input{}.txt".format(str(taxon), str(counter))).is_file():
            continue
        with open("./resistance_data/sources/patric/{}/input{}.txt".format(str(taxon), str(counter)), "w") as outf:
            outf.write("genome.genome_id\n")
            outf.write("\n".join(chunk))
        # I execute the query
        print(" * Executing query for chunk no. {} of {}".format(str(counter),str(len(chunks))))
        cmd="p3-get-genome-drugs -i " + "./resistance_data/sources/patric/{}/input{}.txt".format(str(taxon), str(counter)) + ' --eq "evidence,Laboratory Method" -a antibiotic,biosample_accession,date_inserted,date_modified,evidence,genome_id,laboratory_typing_method,laboratory_typing_method_version,laboratory_typing_platform,measurement,measurement_sign,measurement_unit,measurement_value,owner,pmid,public,resistant_phenotype,source,testing_standard,testing_standard_year,text'
        print("    * executing: {}".format(cmd))
        p=sp.Popen(cmd, shell=True, stdout=sp.PIPE)
        out_cmd,err=p.communicate()
        with open("./resistance_data/sources/patric/{}/amr_data{}.tsv".format(str(taxon), str(counter)),"w") as outf:
            outf.write(out_cmd.decode())
    # I merge the data 
    print("    * Merging the data: {}".format(cmd))
    flag_header=0
    with open("./resistance_data/sources/patric/{}/amr_data_{}.tsv".format(taxon, taxon),"w") as outf:
        for counter in range(1,len(chunks)+1): 
            with open("./resistance_data/sources/patric/{}/amr_data{}.tsv".format(str(taxon), str(counter))) as inp:
                h=inp.readline()
                if flag_header==0:
                    outf.write(h)
                    flag_header=1
                content=inp.read()
                outf.write(content)
