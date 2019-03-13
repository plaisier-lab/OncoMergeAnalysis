# For running the OncoMerge program generally
from subprocess import *
import pandas as pd
import os
import csv
from multiprocessing import Pool, cpu_count
import time
##############################################
### Enter Parameters and File Paths Here
##############################################

# Unique Name for the current run:
Run_Name = 'PanCancer_Final'

# Minimun Mutation Frequency (recommend 0.05):
# if you would like to see results from different mutation frequency values, you may enter multiple values into the list... recommend entering one.
mutation_frequency =  [0.05]

# Paths to Gistic Files and Somatic Mutations File:
#Each must be a string.
#if your path is independent of the cancer's TCGA abbreviation, enter your path as the string
#if your path is not independent, you may enter an * in place of the TCGA abbreviation and OncoMerge will automatically insert it for each cancer
#Lastly if your path includes a TM or TP, you may enter an * in place of the TCGA abbreviation and the TP or TM, first * will be replaced with the TCGA abbreviation and the second * will be replaced with the TM or TP. OncoMerge will automatically fill the information for each cancer. TM will be used for SKCM and TP for all other tumor types.

#Deletions File:
del_path = 'GISTIC_99/gdac.broadinstitute.org_*-*.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt'
    #Amplifications File:
amp_path = 'GISTIC_99/gdac.broadinstitute.org_*-*.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt'
    # Gene Data File:
gene_data_path = 'GISTIC_99/gdac.broadinstitute.org_*-*.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt'
    #all_Thresholded File:
thresh_path = 'GISTIC_99/gdac.broadinstitute.org_*-*.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt'

# Somatic Mutations File:
#Somatic Mutations will only insert a cancer abbreviation and will only do so when a list of two strings is given.
somatic_mutations_path = 'mutations_mc3/*_somMutMC3.csv'

# Path to Output:
#Please give as a string, not a list
output_path = 'output/OncoMerge'
#List of cancers (can be only one cancer in list)
cancers = sorted(['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC'])
##############################################

# Run OncoMerge
def RDProc(a):
    start = time.time()
    print(' '.join(a))
    logFile = open('log/'+a[2]+'_mf_'+a[6]+'_.log','w')
    RDPproc = Popen(' '.join(a), stdout=logFile, stderr=STDOUT, shell = True)
    output = RDPproc.communicate()
    end = time.time()
    print(end-start)



def main():
    runMe = []
    for mf in mutation_frequency:
        for cancer in cancers:
            if not os.path.exists(output_path+'/'+Run_Name+'_'+cancer+'_mf_'+str(mf)+'_OncoMerged.csv'):
                runMe.append(['python OncoMerge.py', '--tumorType', cancer,'-r',Run_Name,' --mf',str(mf),'-d',del_path,'-a',amp_path,'-g',gene_data_path,'-t',thresh_path,'-s',somatic_mutations_path,'-o',output_path])

    cpus = 2
    # cpus = cpu_count()
    print('There are %d CPUs available.' % cpus)
    pool = Pool(processes=cpus)
    pool.map(RDProc, runMe)
    pool.close()
    pool.join()

#for a in runMe:
#    RDProc(a)

if __name__ == "__main__":
    main()
