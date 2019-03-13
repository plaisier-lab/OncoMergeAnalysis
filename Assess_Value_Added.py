#assess the value added from OncoMerge.
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from matplotlib import pyplot as plt
from subprocess import *
import os
import pandas as pd
from scipy.stats import hypergeom as hyg




#1. Reads data from OncoMerged_file into glofmf as the dictionary of dictionaries {Cancer1:gene1:MutType1:Mutation Frequency}

oops = []
muts = {}
glofmf = {}
glofdf = pd.DataFrame([],columns = ['Delta','Cancers', 'Group'])
reinfdf = pd.DataFrame([],columns = ['Delta','Cancers'])
recovdf = pd.DataFrame([],columns = ['Delta','Cancers'])
cancers = sorted(['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC'])
glofs = {}
OncoMerged_File = 'output/GLoF_Violins/MutFile/'+cancer+'_deepFilteredMutFile.csv'
for cancer in cancers:
    print('Pulling in info for ' + cancer)
    #1a. parses OncoMerged_File
    with open(OncoMerged_File, 'r') as inFile:
        reader = [i.split(',') for i in inFile.read().split('\n')[:-1]]
        pats = reader[0][1:]
        glofmf[cancer] = {}
        for line in reader[1:]:
            if not line[0].split('_')[0] in glofmf[cancer].keys():
                glofmf[cancer][line[0].split('_')[0]] = {}
            glofmf[cancer][line[0].split('_')[0]][line[0].split('_')[1]] = sum([int(i) for i in line[1:]])/len(line[1:])
    #1b Create dataframe for violinplot of individual cancers
    tmpgenes = [i for i in glofmf[cancer] if 'GoF' in glofmf[cancer][i] or 'LoF' in glofmf[cancer][i]]
    for gene in tmpgenes:
        # Calculate delta
        if 'GoF' in glofmf[cancer][gene]:
            list_mf = [glofmf[cancer][gene][i] for i in glofmf[cancer][gene].keys() if i!='GoF']
            delta = glofmf[cancer][gene]['GoF'] - max([glofmf[cancer][gene][i] for i in glofmf[cancer][gene].keys() if i!='GoF'])
        elif 'LoF' in glofmf[cancer][gene]:
            list_mf = [glofmf[cancer][gene][i] for i in glofmf[cancer][gene].keys() if i!='LoF']
            delta = glofmf[cancer][gene]['LoF'] - max([glofmf[cancer][gene][i] for i in glofmf[cancer][gene].keys() if i!='LoF'])
        else:
            oops.extend(gene)
        # classify as reinf and recov
        if sum([i>=0.05 for i in list_mf])>0:
            reinfdf.loc[gene+'_'+cancer] = [delta,cancer]
            glofdf.loc[gene+'_'+cancer] = [delta,cancer,'reinf']
        if sum([i>=0.05 for i in list_mf])==0:
            recovdf.loc[gene+'_'+cancer] = [delta,cancer]
            glofdf.loc[gene+'_'+cancer] = [delta,cancer,'recov']



# Create OncoPlotTable
stats = ['Symbol', 'GoF or LoF','Locus', 'GLoF Mut Rate','CNA Mut Rate','PAM Mut Rate','Delta','gistic amp res q value','gistic del res q value','mutsig2cv q value','Group','Cancer']
allstatdf = pd.DataFrame([],columns=stats)
MutSig2CV_Genes = {}
OM_Genes = {}
Gistic_Genes = {}
for cancer in cancers:
    #Load in cancer
    if cancer == 'SKCM':
        path = 'GISTIC_99/gdac.broadinstitute.org_'+cancer+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0'
    else:
        path = 'GISTIC_99/gdac.broadinstitute.org_'+cancer+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0'
    d0 = pd.read_csv(path+'/all_thresholded.by_genes.txt',sep='\t')
    d1 = d0.set_index('Locus ID')
    d1.columns = list(d1.columns[0:2])+[i[:12] for i in d1.columns[2:]]
    n2 = pd.Series(d1['Gene Symbol'], index =d1.index.values)
    somMuts = pd.read_csv('mutations_mc3/'+cancer+'_somMutMC3.csv',index_col=0,header=0).index
    statdf = pd.DataFrame([],columns=stats)
    mutsig = pd.read_csv('sig2cv/'+cancer+'_sig2cv.csv', index_col = 1, header = 0)
    # Load up CNA Loci
    loci = pd.read_csv('output/cna_loci/'+cancer+'_0.05_cna_loci.csv',index_col = 0)
    for i in loci.index:
        loci.loc[i, 'Genes'] = loci.loc[i, 'Genes'].split(' ')
    #Calculate enrichment of recovered genes.
    # insert hygs...
    #Pull in amps from gistic
    amp1 = pd.read_csv(path+'/amp_genes.conf_99.txt',index_col=0,sep='\t')
    ampdict = {}
    for col1 in amp1.columns:
        ampGenes = [i.lstrip('[').rstrip(']') for i in list(amp1[col1].dropna()[3:])]
        for gene in ampGenes:
            ampdict[gene] = float(amp1.loc['residual q value', col1])
    #Pull in dels from gistic
    del1 = pd.read_csv(path+'/del_genes.conf_99.txt',index_col=0,sep='\t')
    deldict = {}
    for col1 in del1.columns:
        delGenes = [i.lstrip('[').rstrip(']') for i in list(del1[col1].dropna()[3:])]
        for gene in delGenes:
            deldict[gene] = float(del1.loc['residual q value', col1])
    # Extract Stats for OncoPlotTable
    for gene in glofdf.loc[glofdf['Cancers']==cancer].index:
        gene0, c0 = gene.split('_')
        locus = 'NA'
        if 'GoF' in glofmf[cancer][gene0]:
            oF = 'GoF'
            CNA = 'CNAamp'
        else:
            oF = 'LoF'
            CNA = 'CNAdel'
        Symbol = n2.loc[int(gene0)]
        if Symbol in mutsig.index:
            sig2cvq =  mutsig.loc[Symbol, 'q']
        else:
            sig2cvq = 'NA'
        if Symbol in ampdict.keys():
            ampq = ampdict[Symbol]
            for i in loci.index:
                if gene in loci.loc[i, 'Genes']:
                    locus = i
        else:
            ampq = 'NA'
        if Symbol in deldict.keys():
            delq = deldict[Symbol]
            for i in loci.index:
                if gene in loci.loc[i, 'Genes']:
                    locus = i
        else:
            delq = 'NA'
        # Create OncoPlot Table
                            #Symbol  GoF/Lof, Locus     GLoF Mut Rate,              'CNA Mut Rate',             'PAM Mut Rate',               'Delta',                amp,  del ,'mutsig q value',       Group,            'Cancer'
        statdf.loc[gene] = [Symbol,      oF,  locus, glofmf[cancer][gene0][oF], glofmf[cancer][gene0][CNA], glofmf[cancer][gene0]['PAM'], glofdf.loc[gene, 'Delta'], ampq, delq, sig2cvq         , glofdf.loc[gene, 'Group'], cancer]
    # Save some data for Enrichment Analyses
    OM_Genes[cancer] = set([n2.loc[int(key)] for key in glofmf[cancer].keys()])
    sigdels = [i for i in deldict.keys() if deldict[i]<=0.05]
    sigamps = [i for i in ampdict.keys() if ampdict[i]<=0.05]
    Gistic_Genes[cancer] = set(sigdels).union(sigamps)
    MutSig2CV_Genes[cancer] = set([i for i in mutsig.index if mutsig.loc[i, 'q']<=0.05])
    statdf.to_csv('output/OncoPlotTables/'+cancer+'_OncoPlotTable.csv')
    allstatdf = allstatdf.append(statdf)


allstatdf.to_csv('output/OncoPlotTables/All_OncoPlotTable.csv')



# Define the Order of each violin plot

means_all = {}
medians_all = {}
#means_reinf= {}
#medians_reinf = {}
#means_recov = {}
#medians_recov = {}
for cancer in cancers:
    #medians_reinf[cancer] = reinfdf.loc[[i for i in reinfdf.index if cancer in i]]['Delta'].median()
    #means_reinf[cancer] = reinfdf.loc[[i for i in reinfdf.index if cancer in i]]['Delta'].mean()
    #medians_recov[cancer] = recovdf.loc[[i for i in recovdf.index if cancer in i]]['Delta'].median()
    #means_recov[cancer] = recovdf.loc[[i for i in recovdf.index if cancer in i]]['Delta'].mean()
    medians_all[cancer] = glofdf.loc[[i for i in glofdf.index if cancer in i]]['Delta'].median()
    means_all[cancer] = glofdf.loc[[i for i in glofdf.index if cancer in i]]['Delta'].mean()

#order_mean_reinf = ['All']
#order_mean_reinf.extend(sorted(means_reinf,key = means_reinf.get))
#order_median_reinf = ['All']
#order_median_reinf.extend(sorted(medians_reinf,key = medians_reinf.get))
acdf_reinf = pd.concat([reinfdf['Delta'],pd.Series(['All' for i in reinfdf.index],name = 'Cancers',index = reinfdf.index)], axis = 1)
reinfAll = reinfdf.append(acdf_reinf)

#order_mean_recov = ['All']
#order_mean_recov.extend(sorted(means_recov,key = means_recov.get))
#order_median_recov = ['All']
#order_median_recov.extend(sorted(medians_recov,key = medians_recov.get))
acdf_recov = pd.concat([recovdf['Delta'],pd.Series(['All' for i in recovdf.index],name = 'Cancers',index = recovdf.index)], axis = 1)
recovAll = recovdf.append(acdf_recov)

order_mean_all = ['All']
order_mean_all.extend(sorted(means_all,key = means_all.get))
order_median_all = ['All']
order_median_all.extend(sorted(medians_all,key = medians_all.get))
acdf_all = pd.concat([glofdf['Delta'],pd.Series(['All' for i in glofdf.index],name = 'Cancers',index = glofdf.index)], axis = 1)
glofAll = glofdf.append(acdf_all)





scale = 'area'
#scale = 'count'

#Make violin plots for each cancer:
with PdfPages('plots/ValueAdded/ValueAddedViolins_meanSorted_area.pdf') as pdf:
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=glofAll, cut = 0, inner = None, order = order_mean_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=reinfAll, cut = 0, inner = None, order = order_mean_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=recovAll, cut = 0, inner = None, order = order_mean_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()

with PdfPages('plots/ValueAdded/ValueAddedViolins_medianSorted_area.pdf') as pdf:
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=glofAll, cut = 0, inner = None, order = order_median_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=reinfAll, cut = 0, inner = None, order = order_median_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()
    ax = sns.violinplot(x = 'Cancers', y = 'Delta', data=recovAll, cut = 0, inner = None, order = order_median_all, scale = scale)
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()




acdf_reinf['Group'] = ['reinf' for i in acdf_reinf.index]
acdf_recov['Group'] = ['recov' for i in acdf_recov.index]


acdf_groups = pd.concat([acdf_reinf,acdf_recov])




with PdfPages('plots/ValueAdded/Compare_GLoF_Groups.pdf') as pdf:
    ax = sns.violinplot(x = 'Group', y = 'Delta', data=acdf_groups, cut = 0, inner = None, scale = 'count')
    plt.ylim(0, .2)
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()

#Add in Violinplot comparing groups's MF as well.


##############################################
### Plot OncoPlots for top 10 in each group
##############################################

print('Setting Up...')

def PlotOncoPlot(name, path,cancer):
    runMe = ['require(maftools)'] #Needed for OncoPlot() and read.maf()
    runMe += ['require(ggplot2)'] #Needed for ylab()
    # Load in the tmpmaf.tsv as maf
    runMe += ['maf <- read.table(\'tmpmaf.tsv\', header =1,sep = \'\\t\', quote = "\\"")']
    runMe += ['del1 <-paste(\''+path+'\', \'del_genes.conf_99.txt\', sep = \'/\')']
    runMe += ['amp1 <- paste(\''+path+'\', \'amp_genes.conf_99.txt\', sep = \'/\')']
    runMe += ['sco1<- paste(\''+path+'\', \'scores.gistic\', sep = \'/\')']
    runMe += ['les1 <- paste(\''+path+'\', \'all_lesions.conf_99.txt\', sep = \'/\')']
    # convert matrix to maf format
    runMe += ['plusG = read.maf(maf = maf, gisticAllLesionsFile = les1, gisticAmpGenesFile = amp1, gisticDelGenesFile = del1, gisticScoresFile = sco1, isTCGA = FALSE, cnLevel = \'deep\')']
    runMe += ['pdf(paste(\'plots/oncoplots_deep/\',\''+name+'\',\'_\',\''+cancer+'\',\'_oncoplot.pdf\',sep=\'\'), onefile = TRUE)']
    runMe += ['plot = oncoplot(maf = plusG, genes = c(\''+'\', \''.join([n2.loc[int(i.split('_')[0])] for i in group.index])+'\'), fontSize = 12)']
    #runMe += ['plot = plot + ylab(\'# Mutations in Frequent Cancers\')']
    #runMe += ['print(plot)']
    runMe += ['dev.off()']
    runMe += ['pdf(paste(\'plots/oncostrips_deep/\',\''+name+'\',\'_\',\''+cancer+'\',\'_oncoplot.pdf\',sep=\'\'), onefile = TRUE)']
    runMe += ['plot = oncostrip(maf = plusG, genes = c(\''+'\', \''.join([n2.loc[int(i.split('_')[0])] for i in group.index])+'\'), fontSize = 12)']
    runMe += ['dev.off()']
    runMe = '\n'.join(runMe)+'\n'
    #print(runMe)
    # Run in R
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out = rProc.communicate(runMe.encode())
    return out #Returns errors and other printed statements

cancers = ['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC']

print('Read in all mafs df')
# Only use cols that are required to make a maf with read.maf()
cols = ['Hugo_Symbol','Protein_Change','Chromosome','Start_position','End_position','Reference_Allele', 'Tumor_Seq_Allele2','Variant_Classification', 'Variant_Type','Tumor_Sample_Barcode']
#Create a df (allmafs) with all mutation Annotation data for all cancers
cancermafs = []
for cancer in cancers:
        tmpmaf = pd.read_csv('MutAnno/All_Samples/'+cancer+'_MutAnnos.tsv', index_col = 0, header = 0, sep = '\t')
        tmpmaf.index = [cancer for i in tmpmaf.index]
        cancermafs.append(tmpmaf)

allmafs = pd.concat([i[cols] for i in cancermafs if 'Protein_Change' in i.columns], axis = 0, ignore_index = False)
fails = ['OV', 'KICH', 'LGG', 'KIRC']
#Plot each gene
for cancer in cancers:
    if cancer not in fails:
        #### Get Genes for cancer with top 5 most value added
        recov_10g = recovdf.loc[recovdf['Cancers']==cancer]['Delta'].sort_values().ix[-10:]
        reinf_10g = reinfdf.loc[reinfdf['Cancers']==cancer]['Delta'].sort_values().ix[-10:]
        recov_10g.name = 'recov'
        reinf_10g.name = 'reinf'
        g5s = [recov_10g,reinf_10g]
        # Make OncoPlot for each group in each cancer
        missinggenes = []
        if cancer == 'SKCM':
            path = 'GISTIC_99/gdac.broadinstitute.org_'+cancer+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0'
        else:
            path = 'GISTIC_99/gdac.broadinstitute.org_'+cancer+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0'
        d0 = pd.read_csv(path+'/all_thresholded.by_genes.txt',sep='\t')
        d1 = d0.set_index('Locus ID')
        d1.columns = list(d1.columns[0:2])+[i[:12] for i in d1.columns[2:]]
        n2 = pd.Series(d1['Gene Symbol'], index =d1.index.values)
        for group in g5s:
            genelist = []
            if len(group)>1:
                for gene in group.index:
                    gene0 = gene.split('_')[0]
                    gene1 = gene.split('_')[1]
                    somemafs = allmafs.loc[gene1]
                    somemafs['Cancer']= somemafs.index
                    somemafs = somemafs.set_index('Hugo_Symbol')
                    Symbol = n2.loc[int(gene0)]
                    if Symbol in somemafs.index and len(somemafs.loc[Symbol].shape)>1 :
                        genelist.append(somemafs.loc[Symbol].reset_index())
                        print(Symbol)
                    else:
                        missinggenes.append(gene)
                if len(genelist)>1:
                    genedf = pd.concat(genelist)
                    genedf.to_csv('tmpmaf.tsv',index = 0, sep = '\t')
                    print('plotting '+cancer+' '+group.name)
                    out = PlotOncoPlot(group.name, path, cancer)
                    print(out)



# Enrichment Analysis

all_HyGs = pd.DataFrame()

for cancer in cancers:
    HyGs = pd.DataFrame([],index = [cancer+'_Enrichment of recovered PAMs in Mutsigcv2',cancer+'_Enrichment of all OncoMerged PAMs in Mutsigcv2',cancer+'_Enrichment of recovered CNAdels in sig gistic',cancer+'_Enrichment of recovered CNAamps in sig gistic',cancer+'_Enrichment of GoFs from amps',cancer+'_Enrichment of LoFs from dels',cancer+'_Enrichment of GoFs from PAMs',cancer+'_Enrichment of LoFs from PAMs'], columns = ['intersection','M','n1','n2','p-value'])
    #read in somMuts
    #read in Mutsigcv2
    #PAM Terms
    intallPAMs = set(mutsig.index).intersection(somMuts)
    omPAMs = [n2.loc[int(i)] for i in glofmf[cancer].keys() if 'PAM' in glofmf[cancer][i]]
    sigPAMs = mutsig.loc[mutsig['q']<=0.05].index
    sigPAMs = [i for i in sigPAMs if i in intallPAMs]
    recovPAMs = [n2.loc[int(i.split('_')[0])] for i in statdf.index if statdf.loc[i,'Group']=='recov' and i in intallPAMs]
# Enrichment of recovered PAMs in Mutsigcv2
    HyGs.loc[cancer+'_Enrichment of recovered PAMs in Mutsigcv2'] = [len(set(sigPAMs).intersection(recovPAMs)) ,len(intallPAMs), len(sigPAMs), len(recovPAMs),hyg.sf(len(set(sigPAMs).intersection(recovPAMs)) ,len(intallPAMs), len(sigPAMs), len(recovPAMs))]
    print('recov PAMS: '+str(HyGs.loc[cancer+'_Enrichment of recovered PAMs in Mutsigcv2'].values))
# Enrichment of all OncoMerged PAMs in Mutsigcv2
    HyGs.loc[cancer+'_Enrichment of all OncoMerged PAMs in Mutsigcv2'] = [len(set(sigPAMs).intersection(omPAMs)) ,len(intallPAMs), len(sigPAMs), len(omPAMs), hyg.sf(len(set(sigPAMs).intersection(omPAMs)) ,len(intallPAMs), len(sigPAMs), len(omPAMs))]
    print('OncoMerged PAMs: '+ str(HyGs.loc[cancer+'_Enrichment of all OncoMerged PAMs in Mutsigcv2'].values))
    #CNA Terms
    intallCNAdels = list(deldict.keys())
    intallCNAamps = list(ampdict.keys())
    sigCNAdels = [i for i in deldict.keys() if deldict[i]>=0.05]
    sigCNAamps = [i for i in ampdict.keys() if ampdict[i]>=0.05]
    recovDels = [n2.loc[int(i.split('_')[0])] for i in statdf.loc[statdf['GoF or LoF']=='LoF'].index if statdf.loc[i,'Group']=='recov']
    #recovAmps =
    recovAmps= [n2.loc[int(i.split('_')[0])] for i in statdf.loc[statdf['GoF or LoF']=='GoF'].index if statdf.loc[i,'Group']=='recov']
# Enrichment of recovered CNAs in sig gistic
    HyGs.loc[cancer+'_Enrichment of recovered CNAdels in sig gistic'] = [len(set(sigCNAdels).intersection(recovDels)) ,len(intallCNAdels), len(sigCNAdels), len(recovDels), hyg.sf(len(set(sigCNAdels).intersection(recovDels)) ,len(intallCNAdels), len(sigCNAdels), len(recovDels))]
    print('recov CNAdels: '+str(HyGs.loc[cancer+'_Enrichment of recovered CNAdels in sig gistic'].values))
    HyGs.loc[cancer+'_Enrichment of recovered CNAamps in sig gistic'] = [len(set(sigCNAamps).intersection(recovAmps)) ,len(intallCNAdels), len(sigCNAamps), len(recovAmps), hyg.sf(len(set(sigCNAamps).intersection(recovAmps)) ,len(intallCNAdels), len(sigCNAamps), len(recovAmps))]
    print('recov CNAamps: '+str(HyGs.loc[cancer+'_Enrichment of recovered CNAamps in sig gistic'].values))
    '''# for GoFs from amps...
    # To be clear... This is confusing.
    # Background: Genes that are Pams or Amps
    # Group 1: GoF and Genes that are amps and not Pam
    # Group 2: Genes that are amps and Pams
    # Intersection: GoFs
#Enrichment of GoFs from amp
    HyGs.loc[cancer+'_Enrichment of GoFs from amps'] = [len(set(statdf.loc[statdf['GoF or LoF']=='GoF'])),len(set(intallCNAamps).union(somMuts)),len(set(intallCNAamps).difference(somMuts).union(set(statdf.loc[statdf['GoF or LoF']=='GoF']))), len(set(intallCNAamps).intersection(somMuts)), hyg.sf(len(set(statdf.loc[statdf['GoF or LoF']=='GoF'])),len(set(intallCNAamps).union(somMuts)),len(set(intallCNAamps).difference(somMuts).union(set(statdf.loc[statdf['GoF or LoF']=='GoF']))), len(set(intallCNAamps).intersection(somMuts)))]
    print('GoFs from amps: '+str(HyGs.loc[cancer+'_Enrichment of GoFs from amps'].values))
    # for LoFs from dels...
    # To be clear... This is confusing.
    # Background: Genes that are Pams or Dels (union)
    # Group 1: LoFs and dels that are not Pams
    # Group 2: dels and Pams (intersection)
    # Intersection: LoFs
#Enrichment of LoFs from dels
    HyGs.loc[cancer+'_Enrichment of LoFs from dels'] = [len(set(statdf.loc[statdf['GoF or LoF']=='LoF'])),len(set(intallCNAdels).union(somMuts)),len(set(intallCNAdels).difference(somMuts).union(set(statdf.loc[statdf['GoF or LoF']=='LoF']))), len(set(intallCNAdels).intersection(somMuts)), hyg.sf(len(set(statdf.loc[statdf['GoF or LoF']=='LoF'])),len(set(intallCNAdels).union(somMuts)),len(set(intallCNAdels).difference(somMuts).union(set(statdf.loc[statdf['GoF or LoF']=='LoF']))), len(set(intallCNAdels).intersection(somMuts)))]
    print('LoFs from dels: '+str(HyGs.loc[cancer+'_Enrichment of LoFs from dels'].values))
    # for GoFs from PAMs...
    # To be clear... This is confusing.
    # Background: Genes that are Pams or Amps
    # Group 1: GoF and Genes that are Pams and not amps
    # Group 2: Genes that are amps and Pams
    # Intersection: GoFs
#Enrichment of GoFs from PAMs
    HyGs.loc[cancer+'_Enrichment of GoFs from PAMs'] = [len(set(statdf.loc[statdf['GoF or LoF']=='GoF'])),len(set(intallCNAamps).union(somMuts)),len(set(somMuts).difference(intallCNAamps).union(set(statdf.loc[statdf['GoF or LoF']=='GoF']))), len(set(intallCNAamps).intersection(somMuts)), hyg.sf(len(set(statdf.loc[statdf['GoF or LoF']=='GoF'])),len(set(intallCNAamps).union(somMuts)),len(set(somMuts).difference(intallCNAamps).union(set(statdf.loc[statdf['GoF or LoF']=='GoF']))), len(set(intallCNAamps).intersection(somMuts)))]
    print('GoFs from amps: '+str(HyGs.loc[cancer+'_Enrichment of GoFs from PAMs'].values))
    # for LoFs from PAMs...
    # To be clear... This is confusing.
    # Background: Genes that are Pams or Dels (union)
    # Group 1: LoFs and PAMs that are not dels
    # Group 2: dels and Pams (intersection)
    # Intersection: LoFs
    #Enrichment of LoFs from PAMs
    HyGs.loc[cancer+'_Enrichment of LoFs from dels'] = [len(set(statdf.loc[statdf['GoF or LoF']=='LoF'])),len(set(intallCNAdels).union(somMuts)),len(set(somMuts).difference(intallCNAdels).union(set(statdf.loc[statdf['GoF or LoF']=='LoF']))), len(set(intallCNAdels).intersection(somMuts)), hyg.sf(len(set(statdf.loc[statdf['GoF or LoF']=='LoF'])),len(set(intallCNAdels).union(somMuts)),len(set(somMuts).difference(intallCNAdels).union(set(statdf.loc[statdf['GoF or LoF']=='LoF']))), len(set(intallCNAdels).intersection(somMuts)))]
    print('LoFs from dels: '+str(HyGs.loc[cancer+'_Enrichment of LoFs from dels'].values))'''
    all_HyGs = pd.concat([all_HyGs,HyGs])


all_HyGs.to_csv('output/HyperGeometric/HyG_1.csv', index = 1)

# Plot heat map:

import numpy as np
from math import log10
cutdels = {}
cutamps = {}
with PdfPages('Heatmap_GLoF_Enrichment.pdf') as pdf:
#residual q must be larger than 0.05, to be considered here. also, we're going to adjust the mf threshold for what dels are considered. How?
    threshs = np.arange(1,11)/100
    headings = [str(i) for i in threshs]
    for cancer in cancers:
        cutdels[cancer] = pd.DataFrame([],index =headings, columns=headings )
        cutdels[cancer].index.name = 'PAMs'
        cutamps[cancer] = pd.DataFrame([],index =headings, columns= headings )
        cutamps[cancer].index.name = 'PAMs'
        CNAdels = [gene for gene in glofmf[cancer] if 'CNAdel' in glofmf[cancer][gene]]
        CNAamps= [gene for gene in glofmf[cancer] if 'CNAamp' in glofmf[cancer][gene]]
        PAMs= [gene for gene in glofmf[cancer] if 'PAM' in glofmf[cancer][gene]]
        for thresh1 in threshs:
            Athresh = thresh1
            Dthresh = thresh1
            for thresh2 in threshs:
                Pthresh = thresh2
                f_CNAdels = [i for i in CNAdels if glofmf[cancer][i]['CNAdel'] >= Dthresh]
                f_CNAamps= [i for i in CNAamps if glofmf[cancer][i]['CNAamp'] >= Athresh]
                f_PAMs= [i for i in PAMs if glofmf[cancer][i]['PAM'] >= Pthresh]
# Max genes should probably only consider all genes from summuts and gistic
                cutdels[cancer].loc[str(Pthresh), str(Dthresh)] = hyg.sf(len(set(f_CNAdels).intersection(f_PAMs)), 22500, len(f_PAMs) ,len(f_CNAdels))
                print('cutdels: ', str(len(set(f_CNAdels).intersection(f_PAMs))), str(22500), str(len(f_PAMs)) ,str(len(f_CNAdels)))
                cutamps[cancer].loc[str(Pthresh), str(Athresh)] = hyg.sf(len(set(f_CNAamps).intersection(f_PAMs)), 22500, len(f_PAMs) ,len(f_CNAamps))
        ax = sns.heatmap(cutamps[cancer].applymap(lambda x: -log10(x) if x>0 else np.nan), vmin = 0, vmax = 100, annot = True)
        plt.title(cancer+ ' amps')
        pdf.savefig()
        plt.close()
        bx = sns.heatmap(cutdels[cancer].applymap(lambda x: -log10(x) if x>0 else np.nan), vmin = 0, vmax = 100,annot = True,fmt = 'g')
        plt.title(cancer+' dels')
        pdf.savefig()
        plt.close()






# Plot Ven Diagram of Genes from OncoMerge, GISTIC, and MutSig2cv:

from matplotlib_venn import venn3

with PdfPages('Mutation_Source_Comparison.pdf') as pdf:
    OncoMerge = pd.Series({cancer: len(OM_Genes[cancer]) for cancer in cancers}, name = 'OncoMerge Genes').sort_values()
    OncoMerge_Gistic = pd.Series({cancer: len(OM_Genes[cancer].intersection(Gistic_Genes[cancer])) for cancer in cancers}, name = 'OncoMerge Overlap with Gistic Genes')
    OncoMerge_MutSig2CV = pd.Series({cancer: len(OM_Genes[cancer].intersection(MutSig2CV_Genes[cancer])) for cancer in cancers}, name = 'OncoMerge Overlap with MutSig2CV genes')
    percents = pd.concat([OncoMerge, OncoMerge_Gistic, OncoMerge_MutSig2CV]).T
    percents.plot.bar()
    plt.ylabel('Percent of Overlap OncoMerge Genes')
    pdf.savefig()
    plt.close()
    for cancer in cancers:
        print( "Making VennDiagram for:  "+str(cancer))
        ax1 = plt.gca()
        l = ['OncoMerge Genes','MutSig2cv Genes','GISTIC2 Genes']
        v = venn3([OM_Genes[cancer],MutSig2CV_Genes[cancer], Gistic_Genes[cancer]],tuple(l))
        #ax.legend(handles=h, labels=l, title="counts")
        plt.title(cancer)
        pdf.savefig()
        plt.close()




