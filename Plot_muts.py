import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import time as t

'''
parser = argparse.ArgumentParser(description='Used to insert a mf value.')
parser.add_argument('--mf', help='Frequency of mutation threshold. Needs to have a CANCER_finalMutFile_deep_filtered_mf_mf. Please input the cooresponding mf as either 0.05 or 0.0', type = str)
args = parser.parse_args()
print(args)
'''

# Set up

start = t.time()
cancers = sorted(['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC'])
sms = ['PAM','amp','GoF','del','LoF']
ind = range(len(sms))
dict1 = {}
mfs = [0.05]
mf = '0.05'
allthings = pd.DataFrame([],columns = ['mf', 'type','Cancers'])

# Plot Bar plot and raw csv of number of each CNAtype in OncoMerged file

for mf in mfs:
    dict1[str(mf)] = {}
    with PdfPages('plots/1_7/mf_'+str(mf)+'_bar.pdf') as pdf:
        for cancer in cancers:
            d1 = pd.read_csv('output/finalMutFile_deep_filtered/region/'+cancer+'_finalMutFile_deep_filtered_mf_'+str(mf)+'.csv', index_col = 0, header = 0)
            d1.index = [i+'_'+cancer for i in d1.index]
            df1 = pd.DataFrame(dict(zip(['mf','type'],[d1.T.mean(),pd.Series(dict(zip(list(d1.index),[i.split('_')[1] for i in list(d1.index)])))])))
            dict1[str(mf)][cancer] = {}
            df1['Cancers'] = [cancer for i in df1.index]
            allthings = allthings.append(df1)
            for sm in sms:
                if sm==sms[-1]:
                    dict1[str(mf)][cancer][sm] = sum(d1.loc[[i for i in d1.index.values if sm in i]].T.mean()>mf) # count of genes for each cna type
                    total = len(set([i for i in d1.index.values if  d1.loc[i].mean()>mf]))
                    #total = len(set([i.split('_')[0] for i in dict1.index.values]))
                    stats = list(dict1[str(mf)][cancer].values())
                    plt.bar(ind,stats, color = ['r','b','g','c','m'])
                    plt.ylim(0, 6000)
                    plt.xticks(ind, sms)
                    plt.title(cancer+', Total Genes: '+str(total))
                    pdf.savefig()
                    plt.close()
                dict1[str(mf)][cancer][sm] = sum(d1.loc[[i for i in d1.index.values if sm in i]].T.mean()>mf)

        pd.DataFrame(dict1).T.reindex(columns = sms).to_csv('output/barraw/'+cancer+'_mf_'+str(mf)+'_raw.csv')

'''
for cancer in cancers:
    d1 = pd.read_csv('output/finalMutFile_deep_filtered/region/'+cancer+'_finalMutFile_deep_filtered_mf_'+args.mf+'.csv', index_col = 0, header = 0)
    d1.index = [i+'_'+cancer for i in d1.index]
    df1 = pd.DataFrame(dict(zip(['mf','type'],[d1.T.mean(),pd.Series(dict(zip(list(d1.index),[i.split('_')[1] for i in list(d1.index)])))])))
    df1['Cancers'] = [cancer for i in df1.index]
    allthings = allthings.append(df1)
    dict1[cancer] = {}
    for sm in sms:
        dict1[cancer][sm] = sum(d1.loc[[i for i in d1.index.values if sm in i]].T.mean()>0.0)
'''

#pd.DataFrame(dict1).T.reindex(columns = sms).to_csv('mf_'+args.mf+'_raw.csv')
endbar = t.time()
print("Time to get data: "+str(endbar-start))

# Plot Histograms to plot the distribution of cancers for the number of different CNAtypes

bins = 11
h = {}
for key in dict1.keys():
    with PdfPages('plots/histograms/'+key+'_histogram.pdf') as pdf:
        c1 = pd.DataFrame(dict1[key]).reindex(['PAM','amp','GoF','del','LoF']).T
        c1.columns = ['PAM','amp','GoF','del','LoF']
        m1 = c1.max()
        m1.name = 'Max'
        mean1 = c1.mean()
        mean1.name = 'Mean'
        median1 = c1.median()
        median1.name = 'Median'
        mode1 = c1.mode().loc[0]
        mode1.name = 'Mode'
        kurt1 = c1.kurtosis()
        kurt1.name = 'Kurtosis'
        skew1 = c1.skew()
        skew1.name = 'Skew'
        c2 = c1.append([mean1,median1,mode1,kurt1,skew1,m1])
        h[key] = {}
        for sm in sms:
            h[key][sm] = plt.hist(c1[sm].T, bins = 11)
            plt.xlabel('Number of CNAs')
            plt.ylim(0,12)
            plt.xlim(0,6000)
            plt.ylabel('Number of Cancers')
            plt.title('Histogram of '+sm+' for mf = '+key)
            pdf.savefig()
            plt.close()
        app1 = pd.concat([pd.Series(h[key][sm][0], index = ['Interval_'+str(i) for i in range(len(h[key][sm][1])-1)],name = sm) for sm in sms],axis = 1, sort=False)
        app1.index = ['Interval_'+str(i) for i in range(len(app1))]
        c3 = c2.append(app1)
        c3.to_csv('mf_'+key+'_raw.csv')



endhist = t.time()
print('TIme to plot Cancer Swarms for each CNA Type: '+str(endbar-endhist))


# Plot Violin Plots

dict2 = {}
for sm in sms:
    dict2[sm] = {}
    for cancer in dict1[str(mf)]:
        dict2[sm][cancer] = dict1[str(mf)][cancer][sm]
#Plots one plot for each CNA type listed in sms. Each plot will have one swarm per cancer type included.

dv = pd.DataFrame(dict2)
dv['Cancers']=dv.index.values
dv.sort_values('PAM')
smsallthings = {sm:allthings.loc[[i for i in allthings.index if sm in i]] for sm in sms}
with PdfPages('plots/1_7/violin_allCancers_mf_'+str(mf)+'_.pdf') as pdf:
    for i in range(len(sms)):
        start_CNA = t.time()
        ax = sns.violinplot(x = 'Cancers', y = 'mf',scale = 'count',inner = None, order = list(dv.sort_values('PAM').index) , data=smsallthings[sms[i]],axis = 1)
        #ax = sns.swarmplot(x = 'Cancers', y = 'mf', order = list(dv.sort_values('PAM').index) , data=smsallthings[sms[i]])#,axis = 1)
        plt.title(sms[i])
        plt.ylim(float(mf), 1)
        plt.xticks(rotation=90)
        pdf.savefig()
        plt.close()
        end_CNA = t.time()
        print('Time to plot '+sms[i]+': '+str(end_CNA-start_CNA))
        print('Number of '+sms[i]+'\'s: '+str(len(smsallthings[sms[i]])))
        '''
        print(sms[i]+': width')
        ax = sns.v1iolinplot(x = 'Cancers', y = 'mf', order = list(dv.sort_values('PAM').index), inner = None, scale = 'width', data=smsallthings[sms[i]],axis = 1)
        plt.title(sms[i]+ ' width')
        plt.xticks(rotation=90)
        pdf.savefig()
        plt.close()
        print(sms[i]+': area')
        ax = sns.v1iolinplot(x = 'Cancers', y = 'mf', order = list(dv.sort_values('PAM').index), inner = None, scale = 'area', data=smsallthings[sms[i]],axis = 1)
        plt.title(sms[i]+ ' area')
        plt.xticks(rotation=90)
        pdf.savefig()
        plt.close()
        '''

endallC = t.time()
print('TIme to plot Cancer Swarms for each CNA Type: '+str(endallC-endhist))


#Plots one plot for each Cancer type listed in cancers. Each plot will have one violin per CNA type included.

with PdfPages('plots/1_7/Swarm_mf_'+str(mf)+'_.pdf') as pdf:
    for cancer in cancers:
        start_cancer = t.time()
        data = allthings.reindex(index = [i for i in allthings.index if cancer in i])
        ax = sns.violinplot(x = 'type', y = 'mf', scale = 'area', inner = None, order = ['PAM','GoF','CNAamp','LoF','CNAdel'],data = allthings.reindex(index = [i for i in allthings.index if cancer in i]))
        #ax = sns.swarmplot(x = 'type', y = 'mf', color = 'k', size=2, order = ['PAM','GoF','CNAamp','LoF','CNAdel'],data=data)#, ax = ax)
        plt.ylim(float(mf), 1)
        plt.title('Frequency of Mutation in ' + cancer)
        pdf.savefig()
        plt.close()
        end_cancer = t.time()
        print('Time to make '+cancer+': '+str(end_cancer-start_cancer))
        print('Number of '+cancer+': '+str(len(data)))

end = t.time()
print('Time to plot CNA Type Swarms for each Cancer: '+str(end-endallC))
print('Total Time for mf '+str(mf)+': '+str(end-start))



