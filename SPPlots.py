import sys
import os
script_path = os.path.dirname(os.path.realpath(__file__)).split('SPTools')[0]
sys.path.append(script_path)
import SPTools as SP
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

##################################################################################
## Tools to use with SP_pipeline and SP_quant scripts. Use quant dataframe      ##
##################################################################################

def draw_diagonal(ax):
    ax_max = max([ax.get_ylim()[1],ax.get_xlim()[1]])
    ax_min = min([ax.get_ylim()[0],ax.get_xlim()[0]])
    ax.plot([ax_min, ax_max],[ax_min, ax_max],'--',color='0.7')
    ax.set_xlim([ax_min, ax_max])
    ax.set_ylim([ax_min, ax_max])
    return ax, (ax_min, ax_max)
    
def SP_pipeline_scatters(df, base_dir, name):
    metrics = [x.split(' avg')[0] for x in df.columns if 'avg' in x]
    colors = ['0.3','royalblue','mediumpurple','deepskyblue','0.5','0.7']

    if len(metrics) <= 4:
            figsize = (8,8)
    else: figsize = (8,12)
    fig, ax = plt.subplots(nrows=len(metrics)/2+len(metrics)%2, ncols=2, figsize=figsize)
    
    for n, metric in enumerate(metrics):
        columns = [x for x in df.columns if metric in x]
        columns = [x for x in columns if 'avg' not in x]

        data = df[columns].apply(np.log2).replace([np.inf, -1*np.inf],np.NaN).dropna(how='any')
        if metric == 'Total Spliceosomal' or metric == 'RNAseq':
            data = data.drop_duplicates(keep='first')
        
        data1 = data[columns[0]]
        data2 = data[columns[1]]
    
        ax[n/2][n%2].scatter(data1, data2, color=colors[n], alpha=0.6, s=10)
        ax[n/2][n%2].set_xlabel(columns[0]+' (log2)')
        ax[n/2][n%2].set_ylabel(columns[1]+' (log2)')
        ax[n/2][n%2], limits = draw_diagonal(ax[n/2][n%2])
        ax[n/2][n%2].set_title(metric)
        
        r, p = stats.pearsonr(data1, data2)
        xy=((limits[1]+limits[0])/2., limits[0]+abs(0.05*limits[0]))
        ax[n/2][n%2].annotate("Pearson R coef: "+"%0.2f" %r, xy=xy, fontsize=12)
        
    if len(metrics) == 3:
        ax[1][1].axis('off')
    
    fig.tight_layout()
    plt.show()
    fig.savefig(base_dir+name+'_scatterplots.pdf',format='pdf')
    
def SP_quant_scatters(final_df, bam_dict, W=False):
    reps = []
    for genotype, sample in bam_dict.iteritems():
        sample1 = (sample['A1'].split('/')[-1].split('-A')[0])
        sample2 = (sample['A2'].split('/')[-1].split('-A')[0])
        reps.append((genotype,sample1,sample2))
    
    # Build combos of replicates and comparisons between genotypes
    combos = [reps[0], reps[1],
             ((reps[0][0],reps[1][0]), reps[0][1], reps[1][1]),
             ((reps[0][0],reps[1][0]), reps[0][2], reps[1][2])]
    
    color_list = ['0.3','royalblue','mediumpurple','deepskyblue','0.5','0.7']
    metrics = [' Intermediate Level',' Total Spliceosomal',' Precursor', '-B Intron Retention']
    if W is True:
        metrics.append(' RNAseq')
    
    for k, metric in enumerate(metrics):
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
        
        for n, (genotype, sample1, sample2) in enumerate(combos):
            if type(genotype) == str:
                final_df.loc[:,(genotype,sample1+' log2'+metric)] = final_df[(genotype,sample1+metric)].apply(np.log2)
                final_df.loc[:,(genotype,sample2+' log2'+metric)] = final_df[(genotype,sample2+metric)].apply(np.log2)
                final_df = final_df.replace([np.inf, np.inf*-1], np.NaN).dropna(how='any')
                
                data1 = final_df.loc[:,(genotype,sample1+' log2'+metric)]
                data2 = final_df.loc[:,(genotype,sample2+' log2'+metric)]
                ax[n/2][n%2].set_title(genotype+' replicates', fontsize=14)
        
            else:
                final_df.loc[:,(genotype[0],sample1+' log2'+metric)] = final_df[(genotype[0],sample1+metric)].apply(np.log2)
                final_df.loc[:,(genotype[1],sample2+' log2'+metric)] = final_df[(genotype[1],sample2+metric)].apply(np.log2)
                final_df = final_df.replace([np.inf, np.inf*-1], np.NaN).dropna(how='any')
                
                data1 = final_df.loc[:,(genotype[0],sample1+' log2'+metric)]
                data2 = final_df.loc[:,(genotype[1],sample2+' log2'+metric)]
                ax[n/2][n%2].set_title(genotype[0]+' vs. '+genotype[1], fontsize=14)
        
            all_data = data1.append(data2)
            all_max = max(all_data)+0.1*max(all_data)
            all_min = min(all_data)-0.1*abs(min(all_data))
        
            ax[n/2][n%2].scatter(data1,data2,color=color_list[k],alpha=0.5,s=10)
            ax[n/2][n%2], limits = draw_diagonal(ax[n/2][n%2])
            ax[n/2][n%2].set_xlabel(sample1+' log2'+metric, fontsize=12)
            ax[n/2][n%2].set_ylabel(sample2+' log2'+metric, fontsize=12)
            
            r, p = stats.pearsonr(data1, data2)
            xy=((limits[1]+limits[0])/2., limits[0]+abs(0.05*limits[0]))
            ax[n/2][n%2].annotate("Pearson R coef: "+"%0.2f" %r, xy=xy, fontsize=12)
        
        fig.tight_layout()    
        plt.show()
        plt.clf()
