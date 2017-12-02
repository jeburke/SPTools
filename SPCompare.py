import sys
import os
script_path = os.path.dirname(os.path.realpath(__file__)).split('SPTools')[0]
sys.path.append(script_path)
import SPTools as SP
import pandas as pd
import numpy as np
import pysam
import copy
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import colors
from sklearn import cluster
from statsmodels.stats import proportion

def main():
    '''Usage: run SPCompare.py Sample_config.txt WT_quantitation.csv organism
    Each line in config file will be : bam_file,genotype,sample (e.g. CM763-A_sorted.bam,WT,A1)
    WT_quantitation.csv is the output from SP_pipeline.py for the wild type sample (indicates where high confidence events are)
    organism is a string - can be 'cyrpto', 'pombe' or 'cerevisiae'
    '''
    
    bam_dict = {}
    with open(sys.argv[1], 'r') as config:
        for line in config:
            info = line.split(',')
            genotype = info[1]
            sample = info[2].strip()
            
            if genotype not in bam_dict:
                bam_dict[genotype] = {}
            
            bam_dict[genotype][sample] = info[0]
    
    prefix = sys.argv[1].split('_config')[0]
    
    organism = sys.argv[3]
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    
    columns = ['5p score','exon size (us)','exon size (ds)','introns in transcript','type','transcript size','intron size',
               'chromosome','position','alt splicing','3p score','transcript','intron position','strand','peak','branch score',
               'branch to 3p distance','percent pPy','seq5','seq3']
    
    quant_df = pd.read_csv(sys.argv[2], index_col=0)
    try:
        quant_df = quant_df[columns]
    except KeyError:
        print "Columns missing from dataframe..."
        print [x for x in columns if x not in quant_df.columns]
        return None
    
    final_df = copy.deepcopy(quant_df)
    final_df.columns = pd.MultiIndex.from_product([['Peaks'], final_df.columns])

    for genotype, samples in bam_dict.iteritems():
        # Determine if whole cell extract samples are present
        Ws = [x for x in samples.keys() if "W" in x]
        if len(Ws) > 1: 
            W = True
        else: 
            W=False
        
        # Quantitate all samples with genotype
        new_df = SP.quantitate_junction_df(samples, quant_df, gff3, W=W)
        
        # Remove original columns and rename new ones with multiindex
        new_columns = [x for x in new_df.columns if x not in columns]
        new_df = new_df[new_columns]
        new_df.columns = pd.MultiIndex.from_product([[genotype], new_df.columns])
        final_df = final_df.join(new_df, how='inner')
        
    final_df.to_csv(prefix+'_quant_df.csv')
    
if __name__ == "__main__":
    main()
    

    
## Heatmap function used for prp43DN in paper

def log_ratio_heatmap(df, all_df, wt_name, mut_name, metrics = ['Precursor','Intermediate Level'], n_clusters=4):
    '''Generates a heatmap clustered using k-means and runs Mann Whitney U or Proportion hypothesis test for intron features in data frame
        Two modes - the first creates an "elbow" graph to determine how many centroids to use (look for inflection point).
                    the second does the final clustering, generates the heatmap and p-values comparing between clusters.
    
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
         Output from main function of SPCompare
    all_df : pandas.core.frame.DataFrame
         Dataframe from finding high confidence events from WT, output of SP_Pipeline
    wt_name : str, name of wild type sample in df (top column name)
    mut_name : str, name of mutant sample in df (top column name)
    metrics : list, default ['Precursor','Intermediate Level']
        Measurements to include in heatmap. Other optionsa are: 'W Intron Retention', 'B Intron Retention', 
        'Total Spliceosomal' and 'RNAseq'
    n_clusters : int, default 4
        Number of centroids for k-means clustering. Run test mode by setting to None.
        
    Returns
    ------
    fig : matplotlib.figure.Figure, Heatmap figure object that can be saved using fig.savefig()
    new_df : pandas.core.frame.DataFrame, New dataframe with cluster numbers
    all_p : pandas.core.frame.DataFrame, dataframe with log10 p-values for intron features for each cluster. 
            Sign of the p-value indicates whether the intron feature scores better (+) or worse (-) than the general population
        '''
    
    new_df = copy.deepcopy(df)
    print len(new_df)
    
    # Filter dataframe for peaks with at least 5 reads in intermediate level and precursor
    filt_cols = [x for x in new_df.columns if ('-A' in x[1]) and (x[0] == 'WT')]
    filt_cols = filt_cols + [x for x in new_df.columns if ('near peak' in x[1]) and (x[0] == 'WT')]
    new_ix = new_df.index
    for col in filt_cols:
        new_ix = set(new_df[new_df[col] >= 5].index).intersection(new_ix)
    
    # Apply filter to dataframe and calculate log ratio averages
    new_df = new_df[new_df.index.isin(new_ix)]
    for metric in metrics:
        new_df[('All',metric+' avg')] = new_df[('All',metric+' log2 ratio1')]+new_df[('All',metric+' log2 ratio2')].divide(2.)
    print len(new_df)

    columns = [x for x in new_df.columns if (x[0] == 'All') & ('change' not in x[1])]
    new_df = new_df[columns]
    new_df = new_df.dropna(how='any')

    # Get min/max for heatmap
    min_max = []
    for column in columns:
        min_max.append(abs(max(new_df[column])))
        min_max.append(abs(min(new_df[column])))
    both_max = max(min_max)-1

    # K-means clustering
    if n_clusters is not None:
        mat = new_df[[x for x in new_df.columns if 'avg' in x[1]]].as_matrix()
        km = cluster.KMeans(n_clusters=n_clusters)
        km.fit(mat)
        labels = km.labels_
        new_df.loc[:,('All','cluster')] = labels
        cluster_names = list(set(labels))
        new_df = new_df.sort_values([('All','cluster')])
    else:
        cluster_range = range(2,13)
        final_n = None
        sos_diff = []
        for n_clusters in cluster_range:
            test_df = copy.deepcopy(new_df)
            mat = test_df[[x for x in test_df.columns if 'avg' in x[1]]].as_matrix()
            km = cluster.KMeans(n_clusters=n_clusters)
            km.fit(mat)
            labels = km.labels_
            test_df.loc[:,('All','cluster')] = labels
            cluster_names = list(set(labels))
            test_df = test_df.sort_values([('All','cluster')])
            
            sos = []
            for index in range(n_clusters):
                for column in [x for x in test_df.columns if 'avg' in x[1]]:
                    n=0
                    reps = []
                    while n<10:
                        reps.append(sum(map(lambda x:x*x, test_df[test_df[('All','cluster')] == index][column])))
                        n += 1
                    sos.append(sum(reps)/float(len(reps)))

            n=0
            cluster_diff = 0
            for n in range(len(sos)-1):
                cluster_diff += abs(sos[n+1]-sos[n])
            sos_diff.append(cluster_diff)
            
        plt.plot(cluster_range,sos_diff, 'o', c='k')
        plt.show()
        plt.clf()
        return None, None, None

    # Calculate p-values for different intron properties and retrieve cluster sizes
    cluster_sizes = []
    all_p = pd.DataFrame(index=set(new_df[('All','cluster')]), columns=['intron size', '5p score', '3p score', 'branch score','percent pPy', 'exon size (ds)','exon size (us)','transcript size',
                      'introns in transcript','branch to 3p distance', ('alt splicing:True'),'alt splicing:False','intron position:First', 
                                                   'intron position:Middle', 'intron position:Last'])
    for index in range(n_clusters):
        cluster_df = new_df[new_df[('All','cluster')] == index]
        cluster_sizes.append(len(cluster_df))
        for x_var in all_p.columns:
            if x_var.split(':')[0] in ['alt splicing','intron position']:
                value = x_var.split(':')[1]
                x_var2 = x_var.split(':')[0]
                count = len(df[(df.index.isin(cluster_df.index)) & (df[('Peaks',x_var2)].apply(str) == value)])
                nobs = len(df[df.index.isin(cluster_df.index)])
                pop_count = len(all_df[~(all_df.index.isin(cluster_df.index)) & (all_df[x_var2].apply(str) == value)])
                pop_obs = len(all_df[~all_df.index.isin(cluster_df.index)])
                pop_prop = float(pop_count)/pop_obs
                z, p = proportion.proportions_ztest(count, nobs, value=pop_prop, alternative='larger')
                if nobs/float(count) > pop_prop:
                    all_p.loc[index, x_var] = -1*np.log10(p)
                else:
                    all_p.loc[index, x_var] = np.log10(p)
            else:
                ks_t, p = stats.mannwhitneyu(all_df[~all_df.index.isin(cluster_df.index)][x_var], 
                                     df[df.index.isin(cluster_df.index)][('Peaks',x_var)])
                if np.median(df[df.index.isin(cluster_df.index)][('Peaks',x_var)]) > np.median(all_df[~all_df.index.isin(cluster_df.index)][x_var]):
                    all_p.loc[index, x_var] = -1*np.log10(p)
                else:
                    all_p.loc[index, x_var] = np.log10(p)
                
    for column in all_p.columns:
        all_p[column] = pd.to_numeric(all_p[column])
    
    new_df= new_df.drop([x for x in new_df.columns if ('avg' in x[1]) or ('Retention' in x[1])], axis=1)

    fig = plt.figure(figsize=(3,6))
    ax = fig.add_subplot(111)
    print len(new_df.columns)
    pcm = ax.pcolor(range(len(new_df.columns)), range(len(new_df.index)), new_df[[x for x in new_df.columns if x[1] != 'cluster']], 
                    cmap='RdBu_r',norm=colors.Normalize(vmin=both_max*-1, vmax=both_max))
    base = 0
    for n in range(n_clusters):
        base += cluster_sizes[n]
        ax.plot([0,len(new_df.columns)],[base,base],'-',color='0.3')
    
    ax.set_xlim(0,len(new_df.columns)-1)
    labels = [x[1].split(' log2 ratio')[0] for x in new_df.columns if x[1] != 'cluster']
    ax.set_xticks(np.arange(len(new_df.columns)-1)+0.5)
    ax.set_xticklabels(labels, rotation='vertical')
    fig.colorbar(pcm)
    
    new_df = df.join(new_df[('All','cluster')], how='left')
    new_df = new_df[new_df[('All','cluster')].apply(str) != 'nan']
    
    return fig, new_df, all_p