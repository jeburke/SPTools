__author__ = 'jordanburke'

import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome')
sys.path.append('/home/jordan/RNA-is-awesome')
import SPTools as SP
import pandas as pd
import pysam
from subprocess import check_output
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import random
from beeswarm import *
from scipy.stats import ks_2samp
from scipy import stats
import re
import operator
import os
import seaborn as sns
sns.set_style('white')

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(dataFrame, samplename, count_type, log=False, base=10, T0only=False):
    values = []
    for name1, name2 in dataFrame.columns:
        if name2 == count_type:
            if name1 == samplename:
                s = pd.Series(dataFrame[(name1,name2)])
                if T0only == True:
                    transcript_list = []
                    for row in dataFrame.iterrows():
                        if isinstance(row[0], tuple):
                            if row[0][0].endswith("T0"):
                                transcript_list.append(row[0][0])
                        else:
                            if row[0].endswith("T0"):
                                transcript_list.append(row[0])
                    s = s[s.index.get_level_values(0).isin(transcript_list)]
                values = s.tolist()
            else:
                continue
    if log==True:
        log_list = []
        for y in values:
            if float(y) == float('Inf'):
                log_list.append(np.NaN)
            elif y == 0:
                log_list.append(np.NaN)
            else:
                log_list.append(math.log(y, base))
        values = log_list
    return values

##########################################################################
## Count splice sites or transcripts with at least some number of reads ##
##########################################################################

def count_above_cutoff(dataframe, sample_tuple, cutoff=5):
    counter1 = 0
    counter2 = 0

    for name1, name2 in dataframe.columns:
        if name1 == sample_tuple[0] and name2 == sample_tuple[1]:
            print name1, name2
            s = pd.Series(dataframe[(name1,name2)])
            #print s
            for row in s.iteritems():
                #print row
                #print row.index
                if isinstance(row[0], tuple):
                    if row[0][0].endswith("T0"):
                        counter2 +=1
                        if row[1] > cutoff:
                            counter1 +=1                                   
                else:
                    if row[0].endswith("T0"):
                        #print row [0]
                        #print row [1]
                        counter2 +=1
                        if row[1] > cutoff:
                            counter1 +=1
                            
    if sample_tuple[1].endswith("prime"):
        print "Number of "+sample_tuple[1]+" splice sites with at least "+str(cutoff)+" reads: "+str(counter1)
        print "Number of "+sample_tuple[1]+" splice sites: "+str(counter2)
        print str(float(counter1)/float(counter2)*100.)+"% of splice sites"
    elif sample_tuple[1] == "Total":
        print "Number of transcripts with at least "+str(cutoff)+" reads: "+str(counter1)
        print "Number of "+sample_tuple[1]+" transcripts: "+str(counter2)
        print str(float(counter1)/float(counter2)*100.)+"% of transcripts"

#################################################################################                    
## Function that determines if values in the 1st sample are higher or lower    ##
##than values in the 2nd sample (2 fold)                                       ##
#################################################################################

def compare_samples(df, SampleTuple1, SampleTuple2, fold_change=10):
    count_type1 = SampleTuple1[1]
    count_type2 = SampleTuple2[1]
    samplename1 = SampleTuple1[0]
    samplename2 = SampleTuple2[0]
    
    for name1, name2 in df.columns:
        if name2 == count_type1 and name1 == samplename1:
            s = pd.Series(df[(name1,name2)])
            dict1 =  s.to_dict()
        elif name2 == count_type2 and name1 == samplename2:
            s = pd.Series(df[(name1,name2)])
            dict2 = s.to_dict()
            
    low_dict = {}
    low_counter = 0
    high_dict = {}
    high_counter = 0
    for name, value in dict1.iteritems():
        cutoff_high = float(dict2[name])*fold_change
        cutoff_low = float(dict2[name])/fold_change
        if name in dict2:
            if float(dict1[name]) > cutoff_high and dict2[name]!=0 and dict1[name]!=0:
                high_dict[name] = []
                high_dict[name].append(dict1[name])
                high_dict[name].append(dict2[name])
                high_dict[name].append(dict1[name]/dict2[name])
                high_counter += 1
            elif float(dict1[name]) < cutoff_low and dict2[name]!=0 and dict1[name]!=0:
                low_dict[name] = []
                low_dict[name].append(dict1[name])
                low_dict[name].append(dict2[name])
                low_dict[name].append(dict1[name]/dict2[name])
                low_counter += 1
    
    with open("{0}_{1}_comparison.tsv".format(SampleTuple1[0],SampleTuple2[0]), "w") as fout:
        fout.write("Transcript\t Intron\t Reads in "+SampleTuple1[0]+"\t Reads in "+SampleTuple2[0]+"\t Ratio\n")
        fout.write("Higher in "+SampleTuple1[0]+" than "+SampleTuple2[0]+"\n")
        for name, value in high_dict.iteritems():
            line_list = [name[0], str(name[1]), str(value[0]), str(value[1]), str(value[2]), "\n"]
            line = "\t".join(line_list)
            fout.write(line)
        fout.write("Lower in "+SampleTuple1[0]+" than "+SampleTuple2[0]+"\n")
        for name, value in low_dict.iteritems():
            line_list = [name[0], str(name[1]), str(value[0]), str(value[1]), str(value[2]), "\n"]
            line = "\t".join(line_list)
            fout.write(line)
    print str(high_counter)+" introns higher in "+SampleTuple1[0]+" than "+SampleTuple2[0]
    print str(low_counter)+" introns lower in "+SampleTuple1[0]+" than "+SampleTuple2[0]
    return (high_dict, low_dict)

def compare_dicts(dict1, dict2):
    both_dict = {}
    for name, values in dict1.iteritems():
        if name in dict2:
            both_dict[name] = values
    print len(both_dict)
    return both_dict
            
###################################################################
## Plot ratios (usually log) of replicates for normalized counts ##
###################################################################


def scatter_plot(xvalues1, yvalues1, xvalues2, yvalues2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)'):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
    ax1.scatter(xvalues2, yvalues2, c='coral', alpha=0.5, label=legend2, edgecolor='coral')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    print np.nanmin(xvalues1)
    print np.nanmax(xvalues1)
    xmax = np.nanmax(xvalues1)+0.5
    xmin = np.nanmin(xvalues1)-0.5
    ax1.set_ylim([xmin,xmax])
    ax1.set_xlim([xmin,xmax])
    ymax = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.set_title(plot_title)
    ax1.plot([xmin, xmax], [xmin,xmax], ls="--", c=".3", lw=1.5)
    plt.show()

#######################################################################################################################
## Scatter plot of normalized counts (log10 is the default).                                                         ##
## SampleTuple1 and 2 are tuples that contain the sample name and read type - e.g. ("CJB66D","Normalized to mature") ##
## DataFrames are from normalize_B_to_mature or normalize_AtoB.                                                      ##
## base is the base for log scale (if you want it to be 2, set base=2)                                               ##
#######################################################################################################################

def scatter_plot2(SampleTuple1, SampleTuple2, DataFrame1, DataFrame2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)', base=10, plot_both=True, log_both=True, scaling="auto"):
    print DataFrame1.columns.lexsort_depth
    #print DataFrame1.index.lexsort_depth
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    if plot_both==True and log_both==True:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
        ax1.scatter(xvalues2, yvalues2, c='coral', alpha=0.5, label=legend2, edgecolor='coral')
        ax1.legend(loc=4)
    elif plot_both==1 and log_both==True:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.legend(loc=4)
    elif plot_both==2 and log_both==True:
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.scatter(xvalues2, yvalues2, c='0.3', alpha=1, label=legend2, edgecolor='0.2')
        ax1.legend(loc=4)
    elif plot_both==True and log_both==False:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=False)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=False)
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=False)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=False)
        ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
        ax1.scatter(xvalues2, yvalues2, c='royalblue', alpha=0.5, label=legend2, edgecolor='darkslateblue')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #print np.nanmin(xvalues1)
    #print np.nanmax(xvalues1)
    
    ax1.set_title(plot_title)
    if scaling == "auto":
        if plot_both == True or plot_both == 1:
            xmax = np.nanmax(xvalues1)+0.5
            xmin = np.nanmin(xvalues1)-0.5
        elif plot_both == 2:
            xmax = np.nanmax(xvalues2)+0.5
            xmin = np.nanmin(xvalues2)-0.5
        ax1.set_ylim([xmin,xmax])
        ax1.set_xlim([xmin,xmax])
        ymax = ax1.get_ylim()
        ax1.plot([xmin, xmax], [xmin, xmax], ls="--", c=".3", lw=1.5)
    plt.show()
    return fig1

def scatter_plot3(DataFrame_x, SampleList_x, DataFrame_y, SampleList_y, plot_title='Plot', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)', base=10):
    a=0
    xvalues = []
    while a < len(SampleList_x):
        xvalues.append(get_ratios(DataFrame_x, SampleList_x[a][0], SampleList_x[a][1], log=True, base=base, T0only=True))
        a += 1
    b=0
    yvalues = []
    while b < len(SampleList_y):
        yvalues.append(get_ratios(DataFrame_y, SampleList_y[b][0], SampleList_y[b][1], log=True, base=base, T0only=True))
        b += 1
    print len(xvalues)
    print len(yvalues)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    n = 0
    color_list = ["royalblue","coral","lightsteelblue"]
    while n < len(yvalues):
        if len(xvalues) == len(yvalues):   
            print len(xvalues[n])
            print len(yvalues[n])
            ax1.scatter(xvalues[n], yvalues[n], c=color_list[n], label=legend1, alpha = 0.5, edgecolor=color_list[n])
            n += 1
        elif len(xvalues) == 1:
            print len(xvalues[0])
            print len(yvalues[n])
            ax1.scatter(xvalues[0], yvalues[n], c=color_list[n], label=legend1, alpha = 0.5, edgecolor=color_list[n])
            n += 1
        else:
            print "X and Y variables do not match"
            break
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #print np.nanmin(xvalues[1])
    #print np.nanmax(xvalues[1])
    #xmax = np.nanmax(xvalues[1])+0.5
    #xmin = np.nanmin(xvalues[1])-0.5
    #ax1.set_ylim([xmin,xmax])
    #ax1.set_xlim([xmin,xmax])
    #ymax = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.set_title(plot_title)
    ax1.plot([xmin, xmax], [xmin,xmax], ls="--", c=".3", lw=1.5)
    plt.show()
    return fig1


#######################################################################################################################
## Bar chart of averaged replicates of normalized counts (error bars are the range).                                 ##
## df is DataFrame from normalize_B_to_mature or normalize_AtoB                                                      ##
## sample1 and 2 are sample names (e.g. CJB66D). These should be in level0 of your DataFrame Columns                 ##
## CNAG_lists are the lists of transcripts you want to plot. 1 and 2 will be different colors. You only need to give ##
## CNAG_list1                                                                                                        ##
#######################################################################################################################

def bar_chart(df, read_type, sample1, sample2, CNAG_list1, CNAG_list2=None, plot_title="Stuff", ylabel="Stuff"):
    if CNAG_list2 is not None:
        CNAG_list = CNAG_list1 + CNAG_list2
    else:
        CNAG_list = CNAG_list1
    gene_list = []
    values1 = []
    values2 = []
    #Check if the dataframe index has only transcript names or both transcripts and exons
    if type(df.index[0]) == str:
        for name1, name2 in df:
            if name1 == sample1 and name2 == read_type:
                for gene in CNAG_list:
                    if gene in df.index:
                        values1.append(df[(name1,name2)][gene])
                        gene_list.append(gene)
            elif name1 == sample2 and name2 == read_type:
                for gene in CNAG_list:
                    if gene in df.index:
                        values2.append(df[(name1,name2)][gene])
    else:
        for gene in CNAG_list:
            for transcript in df.iterrows():
                if transcript[0][0] == gene:
                    values1.append(df[(sample1,read_type)][(transcript[0][0], transcript[0][1])])
                    values2.append(df[(sample2,read_type)][(transcript[0][0], transcript[0][1])])
                    gene_list.append(str(transcript[0][0])+"-"+str(transcript[0][1]))
    n=0
    avg_values = []
    errors = []
    while n < len(values1):
        avg_values.append((values1[n]+values2[n])/2)
        errors.append(abs(values1[n]-values2[n])/2)
        n += 1
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    n_groups = len(avg_values)
    index = np.arange(n_groups)
    bar_width = 0.5
    error_config = {'ecolor': '0.3'}
    ax1.bar(index, avg_values, bar_width,alpha=0.9,color='darkslateblue',yerr=errors,error_kw=error_config)
    plt.xticks(index, gene_list, rotation='vertical')
    plt.xlabel('Transcript')
    plt.ylabel(ylabel)
    plt.title(plot_title)
    print ax1.get_children()
    for child in ax1.get_children()[4:4+len(CNAG_list1)]:
        child.set_color('coral')
    plt.show()
    return fig1


#############################################################################
## Get a random set of transcripts from a file with a list of transcripts. ##
## Pretty self explanatory                                                 ##
#############################################################################

def random_transcripts_from_list(file_with_list, number_of_transcripts):
    gene_list = []
    fin = open(file_with_list, "r")
    for line in fin:
        if line.startswith("CNAG"):
            gene = line.split("\t")[0].strip()
            if gene[-2:-1] == "T":
                gene_list.append(gene)
            else:
                gene_list.append(gene+"T0")
    random_list = random.sample(gene_list, number_of_transcripts)
    return random_list

def beeswarm_plot(DataFrame_list, SampleTuple_list, DataFrame2=None, base=10, color_list="blue", select_random=False):
    print SampleTuple_list
    print len(SampleTuple_list)
    values_list = []
    name_list = []
    median_list = []
    a=0
    while a < len(DataFrame_list):
        n=0
        while n < len(SampleTuple_list):
            print SampleTuple_list[n]
            values = get_ratios(DataFrame_list[a], SampleTuple_list[n][0], SampleTuple_list[n][1], log=True, base=base)
            median_list.append(np.median(np.array(values)))
            if select_random != False:
                values = random.sample(values, select_random)
            values_list.append(values)
            name_list.append(SampleTuple_list[n][0])
            n += 1
        a += 1
    print median_list
    bs, ax = beeswarm(values_list, method="swarm", labels=name_list, col=color_list)
    
def histogram(df1, SampleTuple1, df2=None, SampleTuple2=None, xlabel="Intermediate levels", ylabel="Number of introns", bins=100, legend1='All', legend2=None):
    x1 = get_ratios(df1, SampleTuple1[0], SampleTuple1[1], log=True, base=10)
    x1 = [x for x in x1 if str(x) != 'nan']
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.hist(x1, bins=bins, color='royalblue', edgecolor='royalblue', alpha=0.8, label=legend1 )
    if df2 is not None and SampleTuple2 is None:
        x2 = get_ratios(df2, SampleTuple1[0], SampleTuple1[1], log=True, base=10)
        x2 = [x for x in x2 if str(x) != 'nan']
        ax.hist(x2, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    elif SampleTuple2 is not None and df2 is None:
        y1 = get_ratios(df1, SampleTuple2[0], SampleTuple2[1], log=True, base=10)
        y1 = [x for x in y1 if str(x) != 'nan']
        ax.hist(y1, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    elif df2 is not None and SampleTuple2 is not None:
        y2 = get_ratios(df2, SampleTuple2[0], SampleTuple2[1], log=True, base=10)
        y2 = [x for x in y2 if str(x) != 'nan']
        ax.hist(y2, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    if df2 is not None or SampleTuple2 is not None:    
        ax.legend(loc=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    return fig1
    
##Input is the outuput from the Score_splice_sites.py script (Splice_site_strengths.tsv), can also take score_dict directly
def bin_transcripts(input_file, bin_size, bin_by="length", df_list=None):
    score_dict = {}
    if type(input_file) == str:
        with open(input_file, "r") as fin:
            for line in fin:
                if line.startswith("Transcript"): 
                    continue
                else:
                    columns = re.split(r'\t', line)
                    transcript = columns[0]
                    intron = columns[1]
                    five_p_score = float(columns[2])
                    three_p_score = float(columns[3])
                    length = int(columns[4])
                    key = (transcript, int(intron)+1)
                    if bin_by == "length" and length > 30:
                        score_dict[key] = length
                    elif bin_by == "length" and length <= 30:
                        continue
                    elif bin_by == "five_p_score":
                        score_dict[key] = five_p_score
                    elif bin_by == "three_p_score":
                        score_dict[key] = three_p_score
                    else: print "I don't recognize the bin_by value"
    elif type(input_file) == dict:
        score_dict = input_file
    if df_list is not None:
        score_dict = {key: score_dict[key] for key in score_dict if key in df_list }
    sorted_dict = sorted(score_dict.items(), key=operator.itemgetter(1))

    n = 0
    bin1 = []
    values1 = []
    bin2 = []
    values2 = []
    all_scores = []
    for transcript, score in sorted_dict:
        all_scores.append(score)
        if n < bin_size:
            bin1.append(transcript)
            values1.append(score)
            n += 1
        elif n >= bin_size and n < len(sorted_dict)-bin_size:
            n += 1
        else:
            bin2.append(transcript)
            values2.append(score)
            n += 1
    print "Mean intron length or splice site score:"
    print "%0.1f" % np.mean(all_scores)
    print "Mean intron length or splice site score for bottom bin:"
    print "%0.1f" % np.mean(values1)
    print "Mean intron length or splice site score for top bin:"
    print "%0.1f" % np.mean(values2)
    return (bin1, bin2)
        
    
def bin_by_position(df, first_last=True, first_other=False, other_last=False):
    df_list = df.index.tolist()
    intron_dict = {}
    for intron in df_list:
        if intron[0] not in intron_dict:
            intron_dict[intron[0]]=[]
        intron_dict[intron[0]].append(intron[1])
    
    bin1 = []
    bin2 = []
    for transcript, introns in intron_dict.iteritems():
        introns = sorted(introns)
        first = introns[0]
        last = introns[-1]
        if first_last == True:
            bin1.append((transcript, first))
            bin2.append((transcript, last))
        elif first_other == True:
            bin1.append((transcript, first))
            i=1
            for i in range(len(introns)):
                bin2.append((transcript, introns[i]))
        elif other_last == True:
            bin2.append((transcript, last))
            i=0
            for i in range(len(introns)-1):
                bin1.append((transcript, introns[i]))
    if len(bin1) != len(bin2):
        if len(bin1) < len(bin2):
            bin2 = random.sample(bin2, len(bin1))
        elif len(bin1) > len(bin2):
            bin1 = random.sample(bin1, len(bin2))
    return bin1, bin2
    
def cumulative_function(df, score_file, bin_size, SampleTuple1, bin_by=False, bin_tuple=None, SampleTuple2=None, xlabel="Intermediate levels", ylabel="Number of introns", title="CDF", plot_type="CDF", bins=100):
    #df = df[df[SampleTuple1] > 0]
    if bin_by == 'position':
        transcript_list_1, transcript_list_2 = bin_by_position(df, first_last=False, first_other=False, other_last=True)
    elif bin_by != False:
        df_list = df.index.tolist()
        transcript_list_1, transcript_list_2 = bin_transcripts(score_file, bin_size, bin_by, df_list=df_list)
    elif bin_tuple is not None:
        transcript_list_1= bin_tuple[0]
        transcript_list_2= bin_tuple[1]
    new_df1 = pd.DataFrame(columns=df.columns)
    new_df2 = pd.DataFrame(columns=df.columns)

    new_df1 = df[df.index.map(lambda x: x in transcript_list_1)]
    new_df2 = df[df.index.map(lambda x: x in transcript_list_2)]
    
    fig1 = plt.figure(figsize=(8, 6), dpi=600)
    ax1 = fig1.add_subplot(111)
    if plot_type == "CDF":
        x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=False)
        x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=False)
        values1, base1 = np.histogram(x1, bins=bins)
        values2, base2 = np.histogram(x2, bins=bins)
        cumulative1 = np.cumsum(values1)
        base1 = np.insert(base1, 0, 0)
        cumulative1 = np.insert(cumulative1, 0, 0)
        ax1.plot(base1[:-1], cumulative1, c='blue', linewidth=3.0, label="Low 1")
        cumulative2 = np.cumsum(values2)
        base2 = np.insert(base2, 0, 0)
        cumulative2 = np.insert(cumulative2, 0, 0)
        ax1.plot(base2[:-1], cumulative2, c='orangered', linewidth=3.0, label="High 1")
        if SampleTuple2 is not None:
            y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=False)
            values3, base3 = np.histogram(y1, bins=1000)
            y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=False)
            values4, base4 = np.histogram(y2, bins=1000)
            cumulative3 = np.cumsum(values3)
            base3 = np.insert(base3, 0, 0)
            cumulative3 = np.insert(cumulative3, 0, 0)
            ax1.plot(base3[:-1], cumulative3, color='lightskyblue', linewidth=3.0, label="Low 2")
            cumulative4 = np.cumsum(values4)
            base4 = np.insert(base4, 0, 0)
            cumulative4 = np.insert(cumulative4, 0, 0)
            ax1.plot(base4[:-1], cumulative4, color='coral', linewidth=3.0, label="High 2")
            ax1.legend(loc=4)
    
    elif plot_type == "PDF":
        x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=True)
        x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=True)
        x1 = [x for x in x1 if str(x) != 'nan']
        x2 = [x for x in x2 if str(x) != 'nan']
        #x1 = [1e-15 if str(x) == 'nan' else x for x in x1]
        #x2 = [1e-15 if str(x) == 'nan' else x for x in x2]
        print "Zero values removed:"
        print str(bin_size-len(x2))+" from high bin"
        print str(bin_size-len(x1))+" from low bin"
        #ax1.hist(x2, color='coral', edgecolor='coral', label="High 1", bins=bins, alpha=0.9)
        #ax1.hist(x1, color='royalblue', edgecolor='royalblue', label="Low 1", bins=bins, alpha=0.5)
        sns.distplot(x2, color='orangered', label="High", ax=ax1, bins=bins)
        sns.distplot(x1, color='royalblue', label="Low", ax=ax1, bins=bins)
        ax1.legend(loc=1)
        if SampleTuple2 is not None:
            y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=True)
            y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=True)
            y1 = [x for x in y1 if str(x) != 'nan']
            y2 = [x for x in y2 if str(x) != 'nan']
            ax1.hist(y2, color='orangered',  label="High 2", bins=bins)
            ax1.hist(y1, color='lightskyblue', label="Low 2", bins=bins)

    print "KS_statistic, p_value for replicate 1: "
    ks_x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=False)
    ks_x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=False)
    for value in ks_2samp(ks_x1, ks_x2):
        print "%0.1e" % value
    if SampleTuple2 is not None:
        print "KS_statistic, p_value for replicate 2: "
        ks_y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=False)
        ks_y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=False)
        for value in ks_2samp(ks_y1, ks_y2):
            print "%0.1e" % value
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.set_ylim([0.9*len(new_df1),len(new_df1)+len(new_df1)*0.01])
    #ax1.set_ylim([0*len(new_df1),len(new_df1)+len(new_df1)*0.01])
    xmax = np.nanmax([np.nanmax(x1),np.nanmax(x2)])
    ax1.set_xlim([-1,xmax+.05*xmax])
    #ax1.set_xlim([-2, 50])
    ax1.set_title(title)
    plt.show()
    return fig1


##################################################################################
## Functions for normalizing data and creating a CDF that actually adds up to 1 ##
## Takes a list of data                                                         ##
##################################################################################

def normalize(data):
    data = [float(x)/sum(data) for x in data]
    return data

def cdf_values(data, bins='auto'):
    if bins=='auto':
        bins = len(data)/10
    values, base = np.histogram(data, bins=bins)
    values = normalize(values)
    cumulative = np.cumsum(values)
    base = np.insert(base, 0, min(data))
    cumulative = np.insert(cumulative, 0, 0)
    return cumulative, base

def cdf_for_n_lists(list_of_lists, label_list=None, color_list=None, x_title='Lengths', y_title='Fraction of introns'):
    fig = plt.figure(figsize=(8, 6), dpi=600)
    ax = fig.add_subplot(111)
    if label_list is None:
        label_list = range(len(list_of_lists))
    if color_list is None:
        color_list = ['0.3','cornflowerblue','orangered','0.7','limegreen','mediumvioletred']
    
    all_cdfs = []
    n = 0 
    for n in range(len(list_of_lists)):
        cumulative, base = cdf_values(list_of_lists[n])
        ax.plot(base[1:], cumulative, c=color_list[n], linewidth=3.0, label=label_list[n])
        all_cdfs.append(cumulative)
        
    n=0
    for n in range(len(list_of_lists)-1):
        print "Sample "+str(label_list[n])+" vs. sample "+str(label_list[n+1])
        print "%0.1e" % ks_2samp(list_of_lists[n], list_of_lists[n+1])[1]
        if n == 1:
            print "Sample "+str(label_list[n-1])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(all_cdfs[n-1], list_of_lists[n+1])[1]
        elif n == 2:
            print "Sample "+str(label_list[n-1])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(list_of_lists[n-1], list_of_lists[n+1])[1]
            print "Sample "+str(label_list[n-2])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(list_of_lists[n-2], list_of_lists[n+1])[1]
        elif n == 3:
            print "Sample "+str(label_list[n-1])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(list_of_lists[n-1], list_of_lists[n+1])[1]
            print "Sample "+str(label_list[n-2])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(list_of_lists[n-2], list_of_lists[n+1])[1]
            print "Sample "+str(label_list[n-2])+" vs. sample "+str(label_list[n])
            print "%0.1e" % ks_2samp(list_of_lists[n-2], list_of_lists[n])[1]
            print "Sample "+str(label_list[n-3])+" vs. sample "+str(label_list[n+1])
            print "%0.1e" % ks_2samp(list_of_lists[n-3], list_of_lists[n+1])[1]
            print "Sample "+str(label_list[n-3])+" vs. sample "+str(label_list[n])
            print "%0.1e" % ks_2samp(list_of_lists[n-3], list_of_lists[n])[1]
    
    #ax.legend(bbox_to_anchor=(1, 0), loc='lower left', fontsize=14)
    ax.legend(fontsize=14)
    plt.ylabel(y_title, fontsize=16)
    plt.xlabel(x_title, fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.show()
    
    return fig
        
def sns_distplot_n(list_of_lists, label_list=None, color_list=None, x_title='Lengths', y_title='Franction of introns'):
    fig = plt.figure(figsize=(8, 6), dpi=600)
    ax = fig.add_subplot(111)
    if label_list is None:
        label_list = range(len(list_of_lists))
    if color_list is None:
        color_list = ['0.3','cornflowerblue','orangered','0.7','limegreen','mediumvioletred']
    
    n = 0 
    for n in range(len(list_of_lists)):
        ax = sns.distplot(list_of_lists[n])

        
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


def set_cdf_xlim(df, metric, x_var):
    maxs = []
    mins = []
    for quart in range(0,4):
        data = df[df['quartile'] == quart][metric+' avg']
        try:
            if max(data)-np.percentile(data, 99) < 0.05*max(data):
                maxs.append(max(data)+ 0.2*max(data))
            else:
                maxs.append(np.percentile(data, 99))
            mins.append(np.percentile(data, 1))
        except ValueError:
            pass
    return (min(mins), max(maxs))

def intron_char_plots(df, metric, color="0.8"):
    df['alt splicing'] = df['alt splicing'].apply(bool)
    gridkw = dict(height_ratios=[1,1])
    pal = sns.set_palette(sns.light_palette(color, reverse=True))
    fig, (ax1, ax2) = plt.subplots(2, 4, figsize=(12,8), gridspec_kw=gridkw)
    
    for n, factor in enumerate(['alt splicing', 'intron position','Base 5-4', 'Base 5-5']):
        if 'Base' in factor:
            df = df.sort_values(factor)
        sns.boxplot(x=factor, y=metric+' avg', data=df, ax=ax1[n], palette=pal)
        plt.setp(ax1[n].lines, color=".1")
        ax1[n].tick_params(axis='both', which='major', labelsize=14)
        ax1[n].set_xlabel(factor.capitalize(), fontsize=16)
        ax1[n].set_ylabel('log2 '+metric, fontsize=16)

    quart_info = {0:('Bottom','royalblue'),1:('Lower','skyblue'),2:('Upper','orange'),3:('Top','orangered')}
    for n, x_var in enumerate(['intron size','introns in transcript','W Intron Retention avg','5p score']):
        # Separate explanatory variable into quartiles
        df['quartile'] = 0
        for m, quart in enumerate([25,50,75]):
            df.loc[df[df[x_var] >= np.percentile(df[x_var], quart)].index, 'quartile'] = m+1

        for quart in range(0,4):
            data = df[df['quartile'] == quart][metric+' avg']
            try:
                sns.kdeplot(data, ax=ax2[n], bw=2, cumulative=True, linewidth=3, 
                        color=quart_info[quart][1])
            except ValueError:
                pass
        ax2[n].set_xlim(set_cdf_xlim(df, metric, x_var))
        ax2[n].set_xlabel('log2 '+metric, fontsize=16)
        ax2[n].set_ylabel('Fraction of intermediates', fontsize=16)
        ax2[n].set_title(x_var.capitalize(), fontsize=16)
        ax2[n].tick_params(axis='both', which='major', labelsize=14)
        ax2[n].legend_.remove()

    fig.tight_layout()
    plt.show()
    return fig        

def fig_fix_font(fig, xlabel_list=None, ylabel_list=None, ylabel=None):
    ax_list = fig.axes
    for n, ax in enumerate(ax_list):
        if xlabel_list is None:
            xlabel = ax.get_xlabel()
        else:
            xlabel = xlabel_list[n]
        if ylabel_list is None:
            if ylabel is None:
                ylabel = ax.get_ylabel()
            else:
                ylabel = ylabel
        else:
            ylabel = ylabel_list[n]
            
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
    return fig
        
def mut_compare_boxplots(df, metrics=['Intermediate Level', 'Precursor', 'W Intron Retention'], color='0.8'):
    pal = sns.set_palette(sns.light_palette(color))
    x_list = []
    for metric in metrics:
        fig, ax = plt.subplots(ncols=4, nrows=3, figsize=(10,10))
        
        for n, prop in enumerate(['alt splicing','intron position']):
            sns.boxplot(x=('Peaks',prop), y=('All',metric+' log2 ratio1'), data=df, ax=ax[0][n], palette=pal)
            x_list.append(prop.capitalize())

        for n, pos in enumerate(['Base 3-2','Base 3-3']):
            x_list.append(pos)
            sns.boxplot(x=('Peaks',pos), y=('All',metric+' log2 ratio1'), data=df.sort_values(('Peaks',pos)), ax=ax[0][n+2], palette=pal)

        for n, pos in enumerate([x[1] for x in df.columns if ('Base 5-' in x[1]) and (x[1][-1] in ['4','5','6','7'])]):
            x_list.append(pos)
            sns.boxplot(x=('Peaks',pos), y=('All',metric+' log2 ratio1'), data=df.sort_values(('Peaks',pos)), ax=ax[n/4+1][n%4], palette=pal)

        for n, pos in enumerate([x[1] for x in df.columns if ('branch-' in x[1]) and (x[1][-1] != '3')]):
            x_list.append(pos)
            sns.boxplot(x=('Peaks',pos), y=('All',metric+' log2 ratio1'), data=df.sort_values(('Peaks',pos)), ax=ax[2][n], palette=pal)
        
        fig = fig_fix_font(fig, xlabel_list=x_list, ylabel=metric)
        plt.suptitle(metric, y=1.02, fontsize=16)
        fig.tight_layout()
        plt.show()
        plt.clf()
        
def set_cdf_xlim_Z(df, metric, x_var):
    maxs = []
    mins = []
    for change in ['Other','Up','Down']:
        data = df[df[('All',metric+' change')] == change][('Peaks',x_var)]
        if len(data) > 0:
            if max(data)-np.percentile(data, 99) < 0.05*max(data):
                maxs.append(max(data)+ 0.2*max(data))
            else:
                maxs.append(np.percentile(data, 99))
            mins.append(np.percentile(data, 1))
    return (min(mins), max(maxs))
        
def mut_compare_cdfs(df, metrics=['Intermediate Level', 'Precursor', 'W Intron Retention']):
    groups = {'Other':'0.6','Up':'orangered','Down':'royalblue'}
    for metric in metrics:
        x_list = []
        fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(10,8))

        for n, prop in enumerate(['intron size','5p score','branch score','3p score', 'introns in transcript','transcript size']):
            prop_data = df[('Peaks',prop)].replace([np.inf, np.inf*-1], np.NaN).dropna()
            x_list.append(prop.capitalize())
            for change in ['Other','Up','Down']:
                plot_df = df[df[('All', metric+' change')] == change]
                if len(plot_df) > 10:
                    sns.kdeplot(plot_df[('Peaks',prop)], ax=ax[n/3][n%3], bw=2, cumulative=True, linewidth=3, 
                                    color=groups[change], label=change)
                    ax[n/3][n%3].set_xlim(set_cdf_xlim_Z(df, metric, prop))
            ax[n/3][n%3].legend_.remove()
        
        fig = fig_fix_font(fig, xlabel_list=x_list, ylabel='Fraction of introns')
        fig.tight_layout()
        plt.suptitle(metric, y=1.02, fontsize=16)
        plt.show()
        plt.clf()
        
##################################################################################
## Tools for plotting IGV like traces based on dataframes and bam files         ##
##################################################################################
def config_quant(config_file):
    bam_dict = {}
    with open(config_file,'r') as config:
        for line in config:
            #bam files for quantitation should be file,quant,A1
            if 'quant' in line:
                bam_dict[line.split(',')[-1].strip()] = line.split(',')[0]
    return bam_dict

def gene_patches3(tx, tx_dict, ax, arrow=False):
    iso_list = [x for x in tx_dict if tx in x]
    if len(iso_list) == 0:
        return None
    
    for n, iso in enumerate(iso_list):
        start, end, strand, CDS_start, CDS_end, exons, chrom = SP.tx_info(iso, tx_dict)

        if arrow is False:
            tx_patch = patches.Rectangle((start,0.8-n*0.15),end-start,0.04,edgecolor='0.1',facecolor='0.1')
            ax.add_patch(tx_patch)
        else:
            if strand == '+':
                ax.arrow(start, 0.9, end-start-0.02*(end-start), 0, linewidth=2, head_width=0.1, 
                         head_length=0.02*(end-start), fc='k', ec='k')
            elif strand == '-':
                ax.arrow(end, 0.9, start-end-0.02*(start-end), 0, linewidth=2, head_width=0.1, 
                         head_length=0.02*(end-start), fc='k', ec='k')

        if exons is not None:
            exon_patches = []
            for exon_start, exon_stop in exons:
                exon_patches.append(patches.Rectangle((exon_start, 0.775-n*0.15), exon_stop-exon_start, 0.10,
                                                      edgecolor='0.1',facecolor='0.1'))
            for patch in exon_patches:
                ax.add_patch(patch)
        else:
            CDS_patch = patches.Rectangle((CDS_start, 0.75-n*0.15),CDS_end-CDS_start, 0.10, edgecolor='0.1', facecolor='0.1')
            ax.add_patch(CDS_patch)
        ax.get_yaxis().set_ticks([])
    return strand   

def generate_read_series_A(bam_iterator, chrom, start, end, strand, baseline=0):
    s = pd.Series(baseline, index=range(start, end))
    for read in bam_iterator:
        if read.is_reverse and strand == '+':
            pos = read.reference_end
            if pos not in s.index:
                s[pos] = 0
            s[pos] += 1
        elif not read.is_reverse and strand == '-':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] += 1
    s = s.sort_index()
    return s

def generate_read_series_B(bam_iterator, chrom, start, end, strand, baseline=0):
    s = pd.Series(baseline, index=range(start, end))
    for read in bam_iterator:
        # Check to make sure the read doesn't span an intron
        intron = False
        for cigarop in read.cigar:
            if(cigarop[0]==3):
                intron = True
                break
        # Get reads (flip if read2)
        if intron is False:
            pos_range = range(read.reference_start, read.reference_end)
            if read.is_reverse and strand == '+':
                s[s.index.isin(pos_range)] += 1
            elif not read.is_reverse and strand == '-':
                s[s.index.isin(pos_range)] += 1
    s = s.sort_index()
    return s

def generate_read_series_PE(bam_iterator, chrom, start, end, strand, baseline=0):
    s = pd.Series(baseline, index=range(start, end))
    for read in bam_iterator:
        # Check to make sure the read doesn't span an intron
        intron = False
        for cigarop in read.cigar:
            if(cigarop[0]==3):
                intron = True
                break
        # Get reads (flip if read2)
        if intron is False:
            pos_range = range(read.reference_start, read.reference_end)
            if read.is_read1:
                if read.is_reverse and strand == '+':
                    s[s.index.isin(pos_range)] += 1
                elif not read.is_reverse and strand == '-':
                    s[s.index.isin(pos_range)] += 1
            elif read.is_read2:
                if not read.is_reverse and strand == '+':
                    s[s.index.isin(pos_range)] += 1
                elif read.is_reverse and strand == '-':
                    s[s.index.isin(pos_range)] += 1
    s = s.sort_index()
    return s

def get_introns(open_bam, chrom, start, end, strand):
    introns = {}
    if strand == '+':
        introns.update(open_bam.find_introns((read for read in open_bam.fetch(chrom, start, end) 
                                              if read.is_reverse and r.is_read1)))
        introns.update(open_bam.find_introns((read for read in open_bam.fetch(chrom, start, end) 
                                              if not read.is_reverse and r.is_read2)))
    if strand == '-':
        introns.update(open_bam.find_introns((read for read in open_bam.fetch(chrom, start, end) 
                                              if not read.is_reverse and r.is_read1)))
        introns.update(open_bam.find_introns((read for read in open_bam.fetch(chrom, start, end) 
                                              if read.is_reverse and r.is_read2)))
    print introns

def igv_plots_quant(config_file, df, organism, save_dir=None):
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    bam_files = config_quant(config_file)
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    if save_dir is None:
        save_dir = config_file.split('.txt')[0]+'_plots/'
    elif save_dir[-1] != '/':
        save_dir = save_dir+'/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    open_bams = {}
    totals = {}
    for rep, bam in bam_files.iteritems():
        open_bams[rep] = pysam.Samfile(bam)
        print bam
        total = check_output(['samtools','view','-F 0x04','-c',bam]).strip()
        total = float(total)/1000000.
        print total
        totals[rep] = total

    order = ['A1','B1','W1','A2','B2','W2']
    colors = ['0.1','royalblue','0.6','0.1','royalblue','0.6']
    
    # Get read series for all
    for ix, r in df.iterrows():
        fig, ax = plt.subplots(4, 2, figsize=(12,6), sharex=True)
        fig.subplots_adjust(hspace=0)

        if r['strand'] == '+':
            start = r['position']-100
            end = r['position'] + r['intron size']+100
            
        elif r['strand'] == '-':
            start = r['position']-r['intron size']-100
            end = r['position']+100
        
        max_y = [0,0,0]
        for n, rep in enumerate(order):
            bam_iter = open_bams[rep].fetch(r['chromosome'], int(start), int(end))
            if rep.startswith('A'):
                s = generate_read_series_A(bam_iter, r['chromosome'], int(start), int(end), r['strand'])
                s = s.divide(totals[rep])
                peak_range = range(int(r['position'])-2, int(r['position'])+2)
                max_y[0] = max([max_y[0],sum(s[s.index.isin(peak_range)])])
                
            else:
                s = generate_read_series_B(bam_iter, r['chromosome'], int(start), int(end), r['strand'])
                s = s.divide(totals[rep])
                if rep.startswith('B'):
                    max_y[1] = max([max_y[1], max(s)])
                if rep.startswith('W'):
                    max_y[2] = max([max_y[2], max(s)])
            
            ax[n%3][n/3].bar(s.index, s, linewidth=1, color=colors[n], edgecolor=colors[n], zorder=2)
            ax[n%3][n/3].tick_params(axis='both', which='major', labelsize=14)

        strand = gene_patches3(r['transcript'], tx_dict, ax[3][0], arrow=False)
        strand = gene_patches3(r['transcript'], tx_dict, ax[3][1], arrow=False)
        ax[3][0].set_xlim(start,end)
        ax[3][0].set_xlim(start,end)
        if strand == '-':
            ax[3][0].invert_xaxis()
            ax[3][1].invert_xaxis()
        
        for n in range(6):
            ax[n%3][n/3].set_ylim(0,max_y[n%3]+0.1*max_y[n%3])
            ax[n%3][n/3].set_xlim(start, end)
            if strand == '-':
                ax[n%3][n/3].invert_xaxis()

        ax[0][0].set_ylabel('RPM', fontsize=16)
        ax[0][1].set_ylabel('RPM', fontsize=16)
        ax[0][0].set_title(r['transcript'], fontsize=16)
        ax[0][1].set_title(r['transcript'], fontsize=16)
        ax[3][0].get_xaxis().set_ticks([])
        ax[3][1].get_xaxis().set_ticks([])
        plt.show()
        
        fig.savefig(save_dir+r['transcript']+'_ABW.pdf', format='pdf')
        plt.clf()

def igv_plots_mut(config_file, df, organism, mut_color='limegreen', save_dir=None):
    bam_dict = SP.config_mut_quant(config_file)
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    if save_dir is None:
        save_dir = config_file.split('.txt')[0]+'_plots/'
    elif save_dir[-1] != '/':
        save_dir = save_dir+'/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    open_bams = {}
    geno_list = [None, None, None, None]
    total_list = [0,0,0,0]
    for genotype, bams in bam_dict.iteritems():
        A_reps = {k:v for k,v in bams.items() if k[0] == 'A'}
        for rep, bam in A_reps.iteritems():
            name = bam.split('/')[-1].split('_sorted')[0]
            open_bams[name] = pysam.Samfile(bam)
            
            print name
            total = check_output(['samtools','view','-F 0x04','-c',bam]).strip()
            total = float(total)/1000000.
            print total
            
            if genotype == 'WT':
                if rep == 'A1':
                    ix = 0
                else:
                    ix = 1
            else:
                if rep == 'A1':
                    ix = 2
                else:
                    ix = 3
                    
            geno_list[ix] = name
            total_list[ix] = total

    #total_list = [x/max(total_list) for x in total_list]
    colors = ['0.1','0.1',mut_color,mut_color]
    
    # Get read series for all
    for ix, r in df.iterrows():
        fig, ax = plt.subplots(5, figsize=(12,10), sharex=True)
        fig.subplots_adjust(hspace=0)

        if r[('Peaks','strand')] == '+':
            start = r[('Peaks','position')]-100
            end = r[('Peaks','position')] + r[('Peaks','intron size')]+100
            
        elif r[('Peaks','strand')] == '-':
            start = r[('Peaks','position')]-r[('Peaks','intron size')]-100
            end = r[('Peaks','position')]+100
        
        max_y = 0
        for n, bam in enumerate(geno_list):
            bam_iter = open_bams[bam].fetch(r[('Peaks','chromosome')], int(start), int(end))
            s = generate_read_series_A(bam_iter, r[('Peaks','chromosome')], int(start), int(end), r[('Peaks','strand')])
            s = s.divide(total_list[n])
            ax[n].bar(s.index, s, linewidth=1, color=colors[n], edgecolor=colors[n], zorder=2)
            ax[n].tick_params(axis='both', which='major', labelsize=14)
            
            peak_range = range(int(r[('Peaks','position')])-2, int(r[('Peaks','position')])+2)
            max_y = max([max_y,sum(s[s.index.isin(peak_range)])])
            
        
        strand = gene_patches3(r[('Peaks','transcript')], tx_dict, ax[4], arrow=False)
        ax[4].set_xlim(start,end)
        if strand == '-':
            ax[4].invert_xaxis()
        
        for n in range(4):
            ax[n].set_ylim(0,max_y+0.1*max_y)
            #ax[n].plot([r[('Peaks','position')],r[('Peaks','position')]], [0,max_y+0.1*max_y],'--',color='0.9', zorder=1)
            ax[n].set_xlim(start, end)
            if strand == '-':
                ax[n].invert_xaxis()

        ax[0].set_ylabel('RPM', fontsize=16)
        ax[0].set_title(r[('Peaks','transcript')], fontsize=16)
        ax[0].get_xaxis().set_ticks([])
        plt.show()
        fig.savefig(save_dir+r[('Peaks','transcript')]+'_A.pdf', format='pdf')
        plt.clf()