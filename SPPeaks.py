__author__ = 'jordanburke'


''' These are tools to interpret and filter the output of the ChangePoints peak picking method developed by Jessica Li and Daria Merkurjev.
'''

import sys
import pandas as pd
import random
import re
import numpy as np
from scipy import stats
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import SPTools as SP
from math import log
from matplotlib import pyplot as plt
import pysam

############################################################################################################
## Tools for recovering sequences in the region around peaks. Peak picking is not super accurate and has  ##
## off by 1-2 errors. Currently looks in range around peak for the splice site                            ##
############################################################################################################

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N', 'Y':'R', 'R':'Y'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])
                
#####################################################################################
## Trim bed file reads to just the 5' ends. Helpful for picking peaks in A samples ##
#####################################################################################

def convert_bed_file(bed_file):
    counter = 0
    with open(bed_file, "r") as fin:
        for line in fin:
            counter += 1
            columns = re.split(r'\t', line)
            new_list = []
            new_list.append(columns[0])
            if columns[5].strip() == "+":
                new_list.append(columns[1])
                new_list.append(str(int(columns[1])+1))
            elif columns[5].strip() == "-":
                new_list.append(str(int(columns[2])-1))
                new_list.append(columns[2])
            else:
                print "unknown strand"
                print counter
            new_list = new_list+columns[3:6]
            with open("{0}_5p.bed".format(bed_file.split("/")[-1].split(".")[0]), "a") as fout:
                new_line = "\t".join(new_list)
                fout.write(new_line)

#####################################################################################################
## New simplified code for processing peaks from Jessica Li                                        ##
#####################################################################################################

# Make a bedgraph file from ChangePoints output
def convert_CP_output(cp_output, format2=False):
    line_list = []
    with open(cp_output, 'r') as f:
        for line in f:
            peak = line.split('\t')
            if len(peak[0]) == 0:
                peak = peak[1:]
            if peak[0].startswith('Chr_'):
                chrom = 'chr'+peak[0][4:]
            else:
                chrom = peak[0].strip()
            a = int(peak[1])
            b = int(peak[1])-1
            if format2 is True: height = peak[3]
            else: height = peak[2]
            line_list.append([chrom, str(a), str(b), height])
    name = cp_output.split('/')[-1].split('.')[0]
    with open(name+'.bedgraph', 'w') as fout:
        for line in line_list:
            fout.write('\t'.join(line)+'\n')

#Function to read peak output file
def CP_peaks_by_gene(fin, transcript_dict, cutoff=5):
    genes_by_chr = {}
    for tx, info in transcript_dict.iteritems():
        if info[3] not in genes_by_chr:
            genes_by_chr[info[3]] = []
        genes_by_chr[info[3]].append(tx)
    
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3'}
    
    peak_count = 0
    line_count = 0
    peaks_by_gene = {}
    strand_dict = {'0':'-','1':'+'}
    with open(fin,'r') as f:
        for line in f:
            data = line.split('\t')
            if len(data[0]) > 0 and len(data) > 3:
                chrom = data[0]
                peak = int(data[1])
                strand = strand_dict[data[2]]
                peak_height = float(data[3])
            elif len(data[0]) > 0 and len(data) == 3:
                chrom = data[0]
                if '_' in chrom:
                    chrom = 'chr'+chrom.split('_')[-1]
                peak = int(data[1])
                peak_height = float(data[2])
                strand = None
            else:
                chrom = data[1]
                peak = int(data[2])
                strand = strand_dict[data[3]]
                peak_height = float(data[4])

            if chrom in rom_lat: chrom = rom_lat[chrom]
            
            if peak_height >= cutoff:
                line_count += 1
                if chrom in genes_by_chr:
                    tx_list = genes_by_chr[chrom]
                    for tx in tx_list:
                        start = transcript_dict[tx][0]
                        end = transcript_dict[tx][1]
                        if peak > start and peak < end and (strand == transcript_dict[tx][2] or strand is None):
                            peak_count += 1
                            if tx not in peaks_by_gene:
                                peaks_by_gene[tx] = []
                            peaks_by_gene[tx].append([peak,peak_height,strand])
                            
    #print peak_count
    return peaks_by_gene
 
#Function to compare untagged and 2 replicates and pick peaks that are in both but not in untagged
def CP_compare_reps(untagged, tagged1, tagged2):
    tx_list = list(set(tagged1.keys()).intersection(tagged2.keys()))
    new_peak_dict = {}
    peak_count = 0
    for tx in tx_list:
        new_peak_dict[tx] = []
        for peak,peak_height,strand in tagged1[tx]:
            if peak in zip(*tagged2[tx])[0]:
                if tx not in untagged or peak not in zip(*untagged[tx])[0]:
                    new_peak_dict[tx].append([peak,peak_height,strand])
                    peak_count += 1
    print peak_count
    return new_peak_dict

#Function to check peaks against annotation
def CP_compare_to_annotation(peaks, ss_dict, transcript_dict):
    five_count = 0
    three_count = 0
    other_count = 0
    intronic_count = 0
    peak_count = 0
    for tx, peak_list in peaks.iteritems():
        peak_count += len(peak_list)   
    print peak_count
    compare_df = pd.DataFrame(index = range(peak_count+1), columns=['transcript','chromosome','strand','position','height','type'])
    n=0
    for tx, info in ss_dict.iteritems():
        chrom = transcript_dict[tx][3]
        if tx in peaks:
            if len(info[0]) > 0 and len(peaks[tx]) > 0:
                for peak, height, strand in peaks[tx]:
                    if strand is None:
                        strand = transcript_dict[tx][2]
                    peak_range = range(peak-3,peak+3)
                    try:
                        annotated = False
                        for pos in peak_range:
                            if pos in info[0]:
                                five_count += 1
                                compare_df.ix[n] = [tx[:-2], chrom, strand, int(pos), height, "5prime"]
                                annotated = True
                                break

                            elif pos in info[1]:
                                three_count += 1
                                compare_df.ix[n] = [tx[:-2], chrom, strand, int(pos), height, "3prime"]
                                annotated = True
                                break
                        if annotated is False:
                            other_count += 1
                            intron_flag = False
                            m=0
                            for m in range(len(info[0])):
                                if peak > info[0][m] and peak < info[1][m]:
                                    intronic_count += 1
                                    intron_flag = True
                                    break
                            if intron_flag is True:
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "intronic"]
                            else:
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "other"]
                    except IndexError:
                        print tx
                        print n
                        print peak_count
                        print len(info)
                    n+=1
    print "5prime annotated sites: "+str(five_count)
    print "3prime annotated sites: "+str(three_count)
    print "Unpredicted peaks: "+str(other_count)
    print "Unpredicted peaks in introns: "+str(intronic_count)
    compare_df.dropna(how='all',inplace=True)
    return compare_df


def collapse_unpredicted_peaks(df):
    tx_list = list(set(df['transcript'].tolist()))
    for tx in tx_list:
        tx_df = df[df['transcript'] == tx]
        tx_df = tx_df[tx_df['type'].isin(['other','intronic'])]
        index=tx_df.index
        n=0
        for n in range(len(index)):
            if tx_df['height'][index[n]] < 10:
                if index[n] in df.index:
                    df.drop(index[n], inplace=True)
            else:
                m=0
                for m in range(len(index)):
                    spacing = abs(tx_df['position'][index[n]]-tx_df['position'][index[m]])
                    if spacing > 0 and spacing <= 2:
                        if tx_df['height'][index[n]] > tx_df['height'][index[m]]:
                            if index[m] in df.index:
                                df.drop(index[m], inplace=True)
                        if tx_df['height'][index[n]] < tx_df['height'][index[m]]:
                            if index[n] in df.index:
                                df.drop(index[n], inplace=True)
                            break
                        else:
                            continue
    #print "Number of unpredicted peaks after condensing:"
    #print len(df[df['type'].isin(['other','intronic'])])
    #print "Number of intronic peaks after condensing:"
    #print len(df[df['type'] == 'intronic'])
    df.reset_index(inplace=True)
    return df

#Add sequences and check whether they're splice sites
def add_sequence_to_df(df, fa_dict, flag=False):
    seq_df = df
    sequence = []
    looks_like = []
    for index, row in df.iterrows():
        if row['strand'] == '+':
            if flag is True:
                seq = fa_dict[row['chromosome']][row['position']-6:row['position']+6]
            else:
                seq = fa_dict[row['chromosome']][row['position']-5:row['position']+7]
        elif row['strand'] == '-':
            if flag is True:
                seq = fa_dict[row['chromosome']][row['position']-6:row['position']+6]
            else:
                seq = fa_dict[row['chromosome']][row['position']-6:row['position']+6]
            seq = reverse_complement(seq)
        sequence.append(seq)
        
        if row['type'] == '5prime':
            looks_like.append('5prime')
        elif row['type'] == '3prime':
            looks_like.append('3prime')
        else:
            if seq[6:8] == 'GT' or seq[6:8] == 'GC' or seq[6:8] == 'AT':
                looks_like.append(seq[6:8])
            elif seq[4:6] == 'AG' or seq[4:6] == 'AC':
                looks_like.append(seq[4:6])
            else:
                looks_like.append('')
    seq_df['sequence'] = sequence
    seq_df['looks like'] = looks_like
    return seq_df

def peak_to_seq_pipeline(untagged_peak_file, tagged1_peak_file, tagged2_peak_file, gff3, fasta, junction_df=None, branch_df=None, cutoff=5, name='CP_peaks'):
    
    if 'pombe' in gff3: organism = 'pombe'
    else: organism = None
        
    transcript_dict = SP.build_transcript_dict(gff3, organism=organism)
    print "Finding peaks in transcripts..."
    
    print untagged_peak_file
    untagged = CP_peaks_by_gene(untagged_peak_file, transcript_dict, cutoff=cutoff)
    
    print tagged1_peak_file
    tagged1 = CP_peaks_by_gene(tagged1_peak_file, transcript_dict, cutoff=cutoff)
    
    print tagged2_peak_file
    tagged2 = CP_peaks_by_gene(tagged2_peak_file, transcript_dict, cutoff=cutoff)
    
    print "Comparing peaks between replicates..."
    peaks = CP_compare_reps(untagged, tagged1, tagged2)
    
    print "Checking peaks against annotation..."
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    peak_df = CP_compare_to_annotation(peaks, ss_dict, transcript_dict)
    peak_df = collapse_unpredicted_peaks(peak_df)
    peak_df['genome coord'] = peak_df['chromosome'].str.cat(peak_df['position'].apply(int).apply(str), sep=':')
    
    if type(fasta) == str:
        fasta = SP.make_fasta_dict(fasta)
    print "Adding sequences..."
    peak_seq_df = add_sequence_to_df(peak_df, fasta, flag=flag)
    
    print "Writing bedgraph..."
    with open(name+'.bedgraph', 'w') as fout:
        for ix, r in peak_seq_df.iterrows():
            if r['strand'] == '+':
                position2 = r['position']+1
                height = r['height']
            elif r['strand'] == '-':
                position2 = r['position']-1
                height = r['height']*-1
            line_list = [r['chromosome'], r['position'], position2, height, '\n']
            line_list = [str(x) for x in line_list]
            line = '\t'.join(line_list)
            fout.write(line)
    
    print "Completed"
    return peak_seq_df


### ****Cross-validation**** ###
def compare_peak_junc_df(peak_df, junc_df, organism = None):
    print str(len(peak_df))+' peaks'
    print str(len(junc_df))+' junctions'
    
    new_df = pd.DataFrame(columns=peak_df.columns)

    junc_type = []
    ann_seq1 = []
    ann_seq2 = []
    seq1 = []
    seq2 = []
    junc_size = []
    ann_size = []
    junc_coords = []
    ann_coords = []
    match_count = 0
    for tx in list(set(peak_df['transcript'].tolist())):
        tx_peak = peak_df[peak_df['transcript'] == tx]
        if organism == 'pombe':
            tx_junc = junc_df[junc_df['transcript'] == tx+'.1']
        else:
            tx_junc = junc_df[junc_df['transcript'] == tx+'T0']
        for index, row in tx_peak.iterrows():
            match_flag = False

            peak_range = range(row['position']-1,row['position']+1)
            for pos in peak_range:
                if match_flag is False:
                    tx_juncA = tx_junc[tx_junc['type'] != '5p tethered']
                    tx_juncB = tx_junc[tx_junc['type'] != '3p tethered']
                    
                    if pos in tx_juncA['start'].tolist():
                        match_count += 1
                        match_flag = True
                        for junc_index, junc_row in tx_juncA[tx_juncA['start'] == pos].iterrows():
                            junc_type.append(junc_row['type'])
                            ann_seq1.append(junc_row['annotated sequence1'])
                            ann_seq2.append(junc_row['annotated sequence2'])
                            seq1.append(junc_row['sequence1'])
                            seq2.append(junc_row['sequence2'])
                            junc_size.append(junc_row['size'])
                            ann_size.append(junc_row['annotated intron size'])
                            junc_coords.append((junc_row['start'],junc_row['end']))
                            ann_coords.append((junc_row['annotated intron start'],junc_row['annotated intron end']))
                            new_df = new_df.append(row)
                        break
        
                    elif pos in tx_juncB['end'].tolist():
                        match_count += 1
                        match_flag = True
                        for junc_index, junc_row in tx_juncB[tx_juncB['end'] == pos].iterrows():
                            junc_type.append(junc_row['type'])
                            ann_seq1.append(junc_row['annotated sequence1'])
                            ann_seq2.append(junc_row['annotated sequence2'])
                            seq1.append(junc_row['sequence1'])
                            seq2.append(junc_row['sequence2'])
                            junc_size.append(junc_row['size'])
                            ann_size.append(junc_row['annotated intron size'])
                            junc_coords.append((junc_row['start'],junc_row['end']))
                            ann_coords.append((junc_row['annotated intron start'],junc_row['annotated intron end']))
                            new_df = new_df.append(row)
                        break                   
    print "Overlap:"
    print match_count
    
    new_df['junction type'] = junc_type
    new_df['junction sequence1'] = seq1
    new_df['junction sequence2'] = seq2
    new_df['annotated sequence1'] = ann_seq1
    new_df['annotated sequence2'] = ann_seq2
    new_df['junction size'] = junc_size
    new_df['annotated intron size'] = ann_size
    new_df['junction coords'] = junc_coords
    new_df['annotated intron coords'] = ann_coords
    return new_df


### Function to determine how enriched dinucleotide before and after peak is compared to genome
def peak_seq_enrichment(df, organism):
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    nuc_prob = SP.gc_content(fa_dict)
    p_dict = {'A':nuc_prob[0], 'T':nuc_prob[2], 'C':nuc_prob[1], 'G':nuc_prob[3]}
    
    unpeaks = df[df['type'] == 'other']
    unpeaks = unpeaks.append(df[df['type'] == 'intronic'])
    print "Number of unpredicted peaks:"
    print len(unpeaks)
    nucs = ['G','A','C','T']
    dinucs = set()
    for nuc in nucs:
        for nuc2 in nucs:
            dinucs.add(nuc+nuc2)
    
    five = {}
    three = {}
    for dinuc in dinucs:
        five[dinuc] = len(unpeaks[unpeaks['sequence'].str[6:8].str.contains(dinuc)])
        three[dinuc] = len(unpeaks[unpeaks['sequence'].str[4:6].str.contains(dinuc)])

    five_LO = {}
    three_LO = {}
    for dinuc in five.keys():
        p_dinuc = p_dict[dinuc[0]]*p_dict[dinuc[1]]
        phat_dinuc = five[dinuc]/float(len(unpeaks))
        phat_dinuc2 = three[dinuc]/float(len(unpeaks))

        SE = np.sqrt(phat_dinuc*(1-phat_dinuc)/len(unpeaks))
        SE2 = np.sqrt(phat_dinuc2*(1-phat_dinuc2)/len(unpeaks))
        Z = (phat_dinuc-p_dinuc)/SE
        Z2 = (phat_dinuc2-p_dinuc)/SE2

        pvalue = stats.norm.sf(Z)
        pvalue2 = stats.norm.sf(Z2)
        LO = np.log((1-pvalue)/pvalue)
        LO2 = np.log((1-pvalue2)/pvalue2)

        five_LO[dinuc] = LO
        three_LO[dinuc] = LO2

    fig, ax = plt.subplots(figsize=(12,6))
    width = 0.35
    ind = np.arange(len(five_LO.keys()))
    rects2 = ax.bar(ind, three_LO.values(), width, color='crimson', edgecolor='crimson', label='Before peak')
    rects1 = ax.bar(ind + width, five_LO.values(), width, color='indigo', edgecolor='indigo', label='After peak')
    ax.plot([-1,17],[0,0],'-', color='black')
    ax.plot([-1,17],[2.94,2.94], '--', color='0.7', label='95% CI')
    ax.plot([-1,17],[-2.94,-2.94], '--', color='0.7')

    ax.set_xlim([-1,17])
    ax.set_xticklabels(five_LO.keys(), fontsize=12)
    ax.set_xticks(ind + width / 2)
    ax.set_ylabel('Log odds dinucleotide enrichment', fontsize=14)
    ax.set_title('Unpredicted peaks', fontsize=14)
    ax.legend(fontsize=12)
    
    return fig

def add_intron_size(peaks_df, gff3, organism=None):
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    no_peaks = ss_dict
    intron_sizes = []
    for index, row in peaks_df.iterrows():
        if row['type'] != 'intronic':
            intron_sizes.append(np.NaN)
        else:
            sites = ss_dict[row['transcript']]
            assigned=False
            for pair in sites:
                if pair[0] > pair[1]:
                    if row['position'] >= pair[1] and row['position'] <= pair[0]:
                        intron_sizes.append(pair[0]-pair[1])
                        assigned=True
                        no_peaks[row['transcript']].remove(pair)
                        break
                else:
                    if row['position'] >= pair[0] and row['position'] <= pair[1]:
                        intron_sizes.append(pair[1]-pair[0])
                        assigned=True
                        no_peaks[row['transcript']].remove(pair)
                        break
            if assigned is False:
                intron_sizes.append(np.NaN)
    peaks_df['intron size'] = intron_sizes
    return peaks_df,  no_peaks


def count_reads_at_peaks(A_bam_list, totals_list, peak_df, B_bam_list=None, B_totals_list=None, gff3=None):
    '''Function to compare peak heights (or intermediate levels) between two genotypes give a set of peaks
    
    Parameters
    ----------
    A_bam_list : list
            bam files from cleavage profiling samples, [Mut, Mut, WT, WT]
    totals_list : list
            floats of total million aligned reads in each sample, [Mut, Mut, WT, WT]
    peak_df : pandas.DataFrame
            peak_df from SP_pipeline
    B_bam_list : list, default ``None``
            bam files from total spliceosome samples, [Mut, Mut, WT, WT]
            provide if calculating intermediate level rather than just peak heights.
    B_totals_list : list, default ``None``
            floats of total million aligned reads in each B sample, [Mut, Mut, WT, WT]
    gff3 : str, default ``None``
            provide if using B samples to calculate intermediate level
    
    Returns
    -------
    ratio_df : pandas.DataFrame
            dataframe with ratios and Z scores compared between samples'''
    
    A_bams = []
    for bam_file in A_bam_list:
        A_bams.append(pysam.Samfile(bam_file))

    #Count reads in transcript from B samples if provided
    if B_bam_list is not None:
        B_bams = []
        tx_dict = SPPeaks.build_transcript_dict(gff3)
        for bam_file in B_bam_list:
            B_bams.append(pysam.Samfile(bam_file))   
    
    #Count reads a peaks in dataframe
    all_counts = {}
    names = []
    n=0
    for n in range(len(A_bams)):
        name = A_bam_list[n].split('/')[-1].split('.')[0]
        names.append(name)
        print name
        all_counts[name] = []
        for ix, r in peak_df.iterrows():
            peak_count = 0
            pos = int(r['position'])
            peak_iter = A_bams[n].fetch(r['chromosome'],  pos-50,  pos+50)
            for read in peak_iter:
                if read.is_reverse and r['strand'] == '+':
                    if read.reference_end == pos:
                        peak_count += 1
                elif not read.is_reverse and r['strand'] == '-':
                    if read.reference_start == pos-1:
                        peak_count += 1
            
            if B_bam_list is None:
                all_counts[name].append(float(peak_count)/totals_list[n])
            else:
                tx_count = 0
                start = tx_dict[r['transcript']+'T0'][0]
                end = tx_dict[r['transcript']+'T0'][1]
                peak_iter = B_bams[n].fetch(r['chromosome'], start, end)
                for read in peak_iter:
                    if read.is_reverse and r['strand'] == '+':
                        tx_count += 1
                    elif not read.is_reverse and r['strand'] == '-':
                        tx_count += 1
                tx_rpkm = float(tx_count)/abs(start-end)/B_totals_list[n]
                all_counts[name].append((float(peak_count)/totals_list[n])/tx_rpkm)
    
    ratio_df = peak_df
    for k,v in all_counts.iteritems():
        ratio_df[k] = v
        ratio_df['log2_'+k] = [np.log2(x) for x in v]

    ratio_df['ratio1'] = ratio_df[names[0]]/ratio_df[names[2]]
    ratio_df['log_ratio1'] = ratio_df['ratio1'].apply(np.log)
    ratio_df['log_ratio1'] = ratio_df['log_ratio1'].replace([np.inf, np.inf*-1, np.NaN], 0)
    ratio_df['ratio2'] = ratio_df[names[1]]/ratio_df[names[3]]
    ratio_df['log_ratio2'] = ratio_df['ratio2'].apply(np.log)
    ratio_df['log_ratio2'] = ratio_df['log_ratio2'].replace([np.inf, np.inf*-1, np.NaN], 0)
    ratio_df['Z1'] = pd.Series(stats.mstats.zscore(ratio_df['log_ratio1']), index=ratio_df.index)
    ratio_df['Z2'] = pd.Series(stats.mstats.zscore(ratio_df['log_ratio2']), index=ratio_df.index)
    ratio_df['pvalue1'] = pd.Series(stats.norm.sf(abs(ratio_df['Z1']))*2)
    ratio_df['pvalue2'] = pd.Series(stats.norm.sf(abs(ratio_df['Z2']))*2)
    
    return ratio_df