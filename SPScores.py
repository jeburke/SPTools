import sys
import os
script_path = os.path.dirname(os.path.realpath(__file__)).split('SPTools')[0]
sys.path.append(script_path)
import SPTools as SP
import math
import numpy as np
import pandas as pd
import collections
from math import log
from scipy import stats
from matplotlib import pyplot as plt

def make_fasta_dict(fasta_file):
    fasta_dict = {}
    letters = ['A','C','G','T','N']
    n = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                n += 1
                line_list = line.split(' ')
                chr_num = line_list[0].strip()[1:]
                print chr_num
                fasta_dict[chr_num] = str()
            elif line[0] in letters and "chr"+str(n) in fasta_dict:
                fasta_dict["chr"+str(n)] = fasta_dict["chr"+str(n)]+line.strip()

    fasta_dict = collections.OrderedDict(sorted(fasta_dict.items()))
    return fasta_dict

def gc_content(fasta_dict):
    a = 0
    c = 0
    g = 0
    t = 0
    for chrom, seq in fasta_dict.iteritems():
        a += seq.count('A')
        c += seq.count('C')
        g += seq.count('G')
        t += seq.count('T')
    total = sum([a,c,g,t])
    nucleotide_prob = [float(a)/total, float(c)/total, float(t)/total, float(g)/total]
    return nucleotide_prob

def generate_consensus_matrix(gff3, fasta_dict, PSSM=False):
    #Populate gene dictionary and build genome
    if 'pombe' in gff3.lower():
        transcript_dict = SP.build_transcript_dict(gff3, organism='pombe')
        ss, flag = SP.list_splice_sites(gff3, organism='pombe')
        organism = 'pombe'
    else:
        transcript_dict = SP.build_transcript_dict(gff3)
        ss, flag = SP.list_splice_sites(gff3)
        organism = None
    ss_dict = SP.collapse_ss_dict(ss)
    genome = fasta_dict
    #print genome.keys()
    nuc_prob = gc_content(fasta_dict)
    #print nuc_prob

    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    #First generate a consensus matrix for 5' and 3' splice site, where 1st row is A counts, second row is C, third row is T, fourth row is G.
    pos_matrix_5prime = np.zeros([4,8])
    pos_matrix_3prime = np.zeros([4,8])

    counter1 = 0
    counter2 = 0

    for transcript, introns in ss_dict.iteritems():
        counter2 += 1
        if organism == 'pombe':
            isoform = transcript+'.1'
        else:
            isoform = transcript+'T0'
        strand = transcript_dict[isoform][2]
        chrom = transcript_dict[isoform][3]

        for intron in introns:
            counter1+=1
            if strand == '+':
                seq = fasta_dict[chrom][(intron[0]-1):(intron[0]+7)]
            elif strand == '-':
                seq = fasta_dict[chrom][(intron[0]-6):(intron[0]+2)]
                seq = SP.reverse_complement(seq)

            for a, base in enumerate(seq):
                pos_matrix_5prime[base_dict[base],a] += 1

            if strand == '+':
                seq = fasta_dict[chrom][(intron[1]-5):(intron[1]+3)]
            elif strand == '-':
                seq = fasta_dict[chrom][(intron[1]-2):(intron[1]+6)]
                seq = SP.reverse_complement(seq)
            
            for b, base in enumerate(seq):
                pos_matrix_3prime[base_dict[base],b] += 1
                
    #print counter1
    #print counter2

    float_formatter = lambda x: "%.1f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})
    
    a = 0
    while a < 4:
        b = 0
        while b < 8:
            if PSSM is False:
                pos_matrix_5prime[a,b] = (pos_matrix_5prime[a,b])/float(counter1)
                pos_matrix_3prime[a,b] = (pos_matrix_3prime[a,b])/float(counter1)
            if PSSM is True:
                if pos_matrix_5prime[a,b] == 0: pos_matrix_5prime[a,b] += 1
                if pos_matrix_3prime[a,b] == 0: pos_matrix_3prime[a,b] += 1
                pos_matrix_5prime[a,b] = np.log2((pos_matrix_5prime[a,b]/float(counter1))/nuc_prob[a])
                pos_matrix_3prime[a,b] = np.log2((pos_matrix_3prime[a,b]/float(counter1))/nuc_prob[a])
            b += 1
        a += 1
    
    return (pos_matrix_5prime, pos_matrix_3prime)
            
def score_peaks(df, gff3, fasta_dict):
    pos_matrix_5prime, pos_matrix_3prime = generate_consensus_matrix(gff3, fasta_dict, PSSM=True)
    score = []
    ann_score5 = []
    ann_score3 = []
    
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    for ix, r in df.iterrows():
        strand = r['strand']
        seq = r['sequence']
        ann_seq3 = r['annotated sequence2']
        ann_seq5 = r['annotated sequence1']

        r_score = 0
        r_ann_score5 = 0
        r_ann_score3 = 0
        if r['looks like'] != 'AG' and r['looks like'] != '3prime':
            for a, base in enumerate(seq[4:]):
                r_score += pos_matrix_5prime[base_dict[base],a]
        else:
            for b, base in enumerate(seq[:8]):
                r_score += pos_matrix_3prime[base_dict[base],b]
        score.append(r_score)
        
        if ann_seq5 is not None:
            for a, base in enumerate(ann_seq5):
                r_ann_score5 += pos_matrix_5prime[base_dict[base],a]
            ann_score5.append(r_ann_score5)
        else:
            ann_score5.append(np.NaN)
                
        if ann_seq3 is not None:
            for a, base in enumerate(ann_seq3):
                r_ann_score3 += pos_matrix_3prime[base_dict[base],a]
            ann_score3.append(r_ann_score3)
        else:
            ann_score3.append(np.NaN)
    
    df['score'] = score
    df['annotated 5p score'] = ann_score5
    df['annotated 3p score'] = ann_score3
        
    return df
    
def simple_score_junction(seq5, seq3, PSSM):
    '''Function to score just one intron. 
    Need PSSM from generate_consensus_matrix.
    Sequence lengths should match PSSM (8 each).'''
    
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    score5 = 0
    score3 = 0
    for a, base in enumerate(seq5):
        score5 += PSSM[0][base_dict[base],a]
    for b, base in enumerate(seq3):
        score3 += PSSM[1][base_dict[base],b] 
    
    return score5, score3

## Find and score potential branches
def generate_all_branches():
    branch_options = [('C','T'),
                      ('T','A'),
                      ('G','A','C'),
                      ('A'),
                      ('C','T')]
    branches = {}
    for n, position in enumerate(branch_options):
        branches[n] = set()
        for base in position:
            if n == 0:
                branches[n].add(base)
            else:
                for branch in branches[n-1]:
                    branches[n].add(branch+base)
    branches = branches[n]
    
    return branches

def percent_py(seq):
    py_dict = {'A':0,'G':0,'T':1,'C':1,'N':0}
    score = 0
    for base in seq:
        score += py_dict[base]
    score = float(score)/len(seq)
    return score

def branch_PSSM(peak_branch_df, fa_dict):
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    nuc_prob = gc_content(fa_dict)
    
    pos_matrix_branch = np.zeros([4,5])
    counter = 0
    if type(peak_branch_df) == str:
        with open(peak_branch_df) as f:
            for line in f:
                counter += 1
                seq = line.strip()
                for a, base in enumerate(seq):
                    pos_matrix_branch[base_dict[base],a] += 1
    else:
        for seq in peak_branch_df['Branch seq']:
            counter += 1
            seq = seq[2:7]
            for a, base in enumerate(seq):
                pos_matrix_branch[base_dict[base],a] += 1

    float_formatter = lambda x: "%.1f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})
    
    a = 0
    while a < 4:
        b = 0
        while b < 5:
            if pos_matrix_branch[a,b] == 0: pos_matrix_branch[a,b] += 1
            pos_matrix_branch[a,b] = np.log2((pos_matrix_branch[a,b]/float(counter))/nuc_prob[a])
            b += 1
        a += 1
    
    return pos_matrix_branch

def score_branch(seq, PSSM):
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    branch_score = 0
    for a, base in enumerate(seq):
        branch_score += PSSM[base_dict[base],a]
    return branch_score

def find_score_branches_ppy(quant_df, peak_branch_df, fa_dict):
    #branches = generate_all_branches()
    PSSM = branch_PSSM(peak_branch_df, fa_dict)
    
    branches = []
    if type(peak_branch_df) is not str:
        for ix, branch in peak_branch_df['Branch seq'].iteritems():
            seq = branch[2:7]
            if seq[:-2] != 'A' and 'A' in seq:
                A_ix = branch.rfind('A')
                new_seq = branch[A_ix-3:A_ix+2]
                if len(new_seq) == 5:
                    seq = new_seq
        branches = peak_branch_df['Branch seq'].str[2:7]
    else:
        with open(peak_branch_df) as f:
            for line in f:
                branches.append(line.strip())
    
    # Sort branches by abundance so that the most common ones are first in the search
    br_abund = []
    for branch in set(branches):
        count = len([x for x in branches if x == branch])
        if count > 1:
            br_abund.append((branch, count))
    br_abund = sorted(br_abund, key=lambda x: x[1], reverse=True)
    branches = zip(*br_abund)[0]
    
    branch_dict = collections.OrderedDict()
    for branch in branches:
        branch_dict[branch] = score_branch(branch, PSSM)
    
    branch_3_dist = []
    branch_score = []
    branch_seqs = []
    perc_py = []
    for ix, r in quant_df.iterrows():
        if r['strand'] == '+':
            intron_seq = fa_dict[r['chromosome']][int(r['position']):int(r['position']+r['intron size'])]
            three_site = r['position']+r['intron size']
        elif r['strand'] == '-':
            intron_seq = fa_dict[r['chromosome']][int(r['position']-r['intron size']-1):int(r['position']-1)]
            intron_seq = SP.reverse_complement(intron_seq)
            three_site = r['position']-r['intron size']
        
        if type(peak_branch_df) is not str:
            if ix in peak_branch_df['genome coord']:
                ix_df = peak_branch_df[peak_branch_df['genome coord'] == ix]
                ix_df = ix_df.sort_values('depth', ascending=False)
                best_branch = ix_df.iloc[0,'branch site']
                best_branch = abs(ix_df.iloc[0,'5p splice site']-best_branch)

                seq = ix_df.iloc[0,'Branch seq'][2:7]
                branch_seqs.append(seq)
                branch_score.append(score_branch(seq, PSSM))
                branch_3_dist.append(ix_df.iloc[0,'Branch to 3p distance'])

                if 'N' in intron_seq[best_branch[0]+5:]:
                        print ix
                        print intron_seq
                perc_py.append(percent_py(intron_seq[best_branch[0]+5:]))
                               
            else:
                matches = []
                for branch in branch_dict:
                    if branch in intron_seq:
                        matches.append((intron_seq.index(branch), branch, branch_dict[branch]))

                if len(matches) == 0:
                    # Find the closest A
                    best_ix = intron_seq[:-3].rfind('A')
                    seq = intron_seq[best_ix-3:best_ix+2]
                    score = score_branch(seq, PSSM)
                    best_branch = (best_ix, seq, score)

                elif len(matches) > 1:
                    matches = sorted(matches, key=lambda x: x[2], reverse=True)
                    best_branch = matches[0]
                else:
                    best_branch = matches[0]

                branch_3_dist.append((len(intron_seq)-best_branch[0]-4)/1000.)
                branch_score.append(best_branch[2])
                branch_seqs.append(best_branch[1])

                if len(intron_seq)-best_branch[0]-5 > 1:
                    if 'N' in intron_seq[best_branch[0]+5:]:
                        print ix
                        print intron_seq
                    perc_py.append(percent_py(intron_seq[best_branch[0]+5:]))
                else:
                    perc_py.append(np.NaN)
        else:
            matches = []
            for branch in branch_dict:
                if branch in intron_seq:
                    matches.append((intron_seq.index(branch), branch, branch_dict[branch]))

            if len(matches) == 0:
                # Find the closest A
                best_ix = intron_seq[:-3].rfind('A')
                seq = intron_seq[best_ix-3:best_ix+2]
                score = score_branch(seq, PSSM)
                best_branch = (best_ix, seq, score)

            elif len(matches) > 1:
                matches = sorted(matches, key=lambda x: x[2], reverse=True)
                best_branch = matches[0]
            else:
                best_branch = matches[0]

            branch_3_dist.append((len(intron_seq)-best_branch[0]-4)/1000.)
            branch_score.append(best_branch[2])
            branch_seqs.append(best_branch[1])

            if len(intron_seq)-best_branch[0]-5 > 1:
                if 'N' in intron_seq[best_branch[0]+5:]:
                    print ix
                    print intron_seq
                perc_py.append(percent_py(intron_seq[best_branch[0]+5:]))
            else:
                perc_py.append(np.NaN)
    
    quant_df['branch score'] = branch_score
    quant_df['branch to 3p distance'] = branch_3_dist
    quant_df['percent pPy'] = perc_py
    
    branch_seqs = ['NNNNN' if len(x) < 5 else x for x in branch_seqs ]
    
    for n in range(len(branch_seqs[0])):
        pos = [x[n] for x in branch_seqs]
        quant_df['branch-'+str(n)] = pos
    
    #print str(len(quant_df)-len(quant_df['branch score'].dropna()))+' introns without identifiable branches'
    return quant_df


##########################################################################
## Functions for analyzing splice sites compared to all annotated sites ##
##########################################################################

def generate_all_ss_seqs(gff3, fasta_dict, organism):
    transcript_dict = SP.build_transcript_dict(gff3, organism=organism)
    ss, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss)
    
    all_seq5 = []
    all_seq3 = []
    for transcript, introns in ss_dict.iteritems():
        if organism == 'pombe':
            isoform = transcript+'.1'
        else:
            isoform = transcript+'T0'
        strand = transcript_dict[isoform][2]
        chrom = transcript_dict[isoform][3]

        for intron in introns:
            if strand == '+':
                seq5 = fasta_dict[chrom][(intron[0]-1):(intron[0]+7)]
            elif strand == '-':
                seq5 = fasta_dict[chrom][(intron[0]-6):(intron[0]+2)]
                seq5 = SP.reverse_complement(seq5)

            all_seq5.append(seq5)

            if strand == '+':
                seq3 = fasta_dict[chrom][(intron[1]-5):(intron[1]+3)]
            elif strand == '-':
                seq3 = fasta_dict[chrom][(intron[1]-2):(intron[1]+6)]
                seq3 = SP.reverse_complement(seq3)
            
            all_seq3.append(seq3)
    return all_seq5, all_seq3

def seq_list_to_totals(seq_list):
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    seq_list = [x for x in seq_list if x is not None]
    
    total_array = np.ones([4, len(seq_list[0])])
    
    n=0
    for n in range(len(seq_list)):
        for a, base in enumerate(seq_list[n]):
            total_array[base_dict[base],a] += 1
    return total_array

def position_wise_scores2(seq5_list, seq3_list, organism, title='Intron position strength'):
    '''Uses chi-contingency test to score base proportions at each position in sample against population'''
    
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)

    all_5p, all_3p = generate_all_ss_seqs(gff3, fa_dict, organism)
    
    pop_5p = seq_list_to_totals(all_5p)
    pop_3p = seq_list_to_totals(all_3p)
    samp_5p = seq_list_to_totals(seq5_list)
    samp_3p = seq_list_to_totals(seq3_list)
    print samp_5p.shape

    p5 = []
    for n in range(samp_5p.shape[1]):
        if n == 2 or n == 3:
            p5.append(1)
        else:
            conting = np.array([samp_5p[:,n],pop_5p[:,n]])
            chi2, p, dof, expected = stats.chi2_contingency(conting)
            p5.append(np.log10(p)*-1)
        
    p3 = []
    for n in range(samp_3p.shape[1]):
        if n == 4 or n == 5:
            p3.append(1)
        else:
            conting = np.array([samp_3p[:,n],pop_3p[:,n]])
            chi2, p, dof, expected = stats.chi2_contingency(conting)
            p3.append(np.log10(p)*-1)
    
    fig, ax = plt.subplots(2, 1, figsize=(4,4))
    width = 0.7
    
    max_y = max(p5+p3) + 0.1*max(p5+p3)
    
    ind5 = np.arange(len(p5))
    ax[0].bar(ind5, p5, color='k')
    ax[0].plot([0,8], [2,2], '--', color='0.7')
    ax[0].set_xlim([0,len(p5)])
    ax[0].set_ylabel("5' splice site\n-log10(p-value)")
    ax[0].set_title(title)
    ax[0].set_ylim([0,max_y])

    ind3 = np.arange(len(p3))
    ax[1].bar(ind3, p3, color='k')
    ax[1].plot([0,8], [2,2], '--', color='0.7')
    ax[1].set_xlim([0,len(p3)])
    ax[1].set_ylabel("3' splice site\n-log10(p-value)")
    ax[1].set_ylim([0,max_y])

    ax[0].set_xticks(ind3 + width / 2)
    ax[1].set_xticks(ind3 + width / 2)
    ax[0].set_xticklabels(np.arange(-2,6))
    ax[1].set_xticklabels(np.arange(-5,3))

    fig.tight_layout()
    plt.show()
    return fig