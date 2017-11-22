'''Script for calling 5' and 3' splice sites and relative strengths of each sites 

Usage: python Call_splice_site_motifs.py

This will generate an output spreadsheet with relative splice-site strength for each transcript in the genome'''

import sys
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
import GeneUtility
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
import SPTools as SP
import collections
from math import log
from scipy import stats
from matplotlib import pyplot as plt

def read_juncbase_output(juncbase_output):
    if type(juncbase_output) == list:
        df_list = []
        for output in juncbase_output:
            with open(juncbase_output, 'r') as fin:
                df = pd.read_csv(fin, sep='\t')
                df_list.append(df)
        junc_df = pd.concat(df_list)
        junc_df.replace(to_replace='NA', value=np.NaN, inplace=True)
    elif type(juncbase_output) == str:
        junc_df = pd.read_csv(juncbase_output, sep='\t')
        if not junc_df.columns[0].startswith('#'):
            junc_df = junc_df.drop(junc_df.columns[0], axis=1)
    
    sample_list = []
    junc_df['contained in'] = None
    n_max = len(junc_df.columns)
    n = 11
    while n < n_max:
        if junc_df.columns[n].startswith('avg'): continue
        elif junc_df.columns[n] == 'contained in': n += 1
        elif n == 11:
            sample_list.append(junc_df.columns[n])
            sample_list.append(junc_df.columns[n+1])
            junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]] = junc_df[junc_df.columns[n]]
            wt = junc_df.columns[n]
            n += 2
        else:
            if junc_df.columns[n+1] != 'contained in':
                print 'Averaging '+junc_df.columns[n]+' and '+junc_df.columns[n+1]
                junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]] = junc_df[[(junc_df.columns[n]),(junc_df.columns[n+1])]].mean(axis=1, skipna=False)
                
                sample_list.append(junc_df.columns[n])
                sample_list.append(junc_df.columns[n+1])
                if str(junc_df.columns[n]) == 'nan' or str(junc_df.columns[n+1]) == 'nan':
                    print junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]]
                junc_df[junc_df.columns[n]+'+'+junc_df.columns[n+1]] = ''
            n += 2
    
    #print junc_df.columns
    filt_df = pd.DataFrame(columns=junc_df.columns, index = junc_df.index)
    n = 0
    for n in range(len(junc_df)):
        sample_str = ''
        flag = False
        row = junc_df.iloc[[n]]
        a = 0
        for column in row.columns:
            if column.startswith('avg'):
                if a == 0:
                    wt_name = column
                elif a > 0:
                    if str(junc_df[column][n]) == 'nan':
                        junc_df[column.split('avg_')[1]] = np.NaN
                    else:
                        flag = True
                        sample_str = sample_str+column.split('avg_')[1]+', '
                        if junc_df[column][n] > junc_df[wt_name][n]:
                            junc_df[column.split('avg_')[1]] = 'Up'
                        elif junc_df[column][n] < junc_df[wt_name][n]:
                            junc_df[column.split('avg_')[1]] = 'Down'
                a += 1                   

        filt_df.iloc[[n]] = junc_df.iloc[[n]]
        filt_df.loc[n,'contained in'] = sample_str

    filt_df.replace(to_replace=np.NaN, value='NA', inplace=True)
    filt_df = filt_df[filt_df['strand'] != 'NA']
    filt_df = filt_df.reset_index()
    
    n = 0
    for n in range(len(filt_df)):
        for column in filt_df.columns:
            if column in sample_list:
                if filt_df[column][n] == 'NA':
                    filt_df.loc[n, column] = filt_df[wt][n]
                    
    
    filt_df['coord_1'] = filt_df['exclusion_junctions'].str.split(',').str.get(0).str.split(':').str.get(1).str.split('-').str.get(0)
    filt_df['coord_2'] = filt_df['exclusion_junctions'].str.split(',').str.get(0).str.split(':').str.get(1).str.split('-').str.get(1).str.split(';').str.get(0)
    
    print sample_list
    return (filt_df, sample_list)

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


##Input: dictionary of coordinates by transcript (transcript:[list of coordinates][list of coordinates])
def get_sequence(coord_dict, gff3_file, fasta_file):
    if 'pombe' in gff3_file:
        organism = 'pombe'

    else: organism = None
    
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta_file) is str:
        fasta_dict = make_fasta_dict(fasta_file)
    else:
        fasta_dict = fasta_file
    
    seq_dict = {}
    counter5 = 0
    counter3 = 0
    other = 0
    for transcript, coord_sets in coord_dict.iteritems():
        seq_dict[transcript] = []
        chrom = transcript_dict[transcript][3]
        #if chrom in rom_lat: chrom = rom_lat[chrom]
        strand = transcript_dict[transcript][2]
        for coord in coord_sets[0]:
            seq_type = 'other'
            if strand == "+":
                sequence = fasta_dict[chrom][(coord-9):(coord+11)]
            elif strand == "-":
                sequence = fasta_dict[chrom][(coord-10):(coord+10)]
                sequence = SP.reverse_complement(sequence)

            if sequence[10:12] == 'GT' or sequence[10:12] == 'GC': 
                seq_type = "5'"
                counter5 += 1
            seq_dict[transcript].append((sequence, seq_type))
     
        for coord in coord_sets[1]:
            seq_type = 'other'
            if strand == "+":
                sequence = fasta_dict[chrom][(coord-9):(coord+11)]
            elif strand == "-":
                sequence = fasta_dict[chrom][(coord-10):(coord+10)]
                sequence = SP.reverse_complement(sequence)
                
            if sequence[8:10] == 'AG': 
                seq_type = "3'"
                counter3 += 1
            seq_dict[transcript].append((sequence, seq_type))
    
    #print str(counter5)+" 5' splice sites"
    #print str(counter3)+" 3' splice sites"
    
    return seq_dict
        
def get_junction_sequence(df, gff3_file, fasta_file):
    df = df.sort_values('chr', axis=0)
    
    #transcript_dict[transcript] = [start, end, strand, chromosome, CDS start, CDS end]
    transcript_dict = SP.build_transcript_dict(gff3_file)

    #splice_dict[transcipt] = [[5'sites][3'sites]]
    splice_dict, flag = SP.list_splice_sites(gff3_file)
    
    #fasta_dict[chr] = sequence
    if type(fasta_file) is str:
        fasta_dict = make_fasta_dict(fasta_file)
    else:
        fasta_dict = fasta_file

    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)

    df['Gene'] = "Unknown"
    df['intron'] = "Middle"
    df['sequence1'] = ''
    df['sequence2'] = ''
    df['intron sequence'] = 'No sequence here'

    n = 0
    for n in range(len(df)):
        coord1 = int(df['coord_1'][n].strip())
        coord2 = int(df['coord_2'][n].strip())
        chrom = df['chr'][n].strip()
        strand = df['strand'][n].strip()
        transcripts = transcript_by_chr[chrom]

        for transcript in transcripts:
            tx_strand = transcript_dict[transcript][2]
            start = transcript_dict[transcript][0]
            stop = transcript_dict[transcript][1]
            
            if strand == tx_strand and coord1 >= start and coord2 <= stop:
                df.loc[n,'Gene'] = transcript
               
        if strand == '+':
            sequence1 = fasta_dict[chrom][(coord1-3):(coord1+5)]
            sequence2 = fasta_dict[chrom][(coord2-6):(coord2+2)]
            all_seq = fasta_dict[chrom][(coord1-1):coord2]
        elif strand == '-':
            sequence1 = fasta_dict[chrom][(coord2-6):(coord2+2)]
            sequence1 = SP.reverse_complement(sequence1)
            sequence2 = fasta_dict[chrom][(coord1-3):(coord1+5)]
            sequence2 = SP.reverse_complement(sequence2)
            all_seq = fasta_dict[chrom][(coord1-1):coord2]
            all_seq = SP.reverse_complement(all_seq)
        
        df.loc[n,'sequence1'] = sequence1
        df.loc[n,'sequence2'] = sequence2
        df.loc[n,'intron sequence'] = all_seq

    for transcript in transcripts:
        if transcript in df['Gene'].tolist():
            tx_df = df[df['Gene'] == transcript]
            s = tx_df['coord_1']
            min_idx = s.idxmin()
            first = int(s.min())
            #print transcript_dict[transcript][2]
            #print first
            max_idx = s.idxmax()
            last = int(s.max())
            #print last
        
            if first == last:
                df.loc[min_idx,'intron'] = 'Only'
            else:
                if transcript_dict[transcript][2] == '+':
                    df.loc[min_idx,'intron'] = 'First'
                    df.loc[max_idx,'intron'] = 'Last'
                elif transcript_dict[transcript][2] == '-':
                    df.loc[min_idx,'intron'] = 'Last'
                    df.loc[max_idx,'intron'] = 'First'
            
            for index, coord_1 in s.iteritems():
                if df['intron'][index] == 'Middle':
                    if coord_1 in range(first-10, first+10):
                        df_idx = s[s == coord_1].index[0]
                        if transcript_dict[transcript][2] == '+':
                            df.loc[df_idx, 'intron'] = 'First'
                        elif transcript_dict[transcript][2] == '-':
                            df.loc[df_idx, 'intron'] = 'Last'
                    elif coord_1 in range(last-10, last+10):
                        df_idx = s[s == coord_1].index[0]
                        if transcript_dict[transcript][2] == '+':
                            df.loc[df_idx, 'intron'] = 'Last'
                        elif transcript_dict[transcript][2] == '-':
                            df.loc[df_idx, 'intron'] = 'First'
                
    df = df[df['contained in'] != '']
    df = df.reset_index()
    return df

def mirrored(listtree):
    return [(mirrored(x) if isinstance(x, list) else x) for x in reversed(listtree)]

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


def generate_PSSM(seq_list, fasta_dict):
    #Populate gene dictionary and build genome
    genome = fasta_dict
    nuc_prob = gc_content(fasta_dict)

    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    #First generate a consensus matrix for the sequence, where 1st row is A counts, second row is C, third row is T, fourth row is G.
    PSSM = np.zeros([4,len(seq_list[0])])

    counter = 0
    for seq in seq_list:
        counter += 1
        for a, base in enumerate(seq):
            PSSM[base_dict[base],a] += 1

    float_formatter = lambda x: "%.1f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})
    
    a = 0
    while a < 4:
        b = 0
        while b < len(seq_list[0]):
            if PSSM[a,b] == 0: PSSM[a,b] += 1
            PSSM[a,b] = np.log2((PSSM[a,b]/float(counter))/nuc_prob[a])
            b += 1
        a += 1
    
    return PSSM

def score_new_sites(df, pos_matrix_5prime, pos_matrix_3prime, PSSM=False):
    #Compare each splice site and add scores to dataframe - needs columns named 'coord_1', 'coord_2', 'sequence1', 'sequence2'
    df['5p score'] = ''
    df['3p score'] = ''
    df['intron length'] = 'NA'
    
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    n = 0
    for n in range(len(df)):
        strand = df['strand'][n]
        if strand == '+':
            coord1 = int(df['coord_1'][n])
            coord2 = int(df['coord_2'][n])
        elif strand == '-':
            coord1 = int(df['coord_2'][n])
            coord2 = int(df['coord_1'][n])
        seq5 = df['sequence1'][n]
        seq3 = df['sequence2'][n]
        
        if PSSM is False:
            gene_matrix_5prime = np.zeros([4,8])
            for a, base in enumerate(seq5):
                gene_matrix_5prime[base_dict[base],a] += 1

            gene_matrix_3prime = np.zeros([4,20])
            for b, base in enumerate(seq3):
                gene_matrix_3prime[base_dict[base],b] += 1

            #Calculate Scores (score of 0 is perfect, higher score is worse, 40 is highest possible)
            score_5prime = 0
            score_3prime = 0
            a = 0
            while a < 4:
                b = 0
                while b < 8:
                    score_5prime += abs(pos_matrix_5prime[a,b] - (gene_matrix_5prime[a,b]))
                    score_3prime += abs(pos_matrix_3prime[a,b] - (gene_matrix_3prime[a,b]))
                    b += 1
                a += 1

        elif PSSM is True:
            score_5prime = 0
            score_3prime = 0
            for a, base in enumerate(seq5):
                score_5prime += pos_matrix_5prime[base_dict[base],a]
            for b, base in enumerate(seq3):
                score_3prime += pos_matrix_3prime[base_dict[base],b]                  
                
        intron_length = abs(coord2-coord1)
    
        df.loc[n,'5p score'] = score_5prime
        df.loc[n,'3p score'] = score_3prime
        df.loc[n,'intron length'] = intron_length
        
    return df

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

def reformat_df(df, sample_list):
    new_columns1=['Gene', 'as_event_type', 'chr', 'coord_1', 'coord_2', 'sequence1', 'sequence2', 'intron length', 'intron sequence', '5p score', '3p score', 'strand', 'intron', '#Contains_Novel_or_Only_Known(Annotated)_Junctions', 'contained in']
    new_df1 = df[new_columns1]
    new_df1 = new_df1.sort_values('as_event_type', axis=0)
    new_df1 = new_df1.rename(columns = {'as_event_type':'event','sequence1':'5p sequence','sequence2':'3p sequence','#Contains_Novel_or_Only_Known(Annotated)_Junctions':'known or novel'})
    new_df1 = new_df1.reset_index()
    new_df1 = new_df1.drop('index', axis=1)
    new_df1[['5p score','3p score']] = new_df1[['5p score','3p score']].astype(float)
    
    new_columns2 = ['Gene', 'as_event_type', 'chr', 'coord_1', 'coord_2', 'strand']+sample_list
    new_columns2 = new_columns2+[x for x in df.columns if '+' in x]
    new_df2 = df[new_columns2]
    new_df2 = new_df2.sort_values('as_event_type', axis=0)
    new_df2 = new_df2.rename(columns = {'as_event_type':'event'})
    new_df2 = new_df2.reset_index()
    new_df2 = new_df2.drop('index', axis=1)
    new_df2.replace(to_replace='NA', value=np.NaN, inplace=True)
    new_df2[sample_list] = new_df2[sample_list].astype(float)
    
    return new_df1, new_df2

def write_seq_list_to_file(df, prefix):
    with open('{0}_5pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['sequence1'].tolist():
            f.write(seq+'\n')
    with open('{0}_3pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['sequence2'].tolist():
            f.write(seq+'\n')
            
def percent_py(seq):
    py_dict = {'A':0,'G':0,'T':1,'C':1,'N':0}
    score = 0
    for base in seq:
        score += py_dict[base]
    score = float(score)/len(seq)
    return score
        
def score_PyTract(df, fa_dict, alt_column_name=None, from_branches=False):
    py_score1 = []
    py_score2 = []
    alt_py1 = []
    alt_py2 = []
    
    for ix, r in df.iterrows():
        strand = r['strand']
        chrom = r['chromosome']
        coord = r['annotated intron coords'][1]
        alt_coord = r['junction coords'][1]
        if strand == '+':
            if coord is not None:
                seq1 = fa_dict[chrom][coord-15:coord]
                seq2 = fa_dict[chrom][coord-30:coord-15]
            alt1 = fa_dict[chrom][alt_coord-15:alt_coord]
            alt2 = fa_dict[chrom][alt_coord-30:alt_coord-15]
        if strand == '-':
            if coord is not None:
                seq1 = fa_dict[chrom][coord:coord+15]
                seq2 = fa_dict[chrom][coord+15:coord+30]
                seq1 = SP.reverse_complement(seq1)
                seq2 = SP.reverse_complement(seq2)
            alt1 = fa_dict[chrom][alt_coord:alt_coord+15]
            alt2 = fa_dict[chrom][alt_coord+15:alt_coord+30]
            alt1 = SP.reverse_complement(alt1)
            alt2 = SP.reverse_complement(alt2)

        alt_py1.append(percent_py(alt1))
        alt_py2.append(percent_py(alt2))
        
        if coord is not None:
            py_score1.append(percent_py(seq1))
            py_score2.append(percent_py(seq2))
        else:
            py_score1.append(np.NaN)
            py_score2.append(np.NaN)
    
    df['Py score annotated -15:0'] = py_score1
    df['Py score annotated -30:-15'] = py_score2
    df['Py score alternative -15:0'] = alt_py1
    df['Py score alternative -30:-15'] = alt_py2
    return df
       
    
def position_wise_scores(seq_5p, seq_3p, organism):
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    PSSM_5p, PSSM_3p = generate_consensus_matrix(gff3, fa_dict, PSSM=True)
    base_dict = {"A":0, "C":1, "T":2, "G":3}
    
    seq_5p = [x for x in seq_5p if x is not None]
    seq_3p = [x for x in seq_3p if x is not None]
    
    score_5prime = np.empty([2,len(seq_5p[0])])
    score_3prime = np.empty([2,len(seq_3p[0])])
    all_5p = np.empty([len(seq_5p), len(seq_5p[0])])
    all_3p = np.empty([len(seq_3p), len(seq_3p[0])])
                       
    n=0
    for n in range(len(seq_5p)):
        for a, base in enumerate(seq_5p[n]):
            all_5p[n,a] = PSSM_5p[base_dict[base], a]
    
    a=0
    for a in range(len(score_5prime[0])):
        score_5prime[0,a] = np.median(all_5p[0:,a])
        score_5prime[1,a] = (max(all_5p[0:,a])-min(all_5p[0:,a]))/2.
    print score_5prime
            
    m=0
    for m in range(len(seq_3p)):
        for b, base in enumerate(seq_3p[m]):
            all_3p[m,b] = PSSM_3p[base_dict[base], b]
        
    b=0
    for b in range(len(score_3prime[0])):
        score_3prime[0,b] = np.median(all_3p[0:,b])
        score_3prime[1,b] = (max(all_3p[0:,b])-min(all_3p[0:,b]))/2.
    print score_3prime
    
    return all_5p, all_3p

from scipy import stats

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

def by_pos_plots(df, metrics=['Intermediate Level', 'Precursor']):
    col5 = [x for x in df.columns if 'Base 5' in x[1]]
    col3 = [x for x in df.columns if 'Base 3' in x[1]]
    
    for direction in ['Up','Down']:
        for metric in metrics:
            if len(df[df[('All',metric+' change')] == direction]) > 5:
                for n in range(len(col5)):
                    if n == 0:
                        s5 = df[df[('All',metric+' change')] == direction][col5[n]]
                    else:
                        s5 = s5.str.cat(df[df[('All',metric+' change')] == direction][col5[n]])
                print len(s5)

                for n in range(len(col3)):
                    if n == 0:
                        s3 = df[df[('All',metric+' change')] == direction][col3[n]]
                    else:
                        s3 = s3.str.cat(df[df[('All',metric+' change')] == direction][col3[n]])

                print metric+' '+direction
                fig = SP.position_wise_scores2(s5, s3, 'crypto')

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

                    #branch_3_dist.append(np.NaN)
                    #branch_score.append(np.NaN)
                    #branch_seqs.append('NNNNN')
                    #perc_py.append(percent_py(intron_seq[-30:]))
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

                #branch_3_dist.append(np.NaN)
                #branch_score.append(np.NaN)
                #branch_seqs.append('NNNNN')
                #perc_py.append(percent_py(intron_seq[-30:]))
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
    
    print len(quant_df)
    print len(branch_score)
    
    quant_df['branch score'] = branch_score
    quant_df['branch to 3p distance'] = branch_3_dist
    quant_df['percent pPy'] = perc_py
    
    branch_seqs = ['NNNNN' if len(x) < 5 else x for x in branch_seqs ]
    print len(branch_seqs)
    print len(quant_df)
    
    for n in range(len(branch_seqs[0])):
        pos = [x[n] for x in branch_seqs]
        quant_df['branch-'+str(n)] = pos
    
    print str(len(quant_df)-len(quant_df['branch score'].dropna()))+' introns without identifiable branches'
    return quant_df              