import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome')
sys.path.append('/home/jordan/RNA-is-awesome')
import SPTools as SP
import pandas as pd
import numpy as np
import pysam
import copy
from scipy import stats
import seaborn as sns
from matplotlib import pyplot as plt

def backfill_splice_sites(df, gff3, fa_dict, PSSM, organism=None):
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)

    column_dict = {'position':[],'transcript':[],'alt splicing':[],'type':[],'strand':[], 'introns in transcript':[],
                   'intron size':[],'chromosome':[], '5p score':[], '3p score':[], 'intron position':[], 'exon size (us)':[],
                   'exon size (ds)':[],'transcript size':[], 'peak':[], 'seq5':[],'seq3':[]} 
    new_index = []
    
    for tx in set(df['transcript']):
        strand = df[df['transcript'] == tx].iloc[0]['strand']
        splice_sites = ss_dict[tx]
        if strand == '+':
            splice_sites = sorted(list(splice_sites), key=lambda x:x[0])
        elif strand == '-':
            splice_sites = sorted(list(splice_sites), key=lambda x:x[0], reverse=True)
        
        df_pos = None
        for n, (five_site, three_site) in enumerate(splice_sites):
            # Check if already in dataframe
            in_df = False
            for peak in df[df['transcript'] == tx]['position']:
                if five_site in range(int(peak)-5,int(peak)+5):
                    in_df = True
                    df_pos = peak
                    break
            
            column_dict['transcript'].append(tx)
            if organism == 'pombe':
                iso = tx+'.1'
            else: iso = tx+'T0'
            
            column_dict['intron size'].append(abs(three_site-five_site))
            column_dict['introns in transcript'].append(len(splice_sites))
            column_dict['strand'].append(strand)   
            chrom = df[df['transcript'] == tx].iloc[0]['chromosome']
            column_dict['chromosome'].append(chrom)
            column_dict['transcript size'].append((tx_dict[iso][1]-tx_dict[iso][0])/1000.)

            # Check if first or last intron and add exon size
            if n == 0:
                column_dict['intron position'].append('First')
                if strand == '+':
                    column_dict['exon size (us)'].append((five_site-tx_dict[iso][0])/1000.)
                    if len(splice_sites) > 1:
                        ds_length = (splice_sites[n+1][0] - three_site)/1000.
                        try:
                            if ds_length < 0:
                                ds_length = (splice_sites[n+2][0] - three_site)/1000.
                        except IndexError:
                            ds_length = np.NaN
                    else:
                        ds_length = (tx_dict[iso][1] - three_site)/1000.
                    
                elif strand == '-':
                    column_dict['exon size (us)'].append((tx_dict[iso][1]-five_site)/1000.)
                    if len(splice_sites) > 1:
                        ds_length = (three_site - splice_sites[n+1][0])/1000.
                        try:
                            if ds_length < 0:
                                ds_length = (three_site - splice_sites[n+2][0])/1000.
                        except IndexError:
                            ds_length = np.NaN
                    else:
                        ds_length = (three_site - tx_dict[iso][0])/1000.
                column_dict['exon size (ds)'].append(ds_length)
            
            elif n == len(splice_sites)-1:
                column_dict['intron position'].append('Last')
                column_dict['exon size (us)'].append((abs(five_site-splice_sites[n-1][1])-1)/1000.)
                
                if strand == '+':
                    column_dict['exon size (ds)'].append((tx_dict[iso][1]-three_site)/1000.)
                elif strand == '-':
                    column_dict['exon size (ds)'].append((three_site - tx_dict[iso][0])/1000.)
            else:
                column_dict['intron position'].append('Middle')
                column_dict['exon size (us)'].append((abs(five_site-splice_sites[n-1][1])-1)/1000.)
                column_dict['exon size (ds)'].append(abs(three_site - splice_sites[n+1][0])/1000.)

            if in_df is True:
                peak_index = chrom+':'+str(int(df_pos))
                new_index.append(peak_index)
                column_dict['position'].append(df_pos)
                column_dict['3p score'].append(df.loc[peak_index,'3p score'])
                column_dict['5p score'].append(df.loc[peak_index,'5p score'])
                column_dict['alt splicing'].append(df.loc[peak_index,'alt splicing'])
                column_dict['type'].append(df.loc[peak_index,'type'])
                column_dict['peak'].append(True)
                column_dict['seq5'].append(df.loc[peak_index,'seq5'])
                column_dict['seq3'].append(df.loc[peak_index,'seq3'])

            if in_df is False:
                column_dict['alt splicing'].append(False)
                column_dict['type'].append('5prime')
                column_dict['peak'].append(False)
                
                # Get position, index and sequence for scoring and position code
                if strand == '+':
                    column_dict['position'].append(five_site+1)
                    new_index.append(chrom+':'+str(five_site+1))
                    sequence1 = fa_dict[chrom][(five_site-1):(five_site+7)]
                    sequence2 = fa_dict[chrom][(three_site-5):(three_site+3)]
                
                elif strand == '-':
                    column_dict['position'].append(five_site-1)
                    new_index.append(chrom+':'+str(five_site-1))
                    sequence1 = fa_dict[chrom][(five_site-6):(five_site+2)]
                    sequence1 = SP.reverse_complement(sequence1)
                    sequence2 = fa_dict[chrom][(three_site-2):(three_site+6)]
                    sequence2 = SP.reverse_complement(sequence2)

                column_dict['seq5'].append(sequence1)
                column_dict['seq3'].append(sequence2)
                
                # Score sequences
                score_5, score_3 = SP.simple_score_junction(sequence1, sequence2, PSSM)
                column_dict['3p score'].append(score_3)
                column_dict['5p score'].append(score_5)
    
    # Create new dataframe from column dictionary
    new_df = pd.DataFrame(columns=column_dict.keys(), index=new_index)
    for column, data in column_dict.iteritems():
        new_df[column] = data
    
    return new_df

def make_quant_df(junc_df, branch_df, gff3, fa_dict, organism=None):
    pssm = SP.generate_consensus_matrix(gff3, fa_dict, PSSM=True)
        
    quant_df = junc_df[(junc_df['type'] != '3prime') & (junc_df['looks like'] != 'AG')]
    
    new_intron_size = []
    alt_splice = []
    score_3 = []
    score_5 = []
    seq5 = []
    seq3 = []

    new_quant_df = pd.DataFrame(index=set(quant_df.index), columns=['intron size','alt splicing'])
    for coord in new_quant_df.index:
        coord_df = quant_df[quant_df.index == coord]

        #Determine if multiple junctions come from this peak
        if len(coord_df) > 1: alt_splice.append(True)
        else: alt_splice.append(False)

        if max(coord_df['annotated intron size']) > 0:
            coord_df = coord_df.sort_values('annotated intron size', ascending=False)
            new_intron_size.append(coord_df.ix[0]['annotated intron size']/1000.)
            seq5.append(coord_df.ix[0]['junction sequence1'])
            seq3.append(coord_df.ix[0]['junction sequence2'])
            score_3.append(coord_df.ix[0]['annotated 3p score'])
            score_5.append(coord_df.ix[0]['annotated 5p score'])
            
        else:
            coord_df = coord_df.sort_values('junction size', ascending=False)
            new_intron_size.append(coord_df.ix[0]['junction size']/1000.)
            seq5.append(coord_df.ix[0]['junction sequence1'])
            seq3.append(coord_df.ix[0]['junction sequence2'])
            scores = SP.simple_score_junction(coord_df.ix[0]['junction sequence1'], coord_df.ix[0]['junction sequence2'], pssm)
            score_3.append(scores[1])
            score_5.append(scores[0])
            
    new_quant_df['intron size'] = new_intron_size
    new_quant_df['alt splicing'] = alt_splice
    new_quant_df['5p score'] = score_5
    new_quant_df['3p score'] = score_3
    new_quant_df['seq5'] = seq5
    new_quant_df['seq3'] = seq3

    quant_df = quant_df.sort_values('annotated intron size')
    quant_df = quant_df.reset_index(drop=True).drop_duplicates(subset='genome coord', keep='first').set_index('genome coord')

    new_quant_df = new_quant_df.merge(quant_df[['transcript','chromosome','position','strand','type']], right_index=True, left_index=True)
    
    for coord in set(branch_df['genome coord']):
        if coord not in new_quant_df.index:
            coord_df = branch_df[branch_df['genome coord'] == coord]
            coord_df = coord_df.sort_values('depth')
            best = coord_df.iloc[0]
            coord_dict = {'transcript':best['transcript'][:-2], 
                         'chromosome':best['chromosome'],
                         'position':best['5p splice site'],
                         'strand':best['strand'],
                         'type':best['type'],
                         'intron size':best['intron size'],
                         'alt splicing':np.where(len(coord_df)> 1, True, False),
                         '5p score':np.NaN,
                         '3p score':np.NaN,
                         'seq5':'','seq3':''}

            
            if len(best['5p seq']) > 0:
                coord_dict['seq5'] = best['5p seq']
            else:
                if best['strand'] == '+':
                    coord_dict['seq5'] = fa_dict[best['chromosome']][(int(best['5p splice site'])-1):(int(best['5p splice site'])+7)]
                elif best['strand'] == '-':
                    coord_dict['seq5'] = fa_dict[best['chromosome']][(int(best['5p splice site'])-6):(int(best['5p splice site'])+2)]
                    coord_dict['seq5'] = SP.reverse_complement(coord_dict['seq5'])
                    
            if str(best['3p splice site']) != 'nan':
                three_site = best['3p splice site']
            else:
                if best['strand'] == '+':
                    after_branch = fa_dict[best['chromosome']][best['branch site']:best['branch site']+100]
                elif best['strand'] == '-':
                    after_branch = fa_dict[best['chromosome']][best['branch site']-100:best['branch site']]
                    after_branch = SP.reverse_complement(after_branch)
                for subs in ['TAG','CAG','GAG','AAG']:
                    if subs in after_branch:
                        ix = after_branch.find(subs)+3
                        break
                three_site = best['branch site']+ix
                if best['strand'] == '-':
                    three_site = best['branch site']-ix
                coord_dict['intron size'] = abs(coord_dict['position']-three_site)
            
            if best['strand'] == '+':
                coord_dict['seq3'] = fa_dict[best['chromosome']][int(three_site-5):int(three_site)+3]
            elif best['strand'] == '-':
                coord_dict['seq3'] = fa_dict[best['chromosome']][int(three_site)-2:int(three_site)+6]
                coord_dict['seq3'] = SP.reverse_complement(coord_dict['seq3'])
                    
            coord_dict['5p score'], coord_dict['3p score'] = SP.simple_score_junction(coord_dict['seq5'], coord_dict['seq3'], pssm)
            coord_s = pd.Series(coord_dict, name=coord)
            new_quant_df = new_quant_df.append(coord_s)
    
    new_quant_df = backfill_splice_sites(new_quant_df, gff3, fa_dict, pssm, organism=organism)
    
    for n in range(len(new_quant_df['seq5'].iloc[0])):     
        new_quant_df['Base 5-'+str(n)] = [x[n] for x in new_quant_df['seq5']]
    for n in range(len(new_quant_df['seq3'].iloc[0])):
        new_quant_df['Base 3-'+str(n)] = [x[n] for x in new_quant_df['seq3']]
    new_quant_df = new_quant_df.drop(['seq5','seq3'], axis=1)
    
    lariat_df = junc_df[(junc_df['type'] == '3prime') | (junc_df['looks like'] == 'AG')]
    lariat_df = lariat_df.sort_values(['genome coord','annotated intron size'], ascending=False)
    lariat_df = lariat_df.reset_index(drop=True).drop_duplicates(subset='genome coord', keep='first').set_index('genome coord')
    lariat_df = lariat_df[['transcript','chromosome','position','strand','type']]
    
    return new_quant_df, lariat_df


def quant_from_peak_df(peak_df, gff3, fa_dict, organism=None):
    count1 = 0
    count2 = 0
    
    pssm = SP.generate_consensus_matrix(gff3, fa_dict, PSSM=True)
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    
    quant_df = peak_df[(peak_df['type'] != '3prime') & (peak_df['looks like'] != 'AG')]
    quant_df['genome coord'] = quant_df['chromosome'].str.cat(quant_df['position'].values.astype(str), sep=':')
    quant_df.index = quant_df['genome coord']
    quant_df = quant_df.drop('index', axis=1)
    
    column_dict = {'intron size':[], 'alt splicing':[], '5p score':[], '3p score':[], 'seq5':[], 'seq3':[]}
    new_index = []
    seq5 = []
    seq3 = []

    for coord in quant_df.index:
        coord_df = quant_df[quant_df.index == coord]
        three_site = None
        alt3 = False
        if len(coord_df) > 0:
            coord_df = coord_df.sort_values('height', ascending=False).ix[0]
        introns = ss_dict[coord_df['transcript']]
        if 'prime' in coord_df['type']:
            peak_range = range(coord_df['position']-5,coord_df['position']+5)
            for intron in introns:
                if intron[0] in peak_range:
                    five_site = intron[0]
                    three_site = intron[1]
                    break
            if len(quant_df[(quant_df['transcript'] == coord_df['transcript']) & (quant_df['type'] == 'AG')]) > 0:
                alt3=True
        else:
            if 'AG' in quant_df[quant_df['transcript'] == coord_df['transcript']]['type']:
                five_site = coord_df['position']
                three_df = quant_df[(quant_df['transcript'] == coord_df['transcript']) & (quant_df['type'] == 'AG')]
                three_df = three_df.sort_values('height', ascending=False)
                three_site = three_df.ix[0]['position']
        
        if three_site is not None:
            new_index.append(coord)
            size = abs(three_site-five_site)/1000.
            column_dict['intron size'].append(size)
            column_dict['alt splicing'].append(alt3)
            
            if coord_df['strand'] == '+':
                s5 = fa_dict[coord_df['chromosome']][five_site-2:five_site+6]
                s3 = fa_dict[coord_df['chromosome']][three_site-6:three_site+2]
            elif coord_df['strand'] == '-':
                s5 = fa_dict[coord_df['chromosome']][five_site-6:five_site+2]
                s5 = SP.reverse_complement(s5)
                s3 = fa_dict[coord_df['chromosome']][three_site-2:three_site+6]
                s3 = SP.reverse_complement(s3)
            column_dict['seq5'].append(s5)
            column_dict['seq3'].append(s3)
            scores = SP.simple_score_junction(s5, s3, pssm)
            column_dict['3p score'].append(scores[1])
            column_dict['5p score'].append(scores[0])
            
    new_quant_df = quant_df[quant_df.index.isin(new_index)][['genome coord','chromosome',
                                                             'strand','transcript','position','type']]
    for column, data in column_dict.iteritems():
        new_quant_df[column] = data
    
    new_quant_df = new_quant_df.drop_duplicates(subset='genome coord', keep='first').set_index('genome coord')
    
    new_quant_df = SP.backfill_splice_sites(new_quant_df, gff3, fa_dict, pssm, organism=organism)
    
    #for n in range(len(new_quant_df['seq5'].iloc[0])):     
    #    new_quant_df['Base 5-'+str(n)] = [x[n] for x in new_quant_df['seq5']]
    #for n in range(len(new_quant_df['seq3'].iloc[0])):
    #    new_quant_df['Base 3-'+str(n)] = [x[n] for x in new_quant_df['seq3']]
    #new_quant_df = new_quant_df.drop(['seq5','seq3'], axis=1)
    
    new_quant_df = SP.find_score_branches_ppy(new_quant_df, '/home/jordan/GENOMES/S288C/S288C_branches2.txt', fa_dict)
    
    return new_quant_df


def check_intron_position(transcript, position, gff3, organism):
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    
    first=False
    last=False
    
    introns = ss_dict[transcript]
    
    for n, intron in enumerate(introns):
        if intron[0] in range(position-3,position+3):
            if n == 0:
                first = True
            elif n == len(intron):
                last = True
            break
    return first, last


def count_reads_at_peaks(bam_files, genome_coords, organism=None):
    bams = {}
    for bam_file in bam_files:
        bams[bam_file] = pysam.Samfile(bam_file)
    
    all_reads = {}
    for bam, reader in bams.iteritems():
        print bam
        all_reads[bam] = []

        for coord, strand in genome_coords:
            peak_count = 0
            pos = int(coord.split(':')[1])
            chrom = coord.split(':')[0]
            if organism == 'pombe':
                lat_rom = {'chr1':'I','chr2':'II','chr3':'III'}
                chrom = lat_rom[chrom]
                
            peak_iter = reader.fetch(chrom,  pos-100,  pos+100)
            
            for read in peak_iter:
                if read.is_reverse and strand == '+':
                    if read.reference_end in range(pos-3,pos+3):
                        peak_count += 1
                elif not read.is_reverse and strand == '-':
                    if read.reference_start in range(pos-3,pos+3):
                        peak_count += 1
            
            all_reads[bam].append(peak_count)
    
    return all_reads

def count_reads_in_transcript(bam_files, df, gff3, organism=None):
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    
    bams = {}
    for bam_file in bam_files:
        bams[bam_file] = pysam.Samfile(bam_file)
    
    all_reads = {}

    for bam, reader in bams.iteritems():
        all_reads[bam] = pd.DataFrame(index=df.index, columns=['total','intron'])
        
        for tx in set(df['transcript']):
            tx_df = df[df['transcript'] == tx]
            if organism == 'pombe':
                tx = tx+'.1'
            else:
                tx = tx+'T0'
                
            start, end, strand, CDS_start, CDS_end, exons, chrom = SP.tx_info(tx, tx_dict)
            if organism == 'pombe':
                lat_rom = {'chr1':'I','chr2':'II','chr3':'III'}
                chrom = lat_rom[chrom]
            
            tx_iter = reader.fetch(chrom,  start,  end)
            
            intron_ranges = {}
            for ix, r in tx_df.iterrows():
                if strand == '+':
                    intron_start = int(r['position'])
                    intron_end = int(r['position']+r['intron size'])+1
                elif strand == '-':
                    intron_start = int(r['position']-r['intron size'])
                    intron_end = int(r['position'])+1
                intron_ranges[ix] = [range(intron_start,intron_end),0]

            reads = 0
            for read in tx_iter:
                if read.is_reverse and strand == '+':
                    reads += 1
                    
                    for ix in intron_ranges:
                        if read.reference_end in intron_ranges[ix][0]:
                            intron_ranges[ix][1] += 1
                    
                elif not read.is_reverse and strand == '-':
                    reads += 1
                    for ix in intron_ranges:
                        if read.reference_start in intron_ranges[ix][0]:
                            intron_ranges[ix][1] += 1
                            
            for ix in intron_ranges:
                try:
                    all_reads[bam].loc[ix,'total'] = reads/float(end-start)*1000
                    all_reads[bam].loc[ix,'intron'] = ((intron_ranges[ix][1]/float(tx_df.loc[ix,'intron size'])) /
                                                   (reads/float(end-start)))
                except ZeroDivisionError:
                    all_reads[bam].loc[ix,'total'] = np.NaN
                    all_reads[bam].loc[ix,'intron'] = np.NaN
                    print ix
                    
    return all_reads 

def count_precursor(bam_files, genome_coords, organism=None):
    bams = {}
    for bam_file in bam_files:
        bams[bam_file] = pysam.Samfile(bam_file)
    
    all_reads = {}
    for bam, reader in bams.iteritems():
        print bam
        all_reads[bam] = []

        for coord, strand in genome_coords:
            peak_count = 0
            pos = int(coord.split(':')[1])
            chrom = coord.split(':')[0]
            if organism == 'pombe':
                lat_rom = {'chr1':'I','chr2':'II','chr3':'III'}
                chrom = lat_rom[chrom]
                
            peak_iter = reader.fetch(chrom,  pos-100,  pos+100)
            
            for read in peak_iter:
                if read.reference_end in range(pos+2,pos+55) and read.reference_start in range(pos-55,pos-2):
                    if read.is_reverse and strand == '+':
                        peak_count += 1
                    elif not read.is_reverse and strand == '-':
                        peak_count += 1
            
            all_reads[bam].append(peak_count/0.053)
    
    return all_reads

def quantitate_junction_df(bam_dict, df, gff3, W=True, organism=None):
    quant_df = copy.deepcopy(df)
    quant_df['genome coord'] = quant_df.index
    quant_df = quant_df.reset_index(drop=True).drop_duplicates(subset='genome coord', keep='first').set_index('genome coord')
    
    A_reads = count_reads_at_peaks([bam_dict['A1'], bam_dict['A2']], zip(quant_df.index.tolist(), quant_df['strand'].tolist()),
                                  organism=organism)
    B_reads = count_reads_in_transcript([bam_dict['B1'], bam_dict['B2']], quant_df, gff3, organism=organism)
    
    # read_dicts is a list of dictionaries where each dictionary has bam file names for keys and a 
    # list of read values in the order of the dataframe index
    read_dicts = [A_reads, B_reads]
    
    if W is True:
        W_reads = count_reads_in_transcript([bam_dict['W1'], bam_dict['W2']], quant_df, gff3, organism=organism)
        read_dicts.append(W_reads)
    
    precursor = count_precursor([bam_dict['B1'], bam_dict['B2']], zip(quant_df.index.tolist(), quant_df['strand'].tolist()),
                               organism=organism)
    
    names = []
    for read_dict in read_dicts:
        for bam, read_counts in read_dict.iteritems():
            name = bam.split('/')[-1].split('_sorted')[0]
            names.append(name)
            if name in quant_df:
                pass
            
            if type(read_counts) == list:
                quant_df[name] = read_counts
            else:
                read_counts[name] = read_counts['total']
                read_counts[name+' Intron Retention'] = read_counts['intron']
                quant_df = quant_df.merge(read_counts[[name, name+' Intron Retention']], right_index=True, left_index=True)
                quant_df[name+' Intron Retention'] = pd.to_numeric(quant_df[name+' Intron Retention'])

    for bam, read_counts in precursor.iteritems():
        name = bam.split('/')[-1].split('_sorted')[0]+' near peak'
        names.append(name)
        quant_df[name] = read_counts
        
    quant_df = quant_df.dropna(how='any')
    for name in names:
        if name.endswith('A'):
            sample = name.split('-A')[0]
            quant_df[sample+' Intermediate Level'] = quant_df[name].divide(sum(quant_df[name])/1000000.)/quant_df[sample+'-B'].divide(sum(quant_df[sample+'-B'])/1000000.)
            quant_df[sample+' Intermediate Level'] = pd.to_numeric(quant_df[sample+' Intermediate Level'])
        
        if name.endswith('B'):
            sample = name.split('-B')[0]
            quant_df[sample+' Total Spliceosomal'] = quant_df[name].divide(sum(quant_df[name])/1000000.)
            quant_df[sample+' Total Spliceosomal'] = pd.to_numeric(quant_df[sample+' Total Spliceosomal'])
            
            quant_df[sample+' Precursor'] = quant_df[name+' near peak'].divide(sum(quant_df[name+' near peak'])/1000000.)/quant_df[name].divide(sum(quant_df[name])/1000000.)
            quant_df[sample+' Precursor'] = pd.to_numeric(quant_df[sample+' Precursor'])
                                             
        if name.endswith('W'):
            sample = name.split('-W')[0]
            quant_df[sample+' RNAseq'] = quant_df[name].divide(sum(quant_df[name])/1000000.)
            quant_df[sample+' RNAseq'] = pd.to_numeric(quant_df[sample+' RNAseq'])
    
    if W is True: metrics = ['Intermediate Level','Total Spliceosomal','Precursor','B Intron Retention','RNAseq','W Intron Retention']
    else: metrics = ['Intermediate Level','Total Spliceosomal','Precursor','B Intron Retention']
    for metric in metrics:
        cols = [x for x in quant_df.columns if metric in x]
        quant_df[metric+' avg'] = quant_df[cols].sum(axis=1).divide(float(len(cols))).apply(np.log2)
        quant_df[metric+' avg'] = quant_df[metric+' avg'].replace([np.inf, np.inf*-1], np.NaN)
        
    quant_df = quant_df.dropna(how='any')
    
    return quant_df

### Tools for quantitating new samples

def config_mut_quant(config_file):
    bam_dict = {}
    with open(config_file, 'r') as config:
        for line in config:
            info = line.split(',')
            genotype = info[1]
            sample = info[2].strip()
            
            if genotype not in bam_dict:
                bam_dict[genotype] = {}
            
            bam_dict[genotype][sample] = info[0]
    return bam_dict

def main():
    '''Each line will be : bam_file,genotype,sample
    e.g. CM763-A_sorted.bam,WT,A1'''
    
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
               'chromosome','position','alt splicing','3p score','transcript','intron position','strand','peak',
               'Base 5-0','Base 5-1','Base 5-2','Base 5-3','Base 5-4','Base 5-5','Base 5-6','Base 5-7','Base 3-0',
               'Base 3-1','Base 3-2','Base 3-3','Base 3-4','Base 3-5','Base 3-6','Base 3-7','branch score',
               'branch to 3p distance','percent pPy','branch-0','branch-1','branch-2','branch-3','branch-4']
    
    quant_df = pd.read_csv(sys.argv[2], index_col=0)
    try:
        quant_df = quant_df[columns]
    except KeyError:
        print "Columns missing from dataframe..."
        print columns
        print quant_df.columns
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
        new_df = quantitate_junction_df(samples, quant_df, gff3, W=W)
        
        # Remove original columns and rename new ones with multiindex
        new_columns = [x for x in new_df.columns if x not in columns]
        new_df = new_df[new_columns]
        new_df.columns = pd.MultiIndex.from_product([[genotype], new_df.columns])
        final_df = final_df.join(new_df, how='inner')
        #final_df = final_df.merge(new_df, right_index=True, left_index=True)
        
    final_df.to_csv(prefix+'_quant_df.csv')
    final_df.to_pickle(prefix+'_quant_df.pickle')
        
    SP.SP_quant_scatters(final_df.dropna(how='any'), bam_dict, W=W)
    
def s_log2(s):
    s = pd.to_numeric(s)
    s = s.apply(np.log2)
    s = s.replace([np.inf,-1*np.inf],np.NaN)
    return s

def s_log2_ratio_Zscore(s1,s2):
    s = s2/s1
    s = s_log2(s)
    s = s.replace([np.inf,-1*np.inf],np.NaN).dropna()
    pd.Series(stats.mstats.zscore(s), index=s.index)
    return s   

def log2_Zscore_df(df, wt, mut, metrics=['Intermediate Level', 'Precursor'], Z=2, by_pos_scores=False):
    if type(df) == str:
        new_df = pd.read_pickle(df)
    else:
        new_df = copy.deepcopy(df)
        
    print len(new_df)
    mutA = [x for x in new_df.columns if (x[0] == mut) and (x[1][-2:] == '-A')]
    #wtA = [x for x in new_df.columns if (x[0] == wt) and (x[1][-2:] == '-A')]
    new_df = new_df[new_df[mutA].sum(axis=1) >= 10]
    print len(new_df)
    
    for metric in metrics:
        columns = [x for x in new_df.columns if (metric in x[1]) and ('avg' not in x[1])]
        wt_cols = [x for x in columns if (x[0] == wt) and ('avg' not in x[1])]
        mut_cols = [x for x in columns if (x[0] == mut) and ('avg' not in x[1])]
        
        for column in columns:
            new_df[(column[0], column[1]+' log2')] = s_log2(new_df[column])
            
        if len(wt_cols) != len(mut_cols):
            print "Number of WT reps must match number of mutant reps!"
            print wt_cols
            print mut_cols
            return None
        
        for n, wt_col in enumerate(wt_cols):
            new_df[('All',metric+' log2 ratio'+str(n+1))] = s_log2(new_df[mut_cols[n]]/new_df[wt_col])
            new_index = [x+'-'+str(n) for x in new_df.index]
            
            wt_s = new_df[wt_col]
            wt_s.index = new_index
            mut_s = new_df[mut_cols[n]]
            mut_s.index = new_index
            
            if n == 0:
                wt_s_for_Z = wt_s
                mut_s_for_Z = mut_s
            else:
                wt_s_for_Z = wt_s_for_Z.append(wt_s)
                mut_s_for_Z = mut_s_for_Z.append(mut_s)
            
        Zlist = s_log2_ratio_Zscore(wt_s_for_Z.dropna(), mut_s_for_Z.dropna())
        
        for n, wt_col in enumerate(wt_cols):
            n_up = Zlist[(Zlist.index.str[-1] == str(n)) & (Zlist >= Z)]
            n_up = [x[:-2] for x in n_up.index]
            
            n_down = Zlist[(Zlist.index.str[-1] == str(n)) & (Zlist <= -1*Z)]
            n_down = [x[:-2] for x in n_down.index]
            
            n_other = Zlist[(Zlist.index.str[-1] == str(n)) & (Zlist < Z) & (Zlist > -1*Z)]
            n_other = [x[:-2] for x in n_other.index]
            
            if n == 0:
                up = set(n_up)
                down = set(n_down)
                other = set(n_other)
            
            else:
                up = up.intersection(n_up)
                up = up.difference(n_down).difference(n_other)
                down = down.intersection(n_down)
                down = down.difference(n_up).difference(n_other)
                other = other.intersection(n_other)
                other = other.difference(n_up).difference(n_down)
        
            print len(up)
            print len(down)
        
        new_df[('All',metric+' change')] = None
        new_df.loc[up, ('All',metric+' change')] = 'Up'
        new_df.loc[down, ('All',metric+' change')] = 'Down'
        new_df.loc[other, ('All', metric+' change')] = 'Other'
        
        plot_df = copy.deepcopy(new_df)
        
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8))
        groups = {'Other':'0.8','Up':'tomato','Down':'cornflowerblue'}
        for group in ['Other','Up','Down']:
            gr_df = plot_df[plot_df[('All',metric+' change')] == group]
            if len(gr_df) >= 15:
                for n, wt_col in enumerate(wt_cols):
                    ax[0][n].scatter(s_log2(gr_df[wt_col]), s_log2(gr_df[mut_cols[n]]), 
                                color=groups[group], alpha=0.9, label=group, s=20)

                    ax[0][n].set_xlabel(wt_col[0]+' log2 '+metric, fontsize=12)
                    ax[0][n].set_ylabel(mut_cols[n][0]+' log2 '+metric, fontsize=12)
                    ax[0][n].set_title('Replicate '+str(n+1), fontsize=14)

                    ax[0][n], limits = SP.draw_diagonal(ax[0][n])
                    ax[0][n].legend(fontsize=12)

                sns.kdeplot(gr_df[('Peaks','intron size')], ax=ax[1][0], bw=2, cumulative=True, linewidth=3, 
                            color=groups[group], label=group)
                ax[1][0].set_xlim([30, 400])

                sns.kdeplot(gr_df[('Peaks','5p score')], ax=ax[1][1], bw=2, cumulative=True, linewidth=3, 
                            color=groups[group], label=group)

                ax[1][0].set_xlabel('Intron size (nt)')
                ax[1][0].set_ylabel('Fraction of introns')
                ax[1][1].set_xlabel('5prime splice site score')
                ax[1][1].set_ylabel('Fraction of introns')
        
        ax[1][1].set_xlim([np.percentile(plot_df[('Peaks','5p score')], 0.5),
                              np.percentile(plot_df[('Peaks','5p score')], 99.9)+5])
            
        fig.tight_layout()
        plt.show()
        plt.clf()

    if by_pos_scores is True:
        SP.by_pos_plots(new_df, metrics=metrics)
    
    new_df[('Peaks','predicted')] = True
    new_df.loc[~new_df[('Peaks','type')].str.contains('prime'), ('Peaks','predicted')] = False    
    return new_df

if __name__ == "__main__":
    main()
        
        