import sys
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')

import GeneUtility
import SPTools as SP

import numpy as np
import pandas as pd
import json
import HTSeq
import itertools

from subprocess import call
from multiprocessing import Pool
from functools import partial
import os
from itertools import izip,islice,tee
import re

### Collect all 5p splice sites of interest and get 20 nt of sequence before splice site - pass in splice site dict
def collect_intron_seq(gff3_file, fasta_file, ss_dict=None, junction_bed=None, gene_list=None, peak_df=None, organism=None):
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta_file) == dict:
        fasta_dict = fasta_file
    elif fasta_file.endswith('json'):
        with open(fasta_file, 'r') as f:
            fasta_dict = json.load(f)
    else:
        fasta_dict = make_fasta_dict(fasta_file)
    if ss_dict is not None:
        ss_dict=ss_dict
    elif junction_bed is not None:
        ss_dict = SP.build_junction_dict(junction_bed, gff3_file, transcript_dict, organism=organism)
    elif peak_df is not None:
        ss_dict = {}
        peak_df = peak_df[~peak_df['type'].str.contains('prime')]
        for ix, r in peak_df.iterrows():
            if r['transcript'] not in ss_dict:
                ss_dict[r['transcript']] = []
            if r['strand'] == '+':
                ss_dict[r['transcript']].append((r['position'],r['position']+50))
            elif r['strand'] == '-':
                ss_dict[r['transcript']].append((r['position'],r['position']-50))
                
    else:
        ss_dict, intron_flag = SP.list_splice_sites(gff3_file, gene_list=gene_list, organism=organism)
        ss_dict = SP.collapse_ss_dict(ss_dict)
    
    seq_dict = {}
    for transcript, introns in ss_dict.iteritems():
        if junction_bed is None:
            if organism == 'pombe':
                transcript = transcript+'.1'
            else:
                transcript = transcript+'T0'
        introns = list(introns)
        strand = transcript_dict[transcript][2]
        chrom = transcript_dict[transcript][3]
        n = 0
        for n in range(len(introns)):
            if strand == '+':
                seq_dict[transcript+'-'+chrom+':'+str(introns[n][0]+1)] = fasta_dict[chrom][introns[n][0]+2:introns[n][0]+17]
            elif strand == '-':
                seq = fasta_dict[chrom][introns[n][0]-16:introns[n][0]-1]
                seq_dict[transcript+'-'+chrom+':'+str(introns[n][0])] = SP.reverse_complement(seq)
    return seq_dict

### Samtools view and grep for each sequence
def sub_findre(s,substring,diffnumber):
    sublen=len(substring)
    zip_gen=(izip(substring,islice(s,i,i+sublen)) for i in xrange(len(s)))
    for z in zip_gen:
        l,z=tee(z)
        if sum(1 for i,j in l if i==j)>=sublen-diffnumber:
            new=izip(*z)
            next(new)
            yield ''.join(next(new))
            
def sub_findre2(s,substring):
    if substring not in sub_findre2.cache:
        sub_findre2.cache[substring] = [
            re.compile('.'.join((substring[:x], substring[x+1:])))
            for x in range(len(substring))]
    for regex in sub_findre2.cache[substring]:
        for match in regex.findall(s):
            yield match
sub_findre2.cache = {}

def a_rich(string):
    a_count = 0
    for char in string:
        if char == 'A':
            a_count += 1
    perc_A = float(a_count)/len(string)
    if perc_A >= 0.4:
        high_a = True
    else:
        high_a = False
    return high_a
            
def process_reads((fastq_file, seq_dict), prefix=None, mismatches=0):
    if prefix is None:
        prefix = fastq_file
    fa_line_list = []
    unsplit_fa_line_list = []
    
    if '.fastq' in fastq_file or '.fq0' in fastq_file or '.fq1' in fastq_file or fastq_file.endswith('fq'):
        fastq = HTSeq.FastqReader(fastq_file, "solexa")
    elif '.fasta' in fastq_file or '.fa0' in fastq_file or '.fa1' in fastq_file or fastq_file.endswith('fa'):
        fastq = HTSeq.FastaReader(fastq_file)
    
    if '_unsplit' in fastq_file:
        mismatches = 1
    
    for read in fastq:
        new_read = None
        for intron, seq in seq_dict.iteritems():
            if seq in read.seq:
                read_seq = read.seq.split(seq)[0]
                new_length = len(read_seq)
                new_read = [read.name+':'+intron, read_seq]
                break
            elif seq not in read.seq and mismatches > 0:
                in_read = list(sub_findre2(read.seq, seq))
                if len(in_read) > 0:
                    for n, subseq in enumerate(in_read):
                        read_seq = read.seq.split(subseq)[0]
                        new_length = len(read_seq)
                        new_read = [read.name+':'+intron+'-'+str(n), read_seq]
                    break
        
        if new_read is not None and len(new_read[1]) > 15:
            fa_line_list.append(('>'+new_read[0]+'\n',new_read[1]+'\n'))
        
        if new_read is None and a_rich(read.seq) is False:
            unsplit_fa_line_list.append(('>'+read.name+'\n',read.seq+'\n'))
        
        if len(fa_line_list) == 5000:
            with open('{}_split.fa'.format(prefix), 'a') as fout:
                for pair in fa_line_list:
                    fout.write(pair[0])
                    fout.write(pair[1])
            fa_line_list = []
        
        if len(unsplit_fa_line_list) == 5000:
            with open('{}_unsplit.fa'.format(prefix), 'a') as fout:
                for pair in unsplit_fa_line_list:
                    fout.write(pair[0])
                    fout.write(pair[1])
            unsplit_fa_line_list = []            
            
    with open('{}_split.fa'.format(prefix), 'a') as fout:
        for pair in fa_line_list:
            fout.write(pair[0])
            fout.write(pair[1])
            
    with open('{}_unsplit.fa'.format(prefix), 'a') as fout:
        for pair in unsplit_fa_line_list:
            fout.write(pair[0])
            fout.write(pair[1]) 
            
    with open('fa_list.txt','a') as f:
        f.write('{}_split.fa\n'.format(prefix))
        
    with open('unsplit_fa_list.txt','a') as f:
        f.write('{}_unsplit.fa\n'.format(prefix))
        

def find_split_reads(fastq_file, seq_dict, prefix, threads=0):
    if threads > 0:
        p = Pool(threads)
        fastq_list = split_fastq_file(fastq_file, threads)
        params = []
        for fastq in fastq_list:
            params.append((fastq, seq_dict))
        p.map(process_reads, params)
        p.close()
        p.join()
        
        fa_list = []
        with open('fa_list.txt','r') as f:
            for line in f:
                fa_list.append(line.strip())

        unsplit_fa_list = []
        with open('unsplit_fa_list.txt','r') as f:
            for line in f:
                unsplit_fa_list.append(line.strip())                

        call_list = ["cat",fastq_file+"*_split.fa",">","{0}_split.fa".format(prefix)]
        script = ' '.join(call_list)
        call(script, shell=True)
        
        call_list = ["cat",fastq_file+"*_unsplit.fa",">","{0}_unsplit.fa".format(prefix)]
        script = ' '.join(call_list)
        call(script, shell=True)
        
        for fastq in fastq_list:
            os.remove(fastq)
        try:
            for fa in fa_list:
                os.remove(fa)

            for fa in unsplit_fa_list:
                os.remove(fa)
        except OSError:
            pass
            
        os.remove('fa_list.txt')
        os.remove('unsplit_fa_list.txt')
        
    else:
        process_reads((fastq_file, seq_dict), prefix=prefix)
    
def split_fastq_file(fastq_file, threads):
    with open(fastq_file, 'r') as f:
        for i, l in enumerate(f):
            pass
    file_len = i + 1
    num_files = file_len/20000
    
    call(["split", "-d", "-l 20000" , "-a 5", fastq_file, fastq_file])
    
    n=0
    fastq_list = []
    for n in range(num_files):
        fastq_list.append(fastq_file+format(n, '05'))
    return fastq_list
        
    
def list_branch_points(sorted_bam_file, gff3_file, fasta_dict, organism=None):
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)

    if type(fasta_dict) == str:
        with open(fasta_dict, 'r') as f:
            fasta_dict = json.load(f)
    
    branch_dict = {}
    read_counter = 0
    br_counter = 0
    bam_reader = HTSeq.BAM_Reader(sorted_bam_file)
    for a in bam_reader:
        read_counter += 1
        transcript = a.read.name.split('-chr')[0].split(':')[-1]
        splice_site = a.read.name.split('-')[-1]
        if len(splice_site) < 3:
            splice_site = a.read.name.split('-')[-2]
        if splice_site.startswith('chr'):
            if transcript not in branch_dict:
                branch_dict[transcript] = {}
            if splice_site not in branch_dict[transcript]:
                branch_dict[transcript][splice_site] = []
            if a.iv is not None:
                strand = a.iv.strand
                read_end = a.iv.end
                if strand == '-':
                    read_end = a.iv.start
                if strand == transcript_dict[transcript][2]:
                    branch_dict[transcript][splice_site].append(read_end)
                    br_counter += 1
                
    print "Reads analyzed: "+str(read_counter)
    print "Reads assigned as branches: "+str(br_counter)
    
    new_branch_dict = {}
    for transcript, introns in branch_dict.iteritems():
        new_branch_dict[transcript] = []
        for intron, branches in introns.iteritems():
            new_branch_list = []
            new_branch_counts = []
            for branch in branches:
                flag = False
                if len(new_branch_list) > 0:
                    for pos in range(branch-2,branch+3):
                        if pos in new_branch_list: 
                            flag = True
                            br_id = new_branch_list.index(pos)
                            new_branch_counts[br_id] += 1
                if flag == False: 
                    new_branch_list.append(branch)
                    new_branch_counts.append(1)
            if len(new_branch_list) > 0:
                new_branch_dict[transcript].append([intron, new_branch_list, new_branch_counts])
    
    with open('{0}.bed'.format(sorted_bam_file.split('_sorted.bam')[0]), 'w') as fout:
        fout.write('track name=junctions description="TopHat junctions"\n')
        for transcript, introns in new_branch_dict.iteritems():
            strand = transcript_dict[transcript][2]
            for intron in introns:
                chrom = intron[0].split(':')[0]
                start = int(intron[0].split(':')[1])
                n=0
                for n in range(len(intron[1])):
                    end = intron[1][n]
                    value = intron[2][n]
                    size = abs(end-start)+30
                    if abs(end-start) > 2000:
                        pass
                    elif abs(end-start) > 5 and value >= 5:
                        #[seqname] [start] [end] [id] [score] [strand] [thickStart] [thickEnd] [r,g,b][block_count] [block_sizes] [block_locations]
                        read_id = intron[0]+'-'+str(n)
                        block_size = '0,'+str(size)
                        line_list = [chrom, str(start-1), str(end+1), read_id, str(value), strand, str(start-1), str(end+1), '75,196,213', '2', '1,1', block_size, '\n']
                        line = '\t'.join(line_list)
                        fout.write(line)
    
    return new_branch_dict
            
###################################
##  Script for finding branches  ##
###################################
    
def main():
    '''Usage: run SPBranch.py unmapped1 unmapped2 threads organism [config_file] [untagged]
    
    Parameters
    -----------
    unmapped1 : bam or fastq file of unmapped reads from tophat or bowtie
    unmapped2 : bam or fastq file of unmapped reads from tophat or bowtie
    threads : number of processors to use
    organism : 'pombe or 'crypto'
    config_file : if using peaks to call - list of changepoint output file names and where to find them
    untagged : untagged sample name (must be in file name)
    
    Output
    ------
    bam files with aligned reads. Will be interpreted by SP_pipeline.
    '''
    
    unmapped1 = sys.argv[1]
    unmapped2 = sys.argv[2]
    threads = int(sys.argv[3])
    
    if unmapped1.endswith('bam'):
        btf_args = 'bamToFastq -i {0} -fq {1}'.format(unmapped1, unmapped1.split('.bam')[-1]+'.fq')
        call(btf_args, shell=False)
        unmapped1 = unmapped1.split('.bam')[-1]+'.fq'
    if unmapped2.endswith('bam'):
        btf_args = 'bamToFastq -i {0} -fq {1}'.format(unmapped2, unmapped2.split('.bam')[-1]+'.fq')
        call(btf_args, shell=False)
        unmapped2 = unmapped2.split('.bam')[-1]+'.fq'
        
    cat_args = 'cat {0} {1} > unmapped_all.fq'.format(unmapped1, unmapped2)
    call(cat_args, shell=True)
    
    organism = sys.argv[4]
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
        
    peaks = False
    if len(sys.argv) == 7:
        peaks = True
        with open(sys.argv[5], 'r') as config:
            for line in config:
                if sys.argv[6] in line:
                    CP_untagged = line.strip()
                elif 'changepoint' in line.lower() or 'peak' in line.lower():
                    CP_out.append(line.strip())
        peak_df = SP.peak_to_seq_pipeline(CP_untagged, CP_out[0], CP_out[1], gff3, fa_dict, name='CP_peaks')

    ann_seqs = collect_intron_seq(gff3, fa_dict)
    
    print "Finding unaligned reads with annotated 5' splice sites"
    find_split_reads('unmapped_all.fq', ann_seqs, 'Ann_branches', threads=threads)
    
    print "Aligning split reads to the genome with Bowtie"
    bowtie_args = 'bowtie -p{0} -v1 -M1 --best {1} -f Ann_branches_split.fa --sam Ann_branches.sam'.format(threads, bowtie_index)
    call(bowtie_args, shell=True)
    
    # sort and index
    print "Sorting and indexing bam file"
    samtools1 = 'samtools view -Sbo Ann_branches.bam Ann_branches.sam'
    call(samtools1, shell=True)
    
    samtools2 = 'samtools sort Ann_branches.bam -o Ann_branches_sorted.bam'
    call(samtools2, shell=True)
    
    samtools3 = 'samtools index Ann_branches_sorted.bam'
    call(samtools3, shell=True)
    
    if peaks is True:
        print "Finding unaligned reads with unpredicted splicing events"
        peak_seqs = collect_intron_seq(gff3, fa_dict, peak_df=peak_df)
        find_split_reads('Ann_branches_unsplit.fa', peak_seqs, 'Peak_branches', threads=threads)
        
        print "Aligning split reads to the genome with Bowtie"
        bowtie_args = 'bowtie -p{0} -v1 -M1 --best {1} -f Peak_branches_split.fa --sam Peak_branches.sam'.format(threads, bowtie_index)
        call(bowtie_args, shell=True)
        
        print "Sorting and indexing bam file"
        samtools1 = 'samtools view -Sbo Peak_branches.bam Peak_branches.sam'
        call(samtools1, shell=True)

        samtools2 = 'samtools sort Peak_branches.bam -o Peak_branches_sorted.bam'
        call(samtools2, shell=True)

        samtools3 = 'samtools index Peak_branches_sorted.bam'
        call(samtools3, shell=True)
    
#############################################################################
## Tools for creating a dataframe with branches found by above method      ##
#############################################################################

def create_branch_df(branch_dict, gff3, fa_dict, organism=None):
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    chroms = []
    fives = []
    transcripts = []
    branches = []
    depths = []
    strands = []
    distances = []
    for tx, five_sites in branch_dict.iteritems():
        for five_site in five_sites:
            chrom = five_site[0].split(':')[0]
            pos = int(five_site[0].split(':')[1])
            n=0
            for n in range(len(five_site[1])):
                if abs(five_site[1][n]-pos) > 5 and abs(five_site[1][n]-pos) <= 1000 and five_site[2][n] >= 5:
                    chroms.append(chrom)
                    fives.append(pos)
                    transcripts.append(tx)
                    branches.append(five_site[1][n])
                    depths.append(five_site[2][n])
                    strands.append(tx_dict[tx][2])
                    if tx_dict[tx][2] == '+':
                        distances.append(five_site[1][n]-pos)
                    elif tx_dict[tx][2] == '-':
                        distances.append(pos-five_site[1][n])
    branch_df = pd.DataFrame(index = range(len(fives)))
    branch_df['transcript'] = transcripts
    branch_df['chromosome'] = chroms
    branch_df['5p splice site'] = fives
    branch_df['branch site'] = branches
    branch_df['depth'] = depths
    branch_df['distance'] = distances
    branch_df['strand'] = strands
    
    branch_df = branch_df[branch_df['distance'] > 0]
    branch_df['genome coord'] = branch_df['chromosome'].str.cat(branch_df['5p splice site'].apply(int).apply(str), sep=':')
    branch_df['branch coord'] = branch_df['chromosome'].str.cat(branch_df['branch site'].apply(int).apply(str), sep=':')
    
    branch_df = add_seq(branch_df, fa_dict)
    branch_df = find_3p_site(branch_df, gff3, organism=organism)
    return branch_df

def add_seq(branch_df, fa_dict):
    five_seqs = []
    branch_seqs = []
    for ix, r in branch_df.iterrows():
        five = fa_dict[r['chromosome']][r['5p splice site']-8:r['5p splice site']+8]
        branch = fa_dict[r['chromosome']][r['branch site']-8:r['branch site']+8]
        if r['strand'] == '-':
            five = SP.reverse_complement(five)
            branch = SP.reverse_complement(branch)
        if 'GT' in five[4:11]:
            ix = five.index('GT')
            five = five[ix-2:ix+6]
        else:
            five = five[4:12]
        if 'AG' in branch[4:11]:
            ix = branch.index('AG')
            branch = branch[ix-4:ix+4]
        elif 'AA' in branch[4:11]:
            ix = branch.index('AA')
            branch = branch[ix-4:ix+4]
        elif 'GA' in branch[4:11]:
            ix = branch.index('GA')
            branch = branch[ix-4:ix+4]
        else:
            branch = branch[4:13]
        five_seqs.append(five)
        branch_seqs.append(branch)
    branch_df['5p seq'] = five_seqs
    branch_df['Branch seq'] = branch_seqs
    
    receptors = ['AG', 'AA', 'GA']
    branch_df = branch_df[branch_df['Branch seq'].str[4:6].isin(receptors)]
    return branch_df

def list_alt_branch(branch_df, organism=None):
    alt_branches = []
    for five_site in set(branch_df['genome coord']):
        five_df = branch_df[branch_df['genome coord'] == five_site]
        if len(five_df) > 1:
            alt_branches.append(five_site)
    print len(alt_branches)
    return alt_branches

def find_3p_site(branch_df, gff3, organism=None):
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    
    three_coord = []
    for ix, r in branch_df.iterrows():
        introns = ss_dict[r['transcript'][:-2]]
        matched = False
        for intron in introns:
            if r['5p splice site'] in range(intron[0]-1,intron[0]+2):
                three_coord.append(intron[1])
                matched = True
                break
        if matched is False:
            three_coord.append(np.NaN)
    
    branch_df['3p splice site'] = three_coord
    branch_df['intron size'] = branch_df['5p splice site']-branch_df['3p splice site']
    branch_df['intron size'] = branch_df['intron size'].apply(abs)
    branch_df['Branch to 3p distance'] = branch_df['branch site']-branch_df['3p splice site']
    branch_df['Branch to 3p distance'] = branch_df['Branch to 3p distance'].apply(abs)
    
    return branch_df

def write_seq_list_to_file(df, prefix):
    with open('{0}_5pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['5p sequence'].tolist():
            f.write(seq+'\n')
    with open('{0}_3pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['3p sequence'].tolist():
            f.write(seq+'\n')
    with open('{0}_BPseq_list.txt'.format(prefix),'w') as f:
        for seq in df['Branch sequence'].tolist():
            f.write(seq+'\n')

if __name__ == "__main__":
    main()