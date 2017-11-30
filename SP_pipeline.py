import sys
import os
script_path = os.path.dirname(os.path.realpath(__file__)).split('SPTools')[0]
sys.path.append(script_path)
import SPTools as SP
import json
import numpy as np
import pandas as pd
from subprocess import call

def peak_junction_analysis(peak_df, junc_beds, gff3, fa_dict, organism, base_dir, name):
    # Load in junctions
    junc_df1 = SP.build_junction_df(junc_beds[0], gff3, fa_dict, organism=organism)
    junc_df2 = SP.build_junction_df(junc_beds[1], gff3, fa_dict, organism=organism)
    
    junc_df = SP.combine_junctions(junc_df1, junc_df2)
    #print junc_df

    # Compare peaks and junctions
    peaks_w_junc = SP.compare_peak_junc_df(peak_df, junc_df, organism=organism)
    peaks_w_junc = SP.score_peaks(peaks_w_junc, gff3, fa_dict)
    
    # Reformat dataframe - add index, sort so that the annotated intron is first in each cluster
    peaks_w_junc.index = peaks_w_junc['genome coord']
    peaks_w_junc['type index'] = np.where(peaks_w_junc['junction type'] == 'Annotated', 0, 1)
    peaks_w_junc = peaks_w_junc.sort_values('type index')
    peaks_w_junc.groupby(peaks_w_junc.index).first()
    peaks_w_junc = peaks_w_junc.drop(['index', 'type index'], axis=1)
    peaks_w_junc['intron tuple'] = zip(peaks_w_junc['transcript'].tolist(),peaks_w_junc['annotated intron size'].tolist())
    
    print "\nPeaks with corresponding exon-exon junctions:"
    print len(peaks_w_junc)
    print str(len(set(peaks_w_junc[~peaks_w_junc['type'].str.contains('prime')]['genome coord'])))+" unpredicted"
    
    peaks_w_junc.to_csv(base_dir+name+'_peaks_w_junc.csv')
    peaks_w_junc.to_pickle(base_dir+name+'_peaks_w_junc.pickle')
    
    return peaks_w_junc

def peak_branch_analysis(peak_df, branch_bams, gff3, fa_dict, organism, base_dir, name):
    branches = SP.list_branch_points(branch_bams[0], gff3, fa_dict, organism=organism)
    branch_df = SP.create_branch_df(branches, gff3, fa_dict, organism=organism)
    if len(branch_bams) == 2:
        branches2 = SP.list_branch_points(branch_bams[1], gff3, fa_dict, organism=organism)
        branch_df2 = SP.create_branch_df(branches2, gff3, fa_dict, organism=organism)
        branch_df = branch_df.append(branch_df2)

        bed1 = branch_bams[0].split('_sorted.bam')[0]+'.bed'
        bed2 = branch_bams[1].split('_sorted.bam')[0]+'.bed'
        cat_args = "cat {0} {1} > {2}_all_branches.bed".format(bed1, bed2, name)
        call(cat_args, shell=True)

        os.remove(bed1)
        os.remove(bed2)

    # Compare peaks and branches
    peaks_w_branch = branch_df[branch_df['genome coord'].isin(peak_df['genome coord'])]
    peaks_w_branch = peaks_w_branch.merge(peak_df[['type','genome coord']], right_on='genome coord', left_on='genome coord', how='left')
    peaks_w_branch.index = peaks_w_branch['branch coord']

    print "\nPeaks with corresponding branches:"
    print len(peaks_w_branch)
    print str(len(set(peaks_w_branch['genome coord'])))+" unpredicted"

    peaks_w_branch.to_csv(base_dir+name+'_peaks_w_branch.csv')
    peaks_w_branch.to_pickle(base_dir+name+'_peaks_w_branch.pickle')
    
    return peaks_w_branch

def peaks_only(CP_out, quant_bams, base_dir, name, gff3, fa_dict, organism):
    peak_df = SP.peak_to_seq_pipeline(CP_out[2], CP_out[0], CP_out[1], gff3, fa_dict, name=name+'_CP_peaks')
    peak_df.to_csv(base_dir+name+'_all_peaks.csv')
    
    quant_df = SP.quant_from_peak_df(peak_df, gff3, fa_dict, organism=organism)
    quant_df = SP.quantitate_junction_df(quant_bams, quant_df, gff3, organism=organism)
    
    quant_df.to_csv(base_dir+name+'_quantitation.csv')
    
    scatter = SP.SP_pipeline_scatters(quant_df, base_dir, name)

def main():
    '''Usage: run SP_pipeline.py config_file untagged_sample_name organism
    config file : file that lists all branch, junction and peak files
    untagged_sample_name : prefix for untagged sample
    organism : pombe, crypto or cerevisiae'''
    peak_bams = []
    quant_bams = {}
    junc_beds = []
    branch_bams = []
    untagged = None
    
    # Read configuration file
    organism = sys.argv[2]
    with open(sys.argv[1], 'r') as config:
        for line in config:
            info = line.split(',')
            if info[1].strip().startswith('A'):
                if "untagged" in info[1].lower():
                    untagged = info[0]
                else:
                    peak_bams.append(info[0])
                    quant_bams[info[1].strip()] = info[0]
            elif info[1].strip().startswith('B') or info[1].strip().startswith('W'):
                quant_bams[info[1].strip()] = info[0]
            elif info[1].strip() == 'junctions':
                junc_beds.append(info[0])
            elif info[1].strip() == 'branches':
                branch_bams.append(info[0])
    
    if untagged is None: 
        print "No untagged sample specified"
        return None
    else:
        print "Untagged sample: "+untagged
        
    ## Set up output names and directory
    name = sys.argv[1].split('/')[-1].split('_config')[0]
    base_dir = sys.argv[1].split(name)[0]
    if base_dir == '': base_dir = './'
    print "Output file location and prefix: "+base_dir+name
    
    ## Next run changepoint peak picking tools
    bg_dict = SP.make_bg_files(peak_bams+[untagged], organism)
    print "\n5' bedgraph files created:"
    print bg_dict.values()
    
    CP_out = SP.changepoint_peaks(bg_dict[peak_bams[0]], bg_dict[peak_bams[1]], bg_dict[untagged])
    print "\nChangepoint analysis output:"
    print CP_out

    ## Determine whether branch or junction information has been provided
    if len(junc_beds) >= 2:
        print "\nJunction bed files"
        print junc_beds
        use_junctions = True
        
        if len(branch_bams) >= 1:
            print "\nBranch bam files"
            print branch_bams
            use_branches = True
        elif len(branch_bams) == 0:
            print "\nNo data for branches, continuing with only junctions"
            use_branches = False
    
    elif len(junc_beds) == 0:
        print "\nNo data for junctions, continuing with only peaks"
        use_junctions = False

    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)
    
    if use_junctions:
        #### Generate peak df
        if name+'_peaks_w_branch.csv' not in os.listdir(base_dir) or name+'_peaks_w_junc.csv' not in os.listdir(base_dir):
            if name+'_all_peaks.pickle' not in os.listdir(base_dir):
                peak_df = SP.peak_to_seq_pipeline(CP_out[2], CP_out[0], CP_out[1], gff3, fa_dict, name=name+'_CP_peaks')
                peak_df.to_pickle(base_dir+name+'_all_peaks.pickle')
                peak_df.to_csv(base_dir+name+'_all_peaks.csv')
            else:
                peak_df = pd.read_pickle(base_dir+name+'_all_peaks.pickle')

        #### Junction to peak comparison
        if name+'_peaks_w_junc.csv' not in os.listdir(base_dir):
            print "\nGenerating peaks vs. junctions dataframe..."
            peaks_w_junc = peak_junction_analysis(peak_df, junc_beds, gff3, fa_dict, organism, base_dir, name)

        else: 
            peaks_w_junc = pd.read_pickle(base_dir+name+'_peaks_w_junc.pickle')
            print "\nPeaks vs. junction dataframe already exists"


        #### Branch to peak comparison
        if use_branches is True:
            if name+'_peaks_w_branch.csv' not in os.listdir(base_dir):
                print "\nGenerating peaks vs. branches dataframe..."
                peaks_w_branch = peak_branch_analysis(peak_df, branch_bams, gff3, fa_dict, organism, base_dir, name)
            else: 
                peaks_w_branch = pd.read_pickle(base_dir+name+'_peaks_w_branch.pickle')
                print "\nPeaks vs. branches dataframe already exists"

        #### Clean up dataframe for quantitation
        if name+'_quantitation.csv' not in os.listdir(base_dir):
            quant_df, lariat_df = SP.make_quant_df(peaks_w_junc, peaks_w_branch, gff3, fa_dict, organism=organism)
            quant_df = SP.find_score_branches_ppy(quant_df, peaks_w_branch, fa_dict)
            print "\nCounting reads in transcripts and at peaks..."
            quant_df = SP.quantitate_junction_df(quant_bams, quant_df, gff3, organism=organism)

            quant_df.to_pickle(base_dir+name+'_quantitation.pickle')
            quant_df.to_csv(base_dir+name+'_quantitation.csv')
            lariat_df.to_pickle(base_dir+name+'_lariats.pickle')
            lariat_df.to_csv(base_dir+name+'_lariats.csv')

            scatter = SP.SP_pipeline_scatters(quant_df, base_dir, name)

        else:
            quant_df = pd.read_pickle(base_dir+name+'_quantitation.pickle')
            scatter = SP.SP_pipeline_scatters(quant_df, base_dir, name)
        
    else:
        peaks_only(CP_out, quant_bams, base_dir, name, gff3, fa_dict, organism)
    
    for file in [x for x in os.listdir(base_dir) if x.endswith('pickle')]:
        os.remove(file)
    
    print "\n****Finished****"

if __name__ == "__main__":
    main()
    