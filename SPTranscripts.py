import re
import collections
import json

######################################################################################################
## Build a transcript dictionary with all information from a gff3 file                              ##
######################################################################################################

def build_transcript_dict(gff3_file, organism=None):
    '''Function to build dictionary of all transcripts in gff3 file.
    Parameters
    ----------
    gff3_file : str
                location and name of gff3 file
    organism : str, default ``None``
                only necessary for S. pombe - then use 'pombe'
                
    Returns
    --------
    transcript_dict : dict
                    Dictionary where transcript names are keys and values are list as follows:
                    [start, stop, strand, chromosome, [CDS starts], [CDS stops]]
                    '''
    
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3','IV':'chr4','V':'chr5','VI':'chr6','VII':'chr7','VIII':'chr8','IX':'chr9','X':'chr10','XI':'chr11','XII':'chr12','XIII':'chr13','XIV':'chr14','XV':'chr15','XVI':'chr16','MT':'MT'}
    
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = re.split(r'\t+', line)
            
            if organism == 'pombe' and len(columns) > 1:
                chr_rom = columns[0]
                chrom = rom_lat[chr_rom]
                transcript_types = ['transcript','pseudogene','rRNA','snoRNA','tRNA','snRNA']
                if columns[2] in transcript_types:
                    if columns[8].split(':')[0].split('=')[1] == 'gene': continue
                    transcript = columns[8].split(';')[0].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                elif columns[2] == 'CDS':
                    transcript = columns[8].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                        
            if len(columns) == 9 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8]
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"): transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    #Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]
                    
                    chrom = columns[0]
                    if chrom in rom_lat: chrom = rom_lat[chrom]
                    elif chrom not in rom_lat.keys():
                        if chrom[3:] in rom_lat: chrom = rom_lat[chrom[3:]]
                            
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                    else:
                        transcript_dict[transcript][0] = int(columns[3])
                        transcript_dict[transcript][1] = int(columns[4])
                        transcript_dict[transcript][2] = columns[6]
                        transcript_dict[transcript][3] = columns[0]
                elif columns[2] == "CDS":
                    transcript = columns[8].split("=")[1].split(".")[0].split(';')[0]
                    if 'mRNA' in transcript: transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    if transcript not in transcript_dict:
                        strand = columns[6]
                        chrom = columns[0]
                        transcript_dict[transcript] = [0,0,strand,chrom,[],[]]
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                    
    for tx in transcript_dict:
        if transcript_dict[tx][0] == 0:
            transcript_dict[tx][0] = transcript_dict[tx][4][0]
            transcript_dict[tx][1] = transcript_dict[tx][5][0]
    
    transcript_dict = collections.OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict

def tx_info(tx, tx_dict):
    strand = tx_dict[tx][2]
    start = tx_dict[tx][0]
    end = tx_dict[tx][1]
    chrom = tx_dict[tx][3]
    exons = None
    if len(tx_dict[tx]) > 4 and len(tx_dict[tx][4]) > 0:
        CDS_start = min(tx_dict[tx][4])
        CDS_end = max(tx_dict[tx][5])
        exons = zip(tx_dict[tx][4],tx_dict[tx][5])
    else:
        CDS_start = start
        CDS_end = end
        
    return start, end, strand, CDS_start, CDS_end, exons, chrom

##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome="All", gene_list=None, organism=None):
    '''Function to build dictionary of splice sites in each transcript.
    Parameters
    ----------
    gff3_file : str
                location and name of gff3 file
    chromosome : str, default "All"
               If not "All", will only list transcripts on the selected chromosome
    gene_list : list, default ``None``
               List of genes to limit the dictionary. Useful for optimizing performance of downstream functions.
    organism : str, default ``None``
                only necessary for S. pombe - then use 'pombe'
                
    Returns
    --------
    splice_site_dict : dict
                    Dictionary where transcript names are keys and values are list as follows:
                    [[intron starts], [intron stops]]
    intron flag : bool
                    True if introns rather than exons are defined in the gff3 file.
                    '''
        
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    
    splice_site_dict = {}
    n = 1
    with open(gff3_file,"r") as fin:
        for line in fin:
            columns = re.split(r'\t+', line.strip())

            if organism == 'pombe' and len(columns)>1:
                intron_flag = False
                if columns[2] == 'exon':
                    chr_rom = columns[0]
                    rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                    chrom = rom_lat[chr_rom]
                    transcript = columns[8].split(':')[0].split('=')[1]

                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],chrom]
                    if columns[6] == "+":
                        splice_site_dict[transcript][0].append(int(columns[4])-1)
                        splice_site_dict[transcript][1].append(int(columns[3])-2)
                    elif columns[6] == "-":
                        splice_site_dict[transcript][0].append(int(columns[3])-1)
                        splice_site_dict[transcript][1].append(int(columns[4]))

            if len(columns) > 1 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8].strip()
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"):
                        transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],columns[0]]

                elif columns[2] == "exon":
                    intron_flag=False
                    transcript = columns[8].strip()
                    transcript = transcript[-12:]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if gene_list is None:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[4])-1)
                            splice_site_dict[transcript][1].append(int(columns[3])-2)
                        if columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[3])-1)
                            splice_site_dict[transcript][1].append(int(columns[4]))
                    else:
                        if transcript in gene_list:
                            if columns[6] == "+":
                                splice_site_dict[transcript][0].append(int(columns[4])-1)
                                splice_site_dict[transcript][1].append(int(columns[3])-2)
                            if columns[6] == "-":
                                splice_site_dict[transcript][0].append(int(columns[3])-1)
                                splice_site_dict[transcript][1].append(int(columns[4]))

                #For organisms where introns are annotated instead of exons (e.g. S. cerevisiae)
                elif "intron" in columns[2]: 
                    intron_flag=True
                    transcript = columns[8].strip()
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"):
                        transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T':
                        transcript = transcript+'T0'
                    if transcript not in splice_site_dict:
                        splice_site_dict[transcript] = [[],[],columns[0]]
                    if gene_list is None:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[3])-1)
                            splice_site_dict[transcript][1].append(int(columns[4]))
                        elif columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[4]))
                            splice_site_dict[transcript][1].append(int(columns[3])-1)
                    else:
                        if transcript in gene_list:
                            if columns[6] == "+":
                                splice_site_dict[transcript][0].append(int(columns[3])-2)
                                splice_site_dict[transcript][1].append(int(columns[4])-1)
                            elif columns[6] == "-":
                                splice_site_dict[transcript][0].append(int(columns[4]))
                                splice_site_dict[transcript][1].append(int(columns[3])-1)
    
    
    #Trim to entries in gene list
    if gene_list is not None:
        splice_site_dict = {transcript: splice_site_dict[transcript] for transcript in gene_list}
        
    if chromosome != "All":
        splice_site_dict = dict([(transcript, coords) for transcript, coords in splice_site_dict.iteritems() if coords[2] == chromosome])
    
    if intron_flag is False:
        for transcript, sites in splice_site_dict.iteritems():
            sites5 = sorted(sites[0])
            sites3 = sorted(sites[1])
            if transcript_dict[transcript][2] == "+":
                sites5 = sites5[:-1]
                sites3 = sites3[1:]
            elif transcript_dict[transcript][2] == "-":
                sites5 = sites5[1:]
                sites3 = sites3[:-1]
            splice_site_dict[transcript] = [sites5, sites3]
        
    return (splice_site_dict, intron_flag)    

def collapse_ss_dict(splice_site_dict):
    ss_by_gene = {}
    transcripts = splice_site_dict.keys()
    transcripts = [x[:-2] for x in transcripts]
    genes = list(set(transcripts))
    for gene in genes:
        for transcript, sites in splice_site_dict.iteritems():
            if gene in transcript:
                if gene not in ss_by_gene:
                    ss_by_gene[gene] = set()
                ss_by_gene[gene].update(set(zip(sites[0],sites[1])))
    return ss_by_gene

def read_kallisto_abundance(abundance_file):
    kallisto_dict = {}
    #kallisto_dict[transcript] = [length, eff_length, est_counts, tpm]
    with open(abundance_file, 'r') as fin:
        for line in fin:
            if line.startswith('target'):
                continue
            columns = line.split('\t')
            kallisto_dict[columns[0]] = [int(columns[1]), float(columns[2]), float(columns[3]), float(columns[4].strip())]
    return kallisto_dict

def find_organism_files(organism):
    ''' usage: organism, gff3, fa_dict, bowtie_index = find_organism_files(organism)'''
    if 'crypto' in organism:
        gff3 = '/home/jordan/GENOMES/CNA3_FINAL_CALLGENES_1_gobs.gff3'
        fa = '/home/jordan/GENOMES/H99_fa.json'
        bowtie_index = '/home/jordan/GENOMES/Crypto_for_gobs'
        organism = None
    elif 'pombe' in organism:
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        fa = '/home/jordan/GENOMES/POMBE/Sp_fasta_dict.json'
        bowtie_index = '/home/jordan/GENOMES/POMBE/Spombe'
        organism = 'pombe'
    elif 'cerev' in organism:
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        fa = '/home/jordan/GENOMES/S288C/S288C_fasta_dict.json'
        bowtie_index = '/home/jordan/GENOMES/S288C/S288C'
        organism = None
        
    with open(fa) as f:
        fa_dict = json.load(f)
        
    return organism, gff3, fa_dict, bowtie_index