from rpy2.robjects.packages import importr
base = importr('base')
utils = importr('utils')
cp = importr('changepoint')
from rpy2 import robjects
from matplotlib import pyplot as plt
from subprocess import call
import os

def changepoint_peaks(bg1, bg2, untagged):
    changepoint_test = robjects.r('''
        changepoint_analysis <- function(bedgraph, threshold) {
        data <- read.table(bedgraph, sep="\t", header=F, stringsAsFactors=F)
        data <- data[ data[,4] >= threshold, ]
        gc()
        counter <- 0
        chrs <- unique(data[,1])
        output <- paste0(strsplit(bedgraph, ".bedgraph")[[1]], "_thresh", threshold, "_cps.txt")
        for (chr in chrs) {
            subdata <- data[ data[,1]==chr,]
            if (nrow(subdata) > 1) {
                changepoints <- cpt.mean(subdata[,4], method="PELT", class=F)
            } else { changepoints <- 1
            }
            result <- subdata[changepoints, c(1,2,4)]
            counter <- nrow(result) + counter
        }
        counter
        }
        ''')
    
    changepoint_analysis = robjects.r('''
        changepoint_analysis <- function(bedgraph, threshold) {
        data <- read.table(bedgraph, sep="\t", header=F, stringsAsFactors=F)
        data <- data[ data[,4] >= threshold, ]
        gc()
        chrs <- unique(data[,1])
        output <- paste0(strsplit(bedgraph, ".bedgraph")[[1]], "_thresh", threshold, "_cps.txt")
        for (chr in chrs) {
            subdata <- data[ data[,1]==chr,]
            if (nrow(subdata) > 1) {
                changepoints <- cpt.mean(subdata[,4], method="PELT", class=F)
            } else { changepoints <- 1
            }
            result <- subdata[changepoints, c(1,2,4)]
            write.table(result, file=output, append=T, quote=F, sep="\t", row.names=F, col.names=F)
        }
        output
        }
        ''')
        
    cp_values = {}
    bedgraphs = [bg1, bg2]
    for bedgraph in bedgraphs:
        cp_values[bedgraph] = []
        for thresh in [1,2,3,4,5,10,20,30,40,50,60,70,80,90,100]:
            npoints = changepoint_test(bedgraph, thresh)
            cp_values[bedgraph].append((thresh, int(npoints[0])))
        
    
    # Plot change in number of peaks as a function of changepoint
    fig, ax = plt.subplots(1)
    colors = [('royalblue','skyblue'),('orangered','orange')]
    final_max = 0
    for m, bg in enumerate(cp_values.keys()):
        name = bg.split('/')[-1].split('.bedgraph')[0]
        cps = zip(*cp_values[bg])[0]
        npeaks = zip(*cp_values[bg])[1]
        ax.plot(cps, npeaks, color=colors[m][0], label=name+' #peaks')

        deriv = []
        for n in range(len(npeaks)-1):
            deriv.append(abs(npeaks[n+1]-npeaks[n]))

        ax.plot(cps[:-1], deriv, color=colors[m][1], label=name+' Delta')

        mx_slope = cps[deriv.index(max(deriv))+1]
        if mx_slope > final_max:
            final_thresh = mx_slope
        
        ax.legend()
        ax.set_xlabel("Changepoint Threshold")
        ax.set_ylabel("Number of peaks called")
        
    ax.set_ylim(ax.get_ylim())
    ax.plot([final_thresh, final_thresh], ax.get_ylim(), '--', color='0.7')
    plt.show()
    plt.clf()
    
    print "Threshold chosen: "+str(final_thresh)
    
    output_names = []
    for bedgraph in [bg1, bg2, untagged]:
        name = changepoint_analysis(bedgraph, final_thresh)
        output_names.append(str(name).strip().split('"')[1])
        
    return output_names

def make_bg_files(bam_files, organism):
    path = os.path.dirname(os.path.realpath(__file__))
    if 'pombe' in organism: bg_genome = path+'/Genomes/Sp_for_bg.genome'
    elif 'crypto' in organism: bg_genome = path+'/Genomes/crypto_for_bedgraph.genome'
    elif 'cerev' in organism: bg_genome = path+'/Genomes/S288C_for_bedgraph.genome'
    
    bg_dict = {}
    for bam in bam_files:
        bg_name = bam.split('_sorted.bam')[0]+'.bedgraph'
        call_list = ["genomeCoverageBed -5 -ibam",bam,"-g",bg_genome,"-bg >",bg_name]
        script = ' '.join(call_list)
        call(script, shell=True)
        
        bg_dict[bam] = bg_name
    return bg_dict