__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt
import numpy as np
import random
from beeswarm import *
import operator
import re

#################################################################
## Build dataframes from tables written by crunchBAM           ##
#################################################################

def build_tables(file1, file2=None, header='infer', skiprows=None, low_memory=False, multiIndex=False, int_index=False, sep='\t'):
    #For junction reads, cleavages and feature counts use header=[0,1], skiprows=[2] and multiIndex=True. For all other tables, use default values.
    
    if file2==None:
        fin = open(file1, "r")
        df = pandas.read_table(fin, header=header, skiprows=skiprows, sep=sep)
        index_list = []
        if multiIndex is True:
            for row in df.iterrows():
                index_list.append((row[1][0],row[1][1]))
            df.index = pandas.MultiIndex.from_tuples(index_list)
        elif multiIndex is False and int_index is False:
            for row in df.iterrows():
                index_list.append(row[1][0])
            df.index = index_list

    else:
        fin = open(file1, "r")
        fin2 = open(file2, "r")
        df1 = pandas.read_table(fin, header=header, skiprows=skiprows)
        df2 = pandas.read_table(fin2, header=header, skiprows=skiprows)
        index_list1 = []
        index_list2 = []
        if header != 'infer':
            for row in df1.iterrows():
                index_list1.append((row[1][0],row[1][1]))
            for row in df2.iterrows():
                index_list2.append((row[1][0],row[1][1]))
            df1.index = pandas.MultiIndex.from_tuples(index_list1)
            df2.index = pandas.MultiIndex.from_tuples(index_list2)
        else:
            for row in df1.iterrows():
                index_list1.append(row[1][0])
            df1.index = index_list1
            for row in df2.iterrows():
                index_list2.append(row[1][0])
            df2.index = index_list2
        df1_sorted = df1.sort_index()
        df2_sorted = df2.sort_index()
        df = pandas.merge(df1_sorted, df2_sorted, left_index=True, right_index=True, copy=False, sort=True)
    
    
    columns_list = []
    n = 0
    if ("Transcript0","Transcript0") in df.columns:
        n = 1
    for column in df.columns:            
        if header != 'infer':
            if column[0].startswith("Unnamed: 0"):
                columns_list.append(("Transcript"+str(n),"Transcript"+str(n)))
                n += 1
            elif column[0].startswith("Unnamed: 1"):
                columns_list.append(("Exon"+str(n),"Exon"+str(n)))
                n += 1
            else:
                columns_list.append((column[0].split("_")[0], column[1]))
        else:
            if column.startswith("Unnamed: 0"):
                columns_list.append(("Transcript"+str(n),"Transcript"+str(n)))
                n += 1
            elif column == "Transcript":
                columns_list.append(("Transcript"+str(n),"Transcript"+str(n)))
                n += 1
            else:
                columns_list.append((column.split("_")[0],"Total"))
    df.columns = pandas.MultiIndex.from_tuples(columns_list)
    df = df.sort_index()
    return df

##################################################################
## Build dataframes from tables written by crunchBAM.           ##
## This function concatenates the Transcript and Intron columns ##
## with a "-" to compare to other lists, such as splice site    ##
## scores.                                                      ##
##################################################################

def build_intron_tables(file):
    fin = open(file, "r")
    transcripts = []
    info = []
    header = []
    for line in fin:
        line_list = (line.strip()).split("\t")
        if line.startswith("CNAG"):
            if line_list[0][-2:] == "T0":
                transcripts.append(line_list[0]+"-"+line_list[1])
                a = 2
                info_list = []
                while a < len(line_list):
                    info_list.append(line_list[a])
                    a+=1
                info.append(info_list)
            else:
                transcripts.append(line_list[0]+"T0-"+str(int(line_list[1])+1))
                a = 2
                info_list = []
                while a < len(line_list):
                    info_list.append(line_list[a])
                    a+=1
                info.append(info_list)
        elif line.startswith("Transcript"):
            b = 1
            while b < len(line_list):
                header.append(line_list[b].strip())
                b += 1
    ds_info = pandas.Series(info)
    transcript_dict = dict(zip(transcripts,ds_info))
    df = pandas.DataFrame()
    df = df.from_dict(transcript_dict, orient='index')
    df.columns = header[1:len(df.columns)+1]
    return df

####################################################################
## Normalize counts in A samples to counts in B samples. Also     ##
## divides by transcript length. c1-c4 are read counts for        ##
## internal control RNA. Arguments are output from build_tables   ##
####################################################################

def normalize_AtoB(feature_counts, total_counts, lengths, config, cutoff=0):
    df1 = feature_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    df1_indexed = df1.set_index(df1[("Transcript0","Transcript0")])
    df2_indexed = df2.set_index(df2[("Transcript0","Transcript0")])
    df3_indexed = df3.set_index(df3[("Transcript0","Transcript0")])
    #merged = pd.DataFrame(index=df1.index)
    #for name_tuple, row in 
    merged = pandas.merge(df1_indexed, df2_indexed, left_index=True, right_index=True, how = 'left')
    merged = pandas.merge(merged, df3_indexed, left_index=True, right_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    print merged.columns
    for name, read_type in merged.columns:
        names.append(name)
        if name.startswith("Transcript"):
            print "Reading table"
        elif n < len(configNames) and name == configNames[n][0]+"-A":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            
            #Do the math - (5' cleavages/control counts)/(B sample reads/(Length*control counts))
            merged[(name.split("-")[0],"5prime Normalized")] = pandas.Series((((merged[(name,"5prime")])/c1)/((merged[(name.strip("-A")+"-B","Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)
            merged[(name.split("-")[0],"3prime Normalized")] = pandas.Series((((merged[(name,"3prime")])/c1)/((merged[(name.strip("-A")+"-B","Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)

            n+=1
    #Implement cutoff if indicated
    if cutoff > 0:
        merged = merged[merged[(merged.columns[2][0].split("-")[0]+"-B","Total")] > cutoff]
        filtered_transcripts = merged.index.tolist()
        filtered_transcripts = list(set(filtered_transcripts))
        filtered_transcripts = sorted(filtered_transcripts)
        with open("top_transcripts.txt", "w") as fout:
            for transcript in filtered_transcripts:
                fout.write(transcript+"\n")
        print "Output top transcript list in top_transcripts.txt"
                
    
    merged = merged.fillna(0)
    print "New column names: 5prime Normalized, 3prime Normalized"
    print str(len(merged))+" rows"
    print str(len(merged.columns))+"columns"
    
    #Set up multiIndex for rows
    index_list = []
    for row in merged.iterrows():
        index_list.append((row[1][0],row[1][1]))
    merged.index = pandas.MultiIndex.from_tuples(index_list)
    return merged

#################################################################
## Normalize counts in A samples to counts in B samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA. Arguments are output from build_tables##
#################################################################

def normalize_JunctiontoTotal(feature_counts, total_counts, lengths, config, cutoff=0):
    df1 = feature_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    df1_indexed = df1.set_index(df1[("Transcript0","Transcript0")])
    df2_indexed = df2.set_index(df2[("Transcript0","Transcript0")])
    df3_indexed = df3.set_index(df3[("Transcript0","Transcript0")])
    merged = pandas.merge(df1_indexed, df2_indexed, left_index=True, right_index=True, how = 'left')
    merged = pandas.merge(merged, df3_indexed, left_index=True, right_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name, read_type in merged.columns:
        names.append(name)
        if name.startswith("Transcript"):
            print "Reading table"
        elif n < len(configNames) and name == configNames[n][0]+"-B":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            merged[(name.split("-")[0],"5prime junction")] = pandas.Series((((merged[(name,"5p junct cnts")])/c1)/((merged[(name,"Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)
            merged[(name.split("-")[0],"3prime junction")] = pandas.Series((((merged[(name,"3p junct cnts")])/c1)/((merged[(name,"Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)
            #merged[(name.split("-")[0],"5prime junction")] = pandas.Series((merged[(name, "5p junct cnts")])/c1), index = merged.index)
            #merged[(name.split("-")[0],"3prime junction")] = pandas.Series((merged[(name, "3p junct cnts")])/c1), index = merged.index)
            n+=1
    
    if cutoff > 0:
        merged = merged[merged[(merged.columns[2][0].split("-")[0]+"-B","Total")] > cutoff]
    merged = merged.fillna(0)
    print "5prime junction; 3prime junction"
    print len(merged)
    set_merged_index = set()
    for name1, name2 in merged.iterrows():
        set_merged_index.add(name1[0])
    print len(set_merged_index)
    return merged



#################################################################
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA. Arguments are output from build_tables##
#################################################################

def normalize_B_to_mature(total_counts, lengths, config, cutoff=0, cutoff_sample=None, untagged=None):
    df1 = total_counts
    df2 = lengths
    df3 = config
    df1_indexed = df1.set_index(df1[("Transcript0","Transcript0")])
    df2_indexed = df2.set_index(df2[("Transcript0","Transcript0")])
    df3 = df3.reset_index()
    merged = pandas.merge(df1_indexed,df2_indexed, left_index=True, right_index=True, how = 'left')
    #configNames = []
    configValues = {}
    for name in df3.columns:
        #configNames.append(name)
        s = pandas.Series(df3[name])
        l = s.tolist()
        #controlValues.append(l)
        configValues[name[0]] = l 
    names = []
    for name, read_type in merged.columns:
        names.append(name)
        if name.startswith("Transcript"):
            print "Reading table"
        elif name.endswith("B") and name.split("-")[0] in configValues:
            print "Reading sample " +name
            c1 = configValues[name.split("-")[0]][0]
            print "Control value 1 = " +str(c1)
            c2 = configValues[name.split("-")[0]][1]
            print "Control value 2 = " +str(c2)
            c3 = configValues[name.split("-")[0]][2]
            print "Control value 3 = " +str(c3)
            
            #For no untagged sample: (B totals/control reads)/(W totals/control reads)
            if untagged == None:
                merged[(name.split("-")[0],"Normalized to mature")] = pandas.Series(((merged[(name,"Total")]/c2)/(merged[(name.split("-")[0]+"-W","Total")]/c3)), index = merged.index)
            
            #For including untagged sample: ((B totals/control reads)/(W totals/control reads))/(B totals untagged/W totals untagged)
            elif untagged != None:
                merged[(name.split("-")[0],"Normalized to mature")] = pandas.Series(((merged[(name,"Total")]/c2)/(merged[(name.split("-")[0]+"-W","Total")]/c3))/(merged[(untagged+"-B","Total")]/merged[(untagged+"-W","Total")]/configValues[untagged][2]), index = merged.index)
               
    #Implement cutoff if indicated
    if cutoff > 0:
        merged = merged[merged[(cutoff_sample+"-B","Total")] > cutoff]
    
    merged = merged.fillna(0)
    print "New column name: Normalized to mature"
    print len(merged)
    return merged


####################################################################
## Filter transcripts based on an input table. In this case I'm   ##
## giving it a list of RNAi targets. Can be anything, 1 column    ##
## format with CNAG IDs. mergedtable is from one of the normalize ##
## functions. list_file is the name of a files containing the list##
## of CNAGs you want to filter by.                                ##
####################################################################

def read_list(list_file, no_T0=False):
    transcript_list = []
    with open(list_file, "r") as fin:
        for line in fin:
            if no_T0 is True:
                transcript_list.append(line.split("T0")[0])
            else: transcript_list.append(line.strip())
            
    return transcript_list

def filter_transcripts_by_cnag(mergedtable, list_file):
    geneID = []
    fin = open(list_file, "r")
    df = pandas.DataFrame(mergedtable)
    for line in fin:
        if line.startswith("CNAG"):
            gene = line.strip()
            geneID.append(gene)
            geneID.append(gene+"T0")
            geneID.append(gene+"T1")
            geneID.append(gene+"T2")
            geneID.append(gene+"T3")
            geneID.append(gene+"T4")
            geneID.append(gene+"T5")
    fin.close()
    RNAi_df = pandas.DataFrame(columns=df.columns)
    print len(df.index)
    if type(df.index[0]) == str:
        RNAi_df = df[(df.index.isin(geneID))]
    else:
        index_list = []
        for transcript in df.iterrows():
            if transcript[0][0] in geneID:
                RNAi_df = RNAi_df.append(RNAi_df.loc[transcript[0]])
                index_list.append(transcript[0])
                RNAi_df.index = pandas.MultiIndex.from_tuples(index_list)
    print len(RNAi_df)
    return RNAi_df

def filter_transcripts_by_name(mergedtable, list_file):
    geneID = []
    with open(list_file, "r") as fin:
        for line in fin:
            gene = line.strip()
            if len(gene) > 0 and gene[-2] != "T":
                n = 0
                while n < 6:
                    geneID.append(gene+"T"+str(n))
                    n += 1
    filtered_df = pandas.DataFrame(columns=mergedtable.columns)
    print len(mergedtable.index)
    if type(mergedtable.index[0]) == str:
        filtered_df = mergedtable[(mergedtable.index.isin(geneID))]
    else:
        index_list = []
        for transcript in mergedtable.iterrows():
            if transcript[0][0] in geneID:
                filtered_df = filtered_df.append(mergedtable.loc[transcript[0]])
                index_list.append(transcript[0])
                filtered_df.index = pandas.MultiIndex.from_tuples(index_list)
    print len(filtered_df)
    return filtered_df
            

def merge_tables(dflow, dfhigh, index_list_low=None, index_list_high=None):
    if index_list_low != None:
        dflow = dflow.filter(items=index_list_low, axis=1)
        print dflow.columns
    if index_list_high != None:
        dfhigh = dfhigh[index_list_high]
    merged = pandas.merge(dflow, dfhigh, left_index=True, right_index=True, how = 'left')
    print len(merged.columns)
    print merged.columns
    #if index_list != None:
    #    merged2 = merged.filter(axis=1, items=index_list)
    #    print len(merged)
    #else: merged2 = merged 
    print len(merged.columns)
    return merged

############################################################################
## Metanalysis of introns and exons based on position in gene from 5' end ##
############################################################################

def analyze_transcript_position(feature_counts, feature_type, cutoff=0):
    feature_counts_indexed = feature_counts.set_index(feature_counts[("Transcript0","Transcript0")])
    intron_count = 1
    read_dict = {}
    
    #Pick every row with intron number
    while intron_count < 10:
        read_list = []
        for row in feature_counts_indexed.iterrows():
            if  row[1][1] == intron_count:
                if feature_type == "Exon":
                    read_list.append(row[1][2])
                elif feature_type == "Intron":
                    read_list.append(row[1][3])
                else: print "Unknown feature type"
                
                #Append to dictionary - keys are intron number (or exon) and values are list of read counts
        print feature_type+" position "+str(intron_count)
        print "Before cutoff:"+str(len(read_list))+" "+feature_type+"s analyzed"
        read_list = [math.log(x, 10) for x in read_list if x > cutoff]
        read_dict[intron_count] = read_list
        print "After cutoff:"+str(len(read_list))+" "+feature_type+"s analyzed \n"
        intron_count += 1
            
    #Extract keys and lists for use in beeswarm
    values_list = []
    name_list = []
    for intron_number, read_list in read_dict.iteritems():
        name_list.append(intron_number)
        values_list.append(read_list)
        
    color_list = ["red", "orange", "yellow", "green", "teal", "blue", "purple", "indigo", "grey"]
        
    
    #Make beeswarm - x-axis is intron number and y-axis is counts
    bs, ax = beeswarm(values_list, method="swarm", labels=name_list, col=color_list)
    
    #Do I need to normalize to total RNA?

     
                
