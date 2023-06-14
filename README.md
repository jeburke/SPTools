# SPTools
## Tools for analyzing Spliceosome Profiling data  

# Change testing

Basic usage (Note: requires display environment, e.g. XQuarts on Mac or %matplotlib inline in Jupyter notebooks):  
```python SP_pipeline.py configuration_file untagged_sample organism```  
    configuration_file : see Example/Example_config.txt for format  
    untagged_sample : string (e.g. WT)  
    organism : string - can be 'crypto', 'pombe' or 'cerevisiae'  
    
Output files:  
>*_CP_peaks.bedgraph - Bedgraph formatted file with only reproducibly discovered peaks  
>*_all_peaks.csv - Spreadsheet of all peaks with surrounding sequence and peak height  
>*_peaks_w_junc.csv - Spreadsheet of all junctions with at least on corresponding peak and their properties  
>*_peaks_w_branch.csv - Spreadsheet of all branches that coincide with detected peaks  
>*_quantiation.csv - Spreadsheet with Precursor and Intermediate levels for each high confidence splicing event and all values used for building models  
>*_scatterplots.pdf - Scatter plots of replicates of all measurements  
