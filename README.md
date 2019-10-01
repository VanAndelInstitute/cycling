## Git Large File Support Is Required!

Please note that this repo contains large data files required to reproduce 
key elements of the gene expression analysis described in our paper:

<center>
<b>Isolation of cardiomyocytes undergoing mitosis with complete cytokinesis.</b>

Hsiao-yun Y. Milliron Ph.D1†, Matthew J. Weiland M.S.1†, Eric J. Kort M.D. M.S.1,2†, and Stefan Jovinge MD Ph.D1,3,4*

†Equal contribution

</center>

If you do not have Git LFS installed, you will only get the 
soft links to these files, not the files themselves, when you clone the 
repo. That is fine, unless you want to reproduce the analysis (as opposed 
to just reviewing the code).

Luckily, LFS support is easy to install. See here: https://git-lfs.github.com

## What is here

Supporting analysis scripts and data for the above manuscript enabling regeneration of the key panels from figure 5 and 6 related to single cell RNASeq expression profiling. The figure panes may be reproduced from the provided data files by executing the code blocks in  data_analysis.Rmd (or "knitting" the document from RStudio). 

This repository provides processed data necessary to perform the analysis and generate 
figures. Raw data (both fastqs and 10X genomics mtx files) have been submitted to GEO. We will 
update this page as soon as we have our GEO accession id.

Please see the "pre-requisites" section prior to running analysis to ensure 
required packages are installed. 
