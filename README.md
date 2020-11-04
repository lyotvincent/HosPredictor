# HosPredictor

HosPredictor is an application for virus host prediction.  
The programs we wrote in this work:  
*Data preprocessing  
*Download data & makeblastdb  
*Machine learning trainer  
*A prediction program based on rdf model and sequences alignment  

## Data preprocessing

Data preprocessing programs are in *data_preprocessing* folder.  
Raw data downloaded from NCBI Virus are preprocessed through these program, and then be used to machine learning.  

*take_all_fna_togethor.py* move all .fna file into a folder, because these file are all saved in various sub-folder, which is hard to use.  
*sequences_origin.csv* contains host information.
*filtering_null_items.py*, *filtering_dont_existing_accession.py* invalid data in *sequences_origin.csv*.  
*label_only_all_fna_with_class.py* labels filtered data with "Mammalia|Aves" & "others".  
*count_GC.py* extract features (such as length, GC_content, k-mers) from virus .fna files.  

## Download data & makeblastdb

*make_blast_db.ipynb* is used to download NCBI Virus data (complete reference genomes).  

## Machine learning trainer

*train_random_forest.ipynb* contain the codes which train logistic regression, SVM, random forest, K-Neighbor models using the data generated by data preprocessing.

## A prediction program

*predict_host.py*  

HostPredictor Manual  
An example for simple usage:  
python3 predict_host.py -query XXXX.fasta -blast_db /path/virusdb  

-h -help --h --help     get this manual.  
-query  (needed) the input virus genome sequencing data.  
-blast_db       (needed) the virus database for BLAST program (make_blast_db.ipynb is a simple download & build database program).  
-model   the machine learning model for identifying if the virus data belongs to "Mammalia|Aves" or not. Default "rdf_clf_2020_10_17_17_33_14.pkl."  
-transformer     standardscaler for virus data, data will be processed by it before entered into classifier model. Default "transform_pipeline_standardscaler.pkl"  
-host_file       the file contained host information. Default "sequences_only_all_fna.csv"  
