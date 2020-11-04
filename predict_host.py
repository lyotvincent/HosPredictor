import os, time, joblib, sys
# import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

class HostPredictor:

    def __init__(self, query=None, model='./rdf_clf_2020_10_17_17_33_14.pkl', transformer='transform_pipeline_standardscaler.pkl', blast_db=None, host_file='./sequences_only_all_fna.csv'):
        self.query = query
        self.model = model
        self.transformer = transformer
        self.blast_db = blast_db
        self.host_file = host_file
    
    def run(self):
        start_time = time.time()
        if self.rdf_clf() >= 0.3:
            print("Input sequences maybe belong to \"Mammalia|Aves\"")
            self.seq_align()
        else:
            print("Input sequences maybe belong to \"others\"")
        time_consumption = time.time() - start_time
        print('Predition Use %ss' % time_consumption)

    def rdf_clf(self):
        transform_pipeline_standardscaler = joblib.load("transform_pipeline_standardscaler.pkl")
        rdf_clf = joblib.load("rdf_clf_2020_10_17_17_33_14.pkl")

        features = FeaturesExtractor(self.query).run()
        standard_features = transform_pipeline_standardscaler.transform(features)
        predicted_labels = rdf_clf.predict_proba(standard_features)
        return predicted_labels[0][0]

    def seq_align(self):

        with open(self.host_file, 'r') as f:
            host_lines = f.readlines()
        
        host_dict = dict()
        for line in host_lines:
            line = line.strip().split(',')
            host_dict[line[0]] = line[2]
        
        # accession.version %identity alignment_length
        try:
            blast_line = """blastn -query %s -db %s -outfmt 6 -num_threads 32|cut -f 2,3,4""" % (self.make_fraction_seq(), self.blast_db)
            res = os.popen(blast_line).read()
        except Exception as e:
            print(e)
        finally:
            os.remove('./temp_seq.fasta')
        hits = [list(x.split()) for x in res.splitlines()]
#         print(hits)

        predict_dict = dict()
        for i, hit in enumerate(hits):
            accession = hit[0].split('.')[0]
            if accession in host_dict.keys():
                if host_dict[accession] not in predict_dict.keys():
                    predict_dict[host_dict[accession]] = float(hit[1]) * float(hit[2]) / 100
                else:
                    predict_dict[host_dict[accession]] += float(hit[1]) * float(hit[2]) / 100
        
        predict_list = [[x, round(predict_dict[x], 2)] for x in predict_dict]
        predict_list = sorted(predict_list, key=lambda x:x[1], reverse=True)
        print(predict_list)
        return [x[0] for x in predict_list]

    def make_fraction_seq(self):
        in1 = open(self.query, 'r')
        lines = in1.readlines()
        in1.close()

        fna_seq = ''
        for fna_line in lines:
            if fna_line.startswith('>'):
                continue
            fna_seq += fna_line.strip()

        out1 = open('./temp_seq.fasta', 'w')
        while len(fna_seq) >= 1000:
            out1.write('>seq%s\n' % len(fna_seq))
            out1.write(fna_seq[:1000]+'\n')
            fna_seq = fna_seq[1000:]
        out1.close()

        return './temp_seq.fasta'


class FeaturesExtractor:
    
    def __init__(self, fasta=None):
        self.__fasta = fasta
        
    def run(self):
        
        bases = ['A', 'T', 'C', 'G']
        
        title = ['Length', 'GC_content']
        bases = ['A', 'T', 'C', 'G']
        for i in bases:
            for j in bases:
                for k in bases:
                    field = i+j+k+'_num'
                    title.append(field)
        
        result_list = list()
        
        fna = open(self.__fasta, 'r')
        fna_lines = fna.readlines()
        fna.close()
        fna_seq = ''
        for fna_line in fna_lines:
            if fna_line.startswith('>'):
                continue
            fna_seq += fna_line.strip()
        fna_seq.upper()
        gc_content = str( round( (fna_seq.count('G')+fna_seq.count('C'))/len(fna_seq), 2) )
        result_list.append(len(fna_seq))
        result_list.append(gc_content)
        for i in bases:
            for j in bases:
                for k in bases:
                    result_list.append(str(fna_seq.count(i+j+k)))
        return pd.DataFrame([result_list], columns=title)

def print_help():
    print("HostPredictor Manual")
    print("An example for simple usage:")
    print("python3 predict_host.py -query XXXX.fasta -blast_db /path/virusdb")
    print("---------------------------------------------------------------------")
    print("-h -help --h --help\tget this manual.")
    print("-query\t(needed) the input virus genome sequencing data.")
    print("-blast_db\t(needed) the virus database for BLAST program (make_blast_db.ipynb is a simple download & build database program).")
    print("-model\t the machine learning model for identifying if the virus data belongs to \"Mammalia|Aves\" or not. Default \"rdf_clf_2020_10_17_17_33_14.pkl.\"")
    print("-transformer\t standardscaler for virus data, data will be processed by it before entered into classifier model. Default \"transform_pipeline_standardscaler.pkl\"")
    print("-host_file\t the file contained host information. Default \"sequences_only_all_fna.csv\"")


if __name__ == "__main__":
    parameters = sys.argv
    print(str(parameters))

    if len(parameters) == 1 or "-h" in parameters or "-help" in parameters or "--h" in parameters or "--help" in parameters:
        print_help()
        exit()
    
    if "-query" in parameters:
        query = parameters[parameters.index('-query') + 1]
    else:
        print("-query is needed. Use -h to get manual.")
    
    if "-blast_db" in parameters:
        blast_db = parameters[parameters.index('-blast_db') + 1]
    else:
        print("-blast_db is needed. Use -h to get manual.")
    
    host_predictor = HostPredictor(query=query, blast_db=blast_db)

    if "-model" in parameters:
        model = parameters[parameters.index('-model') + 1]
        host_predictor.model = model
    
    if "-transformer" in parameters:
        transformer = parameters[parameters.index('-transformer') + 1]
        host_predictor.transformer = transformer
    
    if "-host_file" in parameters:
        host_file = parameters[parameters.index('-host_file') + 1]
        host_predictor.host_file = host_file
    
    host_predictor.run()
    
