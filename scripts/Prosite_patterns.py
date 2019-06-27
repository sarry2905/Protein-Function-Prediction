from Bio import SeqIO
from Bio.ExPASy import ScanProsite
import pandas as pd
import re
from Bio._py3k import urlopen
from time import time
import csv

#Input fasta file
with open('Sequence1.fasta', 'r+') as sequence_file:
    sequences = sequence_file.read()
count=0
sequences = sequences.split('>')[1:]
#Output CSV file name
with open("sequence_pattern.csv","a+",newline='') as pattern_file:
    for record in sequences[:]:
        newline = []
        rec = re.sub(r".*\d\n","",record)
        rec = rec.replace("\n","")
        #extract records for each sequence
        handle = ScanProsite.scan(seq=record, noprofile=1)      
        res = ScanProsite.read(handle)
        newline.append(rec)
        for i in range(0,len(res)):
            #separating signatures from the extracted records
            sig = res[i]['signature_ac']
            #make the url for the signature to extract the corresponding consensus patterns
            url ='https://prosite.expasy.org/{0}.txt'.format(sig)
            response = urlopen(url)
            content = response.read()
            content = content.decode("utf-8")
            patterns = re.findall("(\nPA)(.*)(\.\n)?",content)
            pattern =""
            for i in range(len(patterns)):
                pattern = pattern + patterns[i][1].strip()
                pattern = pattern.replace(".","")
            #save the pattern
            newline.append(pattern)
        count += 1
        writer = csv.writer(pattern_file)
        writer.writerow(newline)

