######################################################
#                                                    #
#               unique                               #
#            sample_count                            #
#                                                    #
######################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
import argparse
import gzip
import re
import pandas as pd

parser = argparse.ArgumentParser(description='Input and output')
parser.add_argument('-in', '--input', 
    help="Path to the input file", type=str)

parser.add_argument('-out', '--output', 
    help="Path to the input file", type=str)
    
parser.add_argument('-c', '--mincount', 
    help="minimal total read count to be used", type=int)
    
args = parser.parse_args()
unique_seq=[]

input_sequences=list(SeqIO.parse(args.input, "fasta"))

data_both_df=pd.DataFrame()
data_start_df=pd.DataFrame()
data_end_df=pd.DataFrame()

first_new_read="true"
new_fasta=[]
print("number of reads = "+str(len(input_sequences)))

for i in range(len(input_sequences)):
    read=str(input_sequences[i].seq)
    
    sample_name_start=re.split(":", re.split("sample_name_start=",(input_sequences[i].description))[1])[0]
    sample_name_end=re.split(":", re.split("sample_name_end=",(input_sequences[i].description))[1])[0]  
    
    # detecting unique sequences and count sequences
    if(first_new_read=="true"):
        if(sample_name_start==sample_name_end):
            data_both = {'sample_both':[sample_name_start],
                    'count':[0]}
            first_new_read="false"
            data_both_df=data_both_df.append(pd.DataFrame(data_both))
            
        else:
            data_start = {'sample_start':[sample_name_start],
                        'count':[0]}
            data_start_df=data_start_df.append(pd.DataFrame(data_start))            
            data_end = {'sample_end':[sample_name_end],
                        'count':[0]}
            data_end_df=data_end_df.append(pd.DataFrame(data_end))        
            first_new_read="false"
    last_read="false"
    if ((i+1)==len(input_sequences)):
        i-=1
        print(i)
        last_read="true"
        
    if read==str(input_sequences[i+1].seq):
        if(sample_name_start==sample_name_end):
            added="false"
            for s in range(len(data_both_df)):
                if(sample_name_start==data_both_df.loc[s,'sample_both']):
                    data_both_df.loc[s,'count'] = data_both_df.loc[s,'count']+1
                    added="true"
            if(added=="false"):
                data_both = {'sample_both':[sample_name_start],'count':[1]}
                data_both_df=data_both_df.append(pd.DataFrame(data_both),ignore_index=True)           
        else:
            added="false"
            for s in range(len(data_start_df)):
                if(sample_name_start==data_start_df.loc[s,'sample_start']):
                    data_start_df.loc[s,'count'] = data_start_df.loc[s,'count']+1
                    added="true"
            if(added=="false"):
                data_start = {'sample_start':[sample_name_start],'count':[1]}
                data_start_df=data_start_df.append(pd.DataFrame(data_start),ignore_index=True)
                
            added="false"
            for s in range(len(data_end_df)):
                if(sample_name_end==data_end_df.loc[s,'sample_end']):
                    data_end_df.loc[s,'count'] = data_end_df.loc[s,'count']+1
                    added="true"
            if(added=="false"):
                data_end = {'sample_end':[sample_name_end],'count':[1]}
                data_end_df=data_end_df.append(pd.DataFrame(data_end),ignore_index=True)   
    
    # creating new headers      
    seper=", "
    if read!=str(input_sequences[i+1].seq) or (last_read=="true"):
        count_total=0
        if(len(data_both_df)!=0):
            count_total+=(data_both_df.loc[:,"count"].sum())
                
        if(len(data_start_df)!=0):
            count_total+=(data_start_df.loc[:,"count"].sum())
            
        new_header=[]
        for s in range(len(data_both_df)):
            if(data_both_df.loc[s,'count']!=0):
                new_header.append(str(data_both_df.loc[s,"sample_both"])+":"+str(data_both_df.loc[s,"count"]))
            
        header_join=seper.join(new_header)
        header_join="both_count={"+header_join+"}"
            
        new_header=[]
        for s in range(len(data_start_df)):
            if(data_start_df.loc[s,'count']!=0): 
                new_header.append(str(data_start_df.loc[s,"sample_start"])+":"+str(data_start_df.loc[s,"count"]))      
            
        header_join_start=seper.join(new_header)
        header_join_start="start_count={"+header_join_start+"}"
            
        new_header=[]
        for s in range(len(data_end_df)):
            if(data_end_df.loc[s,'count']!=0):       
                new_header.append(str(data_end_df.loc[s,"sample_end"])+":"+str(data_end_df.loc[s,"count"]))  
            
        header_join_end=seper.join(new_header)
        header_join_end="end_count={"+header_join_end+"}"
            
        new_header=[]
        HEADER="total_count="+str(count_total)+": "+header_join+": "+header_join_start+": "+header_join_end
           
        if (count_total>args.mincount):
            new_seq=SeqRecord(Seq(str(input_sequences[i].seq)),
                id=input_sequences[i].id,
                description=str(HEADER))
                    
            new_fasta.append(new_seq)
                    
            data_both_df=pd.DataFrame()
            data_start_df=pd.DataFrame()
            data_end_df=pd.DataFrame()
            first_new_read="true"
            
#SeqIO.write(new_fasta, "Unique_sample_count.fasta", "fasta")
handle = open(args.output, "w")
for record in new_fasta:
     handle.write(">%s %s\n%s\n" % (record.id, record.description, record.seq))
handle.close()

print(r"""
             _
            / \
           / D \
          /  O  \
         /   N   \
        /___ E ___\
            | |  """)
