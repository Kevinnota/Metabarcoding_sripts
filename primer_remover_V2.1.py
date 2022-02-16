######################################################
#                                                    #
#               remove primers                       #
#            correct orientation                     #
#                                                    #
######################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
import argparse
import gzip
from scipy.spatial.distance import hamming

def hamming_distance(a,b):
    counter=0
    for x,y in zip(a,b):
        counter=counter+(x!=y)  
    return(counter)
    
parser = argparse.ArgumentParser(description='Input and output')
parser.add_argument('-in', '--input', 
    help="Path to the input file", type=str)
parser.add_argument('-out', '--output', 
    help="Path to the output file", type=str) 
parser.add_argument('-l', '--lenght', 
    help="Minimum lenght for writing sequence after primer removal", type=int) 
    
args = parser.parse_args()

pl=len("GGGCAATCCTGAGCCAA")
pr=len("CCATTGAGTCTCTGCACCTATC")
new_fasta=[]
MEM=0
with gzip.open(args.input, "rt") as handle: #gzip
    for recods in SeqIO.parse(handle, "fastq"):
        if(len(recods.seq)>(pl+pr+8+8+10)):
            for i in range(len(recods.seq)-pr):
                S=0
                E=0
                hamdist = hamming_distance(str(recods.seq[(0+i):(pl+i)]), "GGGCAATCCTGAGCCAA")
                if(hamdist<=2):
                    break
                    
                hamdist = hamming_distance(str(recods.seq[(0+i):(pr+i)]), "CCATTGAGTCTCTGCACCTATC")
                if(hamdist<=2):
                    recods.seq=recods.seq.reverse_complement()
                    break   
                    
            for i in range(len(recods.seq)-pr):
                hamdist = hamming_distance(str(recods.seq[(0+i):(pl+i)]), "GGGCAATCCTGAGCCAA")
                if(hamdist<=2):
                    S=pl+i
                    break
            
            read_l=len(recods.seq)
            for i in range(len(recods.seq)-pr):        
                hamdist = hamming_distance(str(recods.seq[(read_l-i-pr):(read_l-i)]), "GATAGGTGCAGAGACTCAATGG")
                if(hamdist<=2):
                    E=(read_l-i-pr)
                    break
            
            if(len(str(recods.seq[(S):(E)]))>=args.lenght):
                new_seq=SeqRecord(Seq(str(recods.seq[(S):(E)])),
                    id=recods.id,
                    description=recods.description)
                new_fasta.append(new_seq)
            
            MEM+=1
            if(MEM==500000):
                print(MEM)
                MEM=0
                with open(args.output, "at") as handle_2:
                    SeqIO.write(new_fasta, handle_2, "fasta")
                new_fasta=[]

with open(args.output, "at") as handle_2:
                SeqIO.write(new_fasta, handle_2, "fasta")
                
print(r"""  
             _
            / \
           / D \
          /  O  \
         /   N   \
        /___ E ___\
            | |  """)
