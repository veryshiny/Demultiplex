#!/usr/bin/env python

#command to run: 

#test file! WRITE THE STATEMENT
#./script_demultiplex.py -i TEST-input_FASTQ/indexes.txt -f1 TEST-input_FASTQ/test_data_R1.fq.gz -f2 TEST-input_FASTQ/test_data_R2.fq.gz -f3 TEST-input_FASTQ/test_data_R3.fq.gz -f4 TEST-input_FASTQ/test_data_R4.fq.gz 


import argparse
import bioinfo
import itertools
import gzip
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

def get_args():
    parser = argparse.ArgumentParser(description="to create a histogram of quality score of reads from fastq file")
    parser.add_argument("-f1", "--fileR1", help="input_fileR1", type= str)
    parser.add_argument("-f2", "--fileR2", help="input_fileR2", type= str)
    parser.add_argument("-f3", "--fileR3", help="input_fileR3", type= str)
    parser.add_argument("-f4", "--fileR4", help="input_fileR4", type= str)
    parser.add_argument("-i", "--index", help="indexfile", type= str)

    return parser.parse_args()

args = get_args()
# R1= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz","rt") #args.fileR1
# R2= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz","rt") #args.fileR2
# R3=gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz","rt") #args.fileR3
# R4= gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz","rt") #args.fileR4


R1=gzip.open(args.fileR1,"rt")
R2=gzip.open(args.fileR2,"rt")
R3=gzip.open(args.fileR3,"rt")
R4=gzip.open(args.fileR4,"rt")




# set_indexes={"GTAGCGTA"   ,	"CGATCGAT"    ,	"GATCAAGG",
# 	"AACAGCGA"    ,	"TAGCCATG"    ,	"CGGTAATC",
# 	"CTCTGGAT" ,   	"TACCGGAT"    ,	"CTAGCTCA",
# 	"CACTTCAC"   ,	"GCTACTCT"    ,	"ACGATCAG",
# 	"TATGGCAC"    ,	"TGTTCCGT"    ,	"GTCCTAAG",
# 	"TCGACAAG"    ,	"TCTTCGAC"    ,	"ATCATGCG",
# 	"ATCGTGGT"    , "TCGAGAGT"    ,	"TCGGATTC",
# "GATCTTGC"    , "AGAGTCCA"   ,	"AGGATAGC"}



#known indexes

indexes= args.index
first_line=True
known_indexes=set()
with open(indexes,"r") as fh:
    #looping through each line & splitting by \t
    for line in fh:
        if first_line==True:
            first_line=False
            continue
        line_parts = line.strip().split('\t')
        #append to list (column 5 (4 bc 0-based) has barcode seq.s)
        known_indexes.add(line_parts[4])


### permutation stuff ###

# Base_array="ATCG"


# permut= list(itertools.product(['A', 'C', 'T','G'], repeat=8))
# print((permut[1]))

# def convertTuple(tup):
#         # initialize an empty string
#     str = ''
#     for item in tup:
#         str = str + item
#     return str
# new_permut=[]

# for tuple in permut:
#     new_permut.append(convertTuple(tuple))

dict_index_pairs={}
for pairs in itertools.product(list(known_indexes), repeat=2): #itertools.product makes permutations
    #for eg if list is [a,b,c] and repeats =2 then output = (a,b) (b,c) (c,a) (a,a) (b,a) (c,c) as tuples
    dict_index_pairs[pairs]=0
    
dict_index_pairs["unknown"]=0

dict_matched_pairs={}
for pairs in list(known_indexes): 
    dict_matched_pairs[pairs]=0
    
dict_index_pairs["unknown"]=0

#reverse complement function

dict_DNA={"A":"T","G":"C","C":"G","T":"A","N":"N"}
def rev_complement(sequence: str) -> str:
    '''Takes a DNA sequence and for each base in the sequence, we complement it with the associated base i.e, A>T, G>C, C>G, T>A. We then reverse the
sequence ('hello world'[::-1]) to output finally the reverse complement'''
    revcomp=""
    sequence=sequence.upper()
    for letter in sequence:
        revcomp=dict_DNA[letter]+revcomp
    return revcomp


#threshold score function


# open all the files to write 

dict_R1={}
dict_R2={}

for index in known_indexes:
    rev=rev_complement(index)
    dict_R1[index]=open(f"output/{index}-{rev}_R1.fq","w")
    
    dict_R2[index]=open(f"output/{index}-{rev}_R2.fq","w")

hopped_R1=open(f"output/hopped_R1.fq","w")
unk_R1=open(f"output/unk_R1.fq","w")
hopped_R2=open(f"output/hopped_R2.fq","w")
unk_R2=open(f"output/unk_R2.fq","w")

record_with_4_lines_R1=""
record_with_4_lines_R4=""

counter=0 #to keep track of the lines to write into a file
i=0 

total_matched=0
total_hopped=0
dict_hopindex_counter={}
for index in known_indexes:
    dict_hopindex_counter[index]=0

for line_r1,line_r2,line_r3, line_r4 in zip(R1,R2,R3,R4):
    i+=1
    line_r1 = line_r1.strip("\n")
    line_r2 = line_r2.strip("\n")
    line_r3 = line_r3.strip("\n")
    line_r4 = line_r4.strip("\n")

    if i%4==2: # means its the sequence line!
        index1=line_r2  
        index2=line_r3
        
        
    if i%4==0: #means its the phred score line
        phred1=line_r2
        phred2=line_r3

    if counter<4: # means its not the last line of record
      
        record_with_4_lines_R1+=str(line_r1)+"\n" #we append the line
        record_with_4_lines_R4+=str(line_r4)+"\n" #we append the line
        # print(record_with_4_lines_R1)
        # print(record_with_4_lines_R4)

        counter+=1
        if counter==4: #after we add 1, if it is the last line then we copy the record into record_copy
            counter=0         
            reverse_index2=rev_complement(index2)
            new_header_r1=record_with_4_lines_R1.split("\n")[0]+" "+ str(index1) + "-" + reverse_index2
            new_header_r4=record_with_4_lines_R4.split("\n")[0]+" "+ str(index1) + "-" + reverse_index2
            record_with_4_lines_R1=new_header_r1+"\n"+record_with_4_lines_R1.split("\n")[1]+"\n"+record_with_4_lines_R1.split("\n")[2]+"\n"+record_with_4_lines_R1.split("\n")[3] +"\n"
            record_with_4_lines_R4=new_header_r4+"\n"+record_with_4_lines_R4.split("\n")[1]+"\n"+record_with_4_lines_R4.split("\n")[2]+"\n"+record_with_4_lines_R4.split("\n")[3] +"\n"

            
            if bioinfo.qual_score(phred1)<26 or bioinfo.qual_score(phred2)<26 or index1 not in known_indexes or reverse_index2 not in known_indexes:
                unk_R1.write(record_with_4_lines_R1)
                unk_R2.write(record_with_4_lines_R4)
                dict_index_pairs["unknown"]+=1
            elif index1 != reverse_index2:
                hopped_R1.write(record_with_4_lines_R1)
                hopped_R2.write(record_with_4_lines_R4)
                dict_index_pairs[(index1,reverse_index2)]+=1
                total_hopped+=1
                dict_hopindex_counter[index1]+=1
                dict_hopindex_counter[reverse_index2]+=1
            else:
                dict_R1[index1].write(record_with_4_lines_R1)
                dict_R2[index1].write(record_with_4_lines_R4)
                dict_index_pairs[(index1,reverse_index2)]+=1
                dict_matched_pairs[index1]+=1

                total_matched+=1
            record_with_4_lines_R1=""
            record_with_4_lines_R4=""  


#WRITING INTO OUTPUT FILES
table=open("output/index_pair_info.txt","w")
output=open("output/Summary_of_Output_run.txt","w")

table.writelines("index1\tRev Comp of index2\tindex2\tNumber of Appearances\tPercentage of Read Appearance\n")
output.writelines("                         Report of Demultiplexing Run\n\nThe overall table which tells us the number of each index pair appearances and their % appearance overall (rounded to 3 decimals) is in the file called index_pair_info.txt in the output folder \n")
output.writelines("There are 5 plots in the output folder: \n  1) Overall Distribution Piechart:\tOverall_distribution.png\n  2) Distribution of Matched Reads:\tmatched_distribution.png\n  3) Distribution of Individual Indexes which Hopped:\thopped_index_global_distribution.png\n  4) Top 20 Index-Pairs which hopped:\thopped20_index_pair_distribution.png\n  5) Heatmap of Number of Appearance of each Read:\theatmap_index_pair_distribution.png\n\n")

total_reads=total_hopped+total_matched+dict_index_pairs["unknown"]
output.writelines("                         Overall Statistics\n\n")
output.writelines(f"Total number of reads = {total_reads}\n")
output.writelines(f"Total number of samples = {len(known_indexes)}\n\n")
output.writelines(f"Overall matched samples % = {round(total_matched/total_reads*100,3)}\n")
output.writelines(f"Overall unknown samples % = {round(dict_index_pairs["unknown"]/total_reads*100,3)}\n")
output.writelines(f"Overall Index-Hopping % = {round(total_hopped/total_reads*100,3)}\n")

top_hop=max(dict_hopindex_counter, key=dict_hopindex_counter.get)
top_match=max(dict_matched_pairs, key=dict_matched_pairs.get)

output.writelines(f"\nIndex which was matched the most is {top_match} and it appears {round(dict_matched_pairs[top_match]/total_reads*100,3)} % overall and is {round(dict_matched_pairs[top_match]/total_matched*100,3)}% of the matched reads \n") 
output.writelines(f"Index which hopped the most is {top_hop} and it appears {round(dict_hopindex_counter[top_hop]/total_reads*100,3)} % overall and is {round(dict_hopindex_counter[top_hop]/total_hopped*100,3)}% of the hopped reads\n\n") 
output.writelines(f"                        Percentage of each matched sample\n\n")

dict_hop_pair={}#for the plot with hopped pairs

dict_index_pairs=dict(sorted(dict_index_pairs.items(), key=lambda item: item[1],reverse=True)) #for sorting the dict in descending order

array_for_heatmap=np.zeros((24,24)) #for heatmap
list_known_indexes=list(known_indexes) #for heatmap

for item in dict_index_pairs:
    
    if item =="unknown":
        table.writelines(f"unknown\tunknown\tunknown\t{dict_index_pairs[item]}\t{round(dict_index_pairs[item]/total_reads*100,3)}\n")
        # output.writelines(f"\nUnknown index pairs = {dict_index_pairs[item]}")
        # output.writelines(f"Unknown index pairs % = {dict_index_pairs[item]/total_reads}\n\n")
        continue
    else:
        table.writelines(f"{item[0]}\t{item[1]}\t{rev_complement(item[1])}\t{dict_index_pairs[item]}\t{round(dict_index_pairs[item]/total_reads*100,3)}\n")
        i=list_known_indexes.index(item[0])
        j=list_known_indexes.index(item[1])
        if i!=j:
            array_for_heatmap[i][j]=dict_index_pairs[item]
    if item[0]==item[1]:
        output.writelines(f"% of reads from Sample {item[0]} = {round(dict_index_pairs[item]/total_reads*100,3)}% \n")#
    else:
        dict_hop_pair[f"{item[0]}-{item[1]}"]=dict_index_pairs[item]


#######PLOTTING TIME!#######

#PLOT 1
dict_matched_pairs=dict(sorted(dict_matched_pairs.items(), key=lambda item: item[1]))
x= dict_matched_pairs.keys()
y= dict_matched_pairs.values()


fig = plt.figure()
plt.barh(x,y,color='crimson', edgecolor="orange")
# plt.xscale('log')
plt.xlabel("Number of Reads")
plt.ylabel("Index")
plt.title(f'Distribution of Matched Reads')
fig.set_size_inches(15, 10, forward=True)
plt.savefig("output/matched_distribution.png", dpi=fig.dpi)

#PLOT 2
dict_hopindex_counter=dict(sorted(dict_hopindex_counter.items(), key=lambda item: item[1]))
x= dict_hopindex_counter.keys()
y= dict_hopindex_counter.values()


fig = plt.figure()
plt.barh(x,y,color='crimson', edgecolor="orange")
# plt.xscale('log')
plt.xlabel("Number of Reads")
plt.ylabel("Index")
plt.title(f'Distribution of Indexes which Hopped')
fig.set_size_inches(15, 10, forward=True)
plt.savefig("output/hopped_index_global_distribution.png", dpi=fig.dpi)


#PLOT 3
labels = [f'Matched: {round(total_matched/total_reads*100,3)}%', f'Hopped: {round(total_hopped/total_reads*100,3)}%', f'Unknown: {round(dict_index_pairs["unknown"]/total_reads*100,3)}%']

sizes = [round(total_matched/total_reads*100,3), round(total_hopped/total_reads*100,3),round(dict_index_pairs["unknown"]/total_reads*100,3) ]
fig = plt.figure()

fig, ax = plt.subplots()
ax.pie(sizes, labels=labels)
fig.set_size_inches(10, 10, forward=True)
plt.title(f'Overall Distribution of Reads')
plt.savefig("output/Overall_distribution.png", dpi=fig.dpi)

#PLOT 4
dict_hop_pair=dict(sorted(dict_hop_pair.items(), key=lambda item: item[1],reverse=True)[:20])
dict_hop_pair=dict(sorted(dict_hop_pair.items(), key=lambda item: item[1]))

x= dict_hop_pair.keys()
y= dict_hop_pair.values()


fig = plt.figure()
plt.barh(x,y,color='crimson', edgecolor="orange")
# plt.xscale('log')
plt.xlabel("Number of Reads")
plt.ylabel("Index")
plt.title(f'Distribution of Top 20 Index-Pairs which Hopped')
fig.set_size_inches(15, 10, forward=True)
plt.savefig("output/hopped20_index_pair_distribution.png", dpi=fig.dpi)

#PLOT 5




fig, ax = plt.subplots()
im = ax.imshow(array_for_heatmap,cmap="RdPu",norm=colors.LogNorm())

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(list_known_indexes)), labels=list_known_indexes)
ax.set_yticks(np.arange(len(list_known_indexes)), labels=list_known_indexes)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
# textcolors=("black", "white")
# for i in range(len(list_known_indexes)):
#     for j in range(len(list_known_indexes)):
#         if array_for_heatmap[i, j] > 2:
#             text = ax.text(j, i, array_for_heatmap[i, j],
#                        ha="center", va="center", color=textcolors[1])
#         else:
#             text = ax.text(j, i, array_for_heatmap[i, j],
#                        ha="center", va="center", color=textcolors[0])

plt.colorbar(im,label="log(Number of hopped reads)") 
ax.set_title("Distribution of reads of each hopped Read Pair")
fig.tight_layout()
fig.set_size_inches(15, 10, forward=True)
plt.savefig("output/heatmap_index_pair_distribution.png", dpi=fig.dpi)



table.close()
output.close()


for index in known_indexes:
   
    dict_R1[index].close()    
    dict_R2[index].close()

hopped_R1.close()
unk_R1.close()
hopped_R2.close()
unk_R2.close()
R1.close()
R2.close()
R3.close()
R4.close()




