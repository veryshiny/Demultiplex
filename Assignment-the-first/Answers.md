# Assignment the First

## Part 1
1. Link to python script: [histogram.py](histogram.py)


| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read 1  | 101 | phred +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1  | 8 | phred +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | phred +33  |
| 1294_S1_L008_R4_001.fastq.gz | read 2 | 101 | phred +33  |

2. Per-base NT distribution
    ![](https://github.com/veryshiny/Demultiplex/blob/master/Assignment-the-first/R1_hist.png)
    ![](https://github.com/veryshiny/Demultiplex/blob/master/Assignment-the-first/R2_hist.png)
    ![](https://github.com/veryshiny/Demultiplex/blob/master/Assignment-the-first/R3_hist.png)
    ![](https://github.com/veryshiny/Demultiplex/blob/master/Assignment-the-first/R4_hist.png)

3.  
    3.	What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis (**RNA-SEQ experiment**), respectively? Justify your answer.

    Illumina sequencing uses bins for their PHRED scores to minimize the data used. So each bin has 3 PHRED score values rounded down to the nearest multiple of 3. For eg, 28 would be 27, 4 would be 3, etc.

    For sample identification using indexes: we used a average quality score cut-off of 26 for all 8 positions on the index. The average of values for each position of the index is mostly around 30. We choose the bin right below this so that we provide some leeway for the scores, while still maintaining a good quality for most of the positions. We also did not want to go through and check every position of the index when the files were so large, as it would increase the processing time to demultiplex.

    For downstream analysis using : we can use an average quality score cut-off of 30, because the average of values for each position of the index is always above 30. We can afford to make this threshold more stringent or not, depending on the nature of our analyses.


    3.	How many indexes have undetermined (N) base calls? 
    for R2 (index 1): 3976613
    for R3 (index 2):3328051

    
## Part 2
1. Define the problem

We are given 4 files: R1, R2, R3 and R4. R1 and R4 are read files, while R2 is the index which the researchers used to match the reads from each sample and R3 is the reverse complement of the same. We are also given the list of 24 indexes which were used to generate the RNA-seq data. The goal is demultiplex the data, which means to get the R1 and R4 reads from each of the 24 samples into a pair of files each (48 in total: called <index>_R1.fq,<index>_R2.fq). We also want to check for the number of index-hopped samples, i.e, the ones which have correct indexes, but they do not match for a particular read between R2 and R3. These go in a pair of files (hopped_R1.fq,hopped_R2.fq) as well. Any samples which do not match the list of indexes we're given AND/OR do not cross our quality score threshhold need to be put into a pair of files as well (unk_R1.fq, unk_R2.fq).


   
2. Describe output
   The output consists of a total of 52 fastq files (26 pairs of R1 and R2), each header also includes the index sequence from R2 and the reverse complement of the index sequence from R3. The files are:
   - One pair for each pair of the 24 indexes (48)
   - One pair of index-hopped files (2)
   - One pair of unknown index files (2)
   
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
   ![](https://github.com/veryshiny/Demultiplex/blob/master/Assignment-the-first/pseudocode.png?raw=true)
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

```python
def rev_complement(sequence: str) -> str:
    '''Takes a DNA sequence and for each base in the sequence, we complement it with the associated base i.e, A>T, G>C, C>G, T>A. We then reverse the
sequence ('hello world'[::-1]) to output finally the reverse complement'''
    return revcomp
Input: ATTGC
Expected output: GCAAT
```

I also used bioinfo.py functions such as qual_score!
