Please upload your 4 *properly formatted* input test FASTQ files to this directory.

Make sure they cover all three categories (dual matched, index-hopped, unknown index).

I have zipped and unzipped versions to see the contents of the file clearly

Script: 
```bash
usr/bin/time  -v ./script_add_header.py -i /projects/bgmp/shared/2017_sequencing/indexes.txt -f1 test_data_R1.fq.gz -f2 test_data_R2.fq.gz -f3 test_data_R3.fq.gz -f4 test_data_R4.fq.gz 
```