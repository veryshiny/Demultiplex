# Lab Notebook :rose:

# Where the files are: /projects/bgmp/shared/2017_sequencing/

I cd-ed to there

```
$ ls -lah
total 46G
drwxrwsr-x+  3 coonrod  is.racs.pirg.bgmp 8.0K Apr 23 13:48 .
drwxrwsr-x+ 42 sdwagner is.racs.pirg.bgmp 8.0K Jan 22  2024 ..
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz
drwxrws---+  2 coonrod  is.racs.pirg.bgmp 8.0K Jul  1  2022 demultiplexed
-rwxrwxr-x+  1 sdwagner is.racs.pirg.bgmp  631 Aug  9  2021 indexes.txt
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  327 Aug 16  2017 README.txt
```

The files are BIG 😧 DON'T UNZIP 🤐

# my slurm input

```
srun -A bgmp -p bgmp --mem=100gb -c 8 --pty bash
```

# Initial data exploration

## Head-ing all the files
```
$ zcat 1294_S1_L008_R1_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ

$ zcat 1294_S1_L008_R2_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ

$ zcat 1294_S1_L008_R3_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF

$ zcat 1294_S1_L008_R4_001.fastq.gz | head -4
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--
```

## Checking the lengths of the sequences in each file

```
zcat 1294_S1_L008_R1_001.fastq.gz | sed -n '2~4p' | head -100000 | awk '{print l
ength}' | uniq  
101

zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | head -100000 | awk '{print l
ength}' | uniq  
8

zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | head -100000 | awk '{print l
ength}' | uniq  
8

zcat 1294_S1_L008_R4_001.fastq.gz | sed -n '2~4p' | head -100000 | awk '{print l
ength}' | uniq  
101

```

## looking for how a read looks when it has an actual index (B1	GTAGCGTA in this case)

from an index file! R2 which means that R2 is technically index 1 because the index sequence matches the list we were given
```
zcat 1294_S1_L008_R2_001.fastq.gz |grep -A 2 -B 1 "GTAGCGTA" | head -4
@K00337:83:HJKJNBBXX:8:1101:1438:1701 2:N:0:1
GTAGCGTA
+
AAFFFJJF
```

so I take the header of this sequence and search for it in the other files to check how its supposed to look??

**Remember to ctrl+C KILL IT after the first thing printed**
```
zcat 1294_S1_L008_R3_001.fastq.gz | grep -A 3 "@K00337:83:HJKJNBBXX:8:1101:1438:
1701" 
@K00337:83:HJKJNBBXX:8:1101:1438:1701 3:N:0:1
TACGCTAC
+
-AFAFJJJ
^C
```
The quality looking good (mayhaps?) THE INDEXES ARE REVERSE COMPLEMENTS YAY! no Ns in sight 👀

Now to check the reads!!

```
$zcat 1294_S1_L008_R1_001.fastq.gz | grep -A 3 "@K00337:83:HJKJNBBXX:8:1101:1438:1701"

@K00337:83:HJKJNBBXX:8:1101:1438:1701 1:N:0:1
GCCTTTCGTGCCTCCTGGGCTTCTAATTCCTTCTCATTCATGTCACGGTGTTTGACGACAAGCAAGGCATTGGTACCTGACTGAACCCCAGCCTTGGCCCG
+
AAA<F-FAJJAJJJJFFFAJFAJF<FAFJJJJJAJJJJJJJFFJJFJJJJJJJJJJJJJJJJJFJJJJJJJFFFFJJJJFJJJF-FFJFFJJJJJFJJJJJ

$zcat 1294_S1_L008_R4_001.fastq.gz | grep -A 3 "@K00337:83:HJKJNBBXX:8:1101:1438:1701"

@K00337:83:HJKJNBBXX:8:1101:1438:1701 4:N:0:1
TGTGGCTTATTTTCTGCCTGTGGAAGAAACACTAAAGAAACGAAAGCGGGACCAAGAGGAAGAGATGGACTATGCACCAGATGATGTGTATGACTATAAGA
+
AAA7AAFF<JJ<<FAFJJFJFAAAFJJFJJJJAFJJFJJJJFFJJFJFJFFJF<JJ<AFFJ<FAFJJJJJFJ7JJFFFJJJFJJJ7<F7JJJFF-AAFFAJ
^C
```

