# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. ![R1_hist](//projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/Assignment-the-first/R1_hist.png "R1 QS Histogram")
    ![R2_hist](//projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/Assignment-the-first/R2_hist.png "R2 QS Histogram")
    ![R3_hist](//projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/Assignment-the-first/R1_hist.png "R3 QS Histogram")
    ![R4_hist](//projects/bgmp/nkearns/bioinformatics/Bi622/Demultiplex/Demultiplex/Assignment-the-first/R1_hist.png "R4 QS Histogram")
    
    2. 
    ```
    The ultimate use of this RNAseq data is for differential gene expression. The next step in this experiment would be to compare the different sequencing libraries. A quality score cutoff of 20 would likely be appropriate for this purpose because the reads only need to be aligned back to a reference genome (so a highly quality genome has already been assembled). However, you would likely want a higher quality score cutoff (like 30) for the indexes to allow appropriate separation/demultiplexing of the data.```
    ```
    3. 
    ```
    From the index 1 file, there were 3976613 indexes that contained Ns
    I used the following command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | grep -F "N" | wc -l
    From the index 2 file, there were 3328051 indexes that contained Ns
    I used the following command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep -F "N" | wc -l
    ```
    
## Part 2
1. Define the problem
All the reads are in the same file, and we want to separate them into 24 groups by index. The indexes and reads are at the same positions in different files, but the indexes are not directly attached to each read. Another problem is that some of the indexes are swapped, and these incorrectly tagged reads need to be sorted out. Some indexes will also not match the original set, or will have "N" bases in them. Reads with these indexes will also need to be sorted out. Finally, not all reads will have sufficient quality scores. Low quality reads (with an average quality score of below 20) will have to be be sorted out. A relatively low QS (20, as opposed to 30 or 40) was selected for the reads because the reads are not being used to assemble a genome, the fragments only have to align to each other. A higher QS (30) was selected for the index because the indexes are essential for organizing/demultiplexing.

2. Describe output

48 files = one Read1 FASTQ file and one Read2 FASTQ file per matching index-pair (24 index pairs)
    - each read file needs to have the header adjusted to have the indexes on it (the read 1 and read 2 index)
2 files = two FASTQ files with index-hopped reads-pairs
    - If the indexes in the header are not the same, bin it into the index-hopped reads-pairs file
2 files = files undetermined (non-matching or low quality) index-pairs
    - The indexes don't match any of the 24 set indexes
    - Or if the indexes are low-quality (Ns in the index)

The algorithm should report:
    - the number of read-pairs with properly matched indexes (per index-pair) 
    - the number of read pairs with index-hopping observed 
    - the number of read pairs with unknown index(es). 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

       

4. Pseudocode

Open Read1, Index1, Index2, and Read2 - don't use with open
Make a name for each file that makes sense:
Index1 = 1294_S1_L008_R2_001.fastq.gz
Index2 = 1294_S1_L008_R3_001.fastq.gz 
Read1 = 1294_S1_L008_R1_001.fastq.gz
Read2 = 1294_S1_L008_R4_001.fastq.gz

Make an array of all the indexes
For each index in that array, create and open two files: (index#_R1.out, index#_R2.out)

Initialize a line counter
For line in file:
    Increment the counter
    Create a sliding window of three lines corresponding to the header (as done in PS5) so only the current record is held in memory as an array (record_array)
    Use readline() to read in a line and store it as a variable, then store those variables in the array
    Index1_record_array = [header, sequence, qscore_line]
    Index2_record_array = [header, sequence, qscore_line]
    Read1_record_array = [header, sequence, qscore_line]
    Read2_record_array = [header, sequence, qscore_line]

    Use Rev_Comp to create the reverse complement

    Use Convert_Phred to turn the phred score line into quality scores

    Use Average_QS to fine the average quality score from each line

    Create the total index from Index and Index 2. Total_Index = Index1-Index2_RC

    For Read1: 
        create new_header:
        new_header = current_header + total_index 

        Check if indexes contain N:
        If Index1 or Index2 contains an N > write the corresponding record from Read1 to R1_bad.out

        Check if index Qscore is below the cutoff(20), or if Read quality line is below the cutoff(30):
        Else if either index is below the cutoff (20) or read quality is below cutoff (30) > write the corresponding record from Read1 to R1_bad.out

        Use Index_Match to check if indexes match:
        If both indexes are in the known set of indexes: 
            If Index1 is the same as (==) Index2 > return TRUE and write the corresponding record from Read1 to index#_R1.out
            else if Index1 is different from (!=)  Index2 > return FALSE and write the corresponding record from Read1 to R1_swapped.out
        Else (one or both of the indexes is not in the known set of indexes):
            Write the corresponding record to R1_bad.out


    For Read2:
        create new_header:
        new_header = current_header + total_index 

        Check if indexes contain N:
        If Index1 or Index2 contains an N > write the corresponding record from Read2 to R2_bad.out

        Check if index Qscore is below the cutoff(20), or if Read quality line is below the cutoff(30):
        Else if either index is below the cutoff (20) or read quality is below cutoff (30) > write the corresponding record from Read2 to R2_bad.out

        Use Index_Match to check if indexes match:
        If both indexes are in the known set of indexes: 
            If Index1 is the same as (==) Index2 > return TRUE and write the corresponding record from Read2 to index#_R2.out
            else if Index1 is different from (!=)  Index2 > return FLASE and write the corresponding record from Read1 to R2_swapped.out
        Else (one or both of the indexes is not in the known set of indexes):
            Write the corresponding record to R2_bad.out

Close all the files!

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
```
Convert_phred(phred_score: str) -> int
```Takes as it's input each individual Phred score (character), and uses Phred33 encoding to convert it to an integer Quality Score.```
    Input: E
    Output: 36
    Return(Quality_Score)
```    
```
Average_QS(quality_line: str) -> float
```Takes as it's inpue an entire quality score line and finds the average quality score of that line```
    Input: GHF@EIJ (A quality line)
    Output: 37.42857142857143 (The average quality score of that line)
    Return(Average_QS)
```
```
Rev_Comp(index2: str) -> str
```Takes as it's input index2 and creates the reverse complement```
    Input: TACGCTAC
    Output: GTAGCGTA
    Return(index2_RC)
```
```
Index_Match(index1, index2_RC) -> bool
```Takes as it's input index 1 and index 2, compares them, and if they are the same, returns True. If the indexes don't match, it returns False.```
    Input: TACGCTAC, TACGCTAC
    Output: True
    Return(True or False)
```