# Final-thesis - GAIA 
In this repository are located all the codes I used for my TFM. 

Here are present the codes I used in the pre-processing of the data using Bash. 

METAGENOMICS 

The aim of the project is to improve GAIA capacity to detect and classify strains. We’ll work on a machine learning model to try to develop a more accurate system. The first step is to find data, we used NCBI with SRA and ENA repositories. We used dataset with precise characteristics: 

-	Sequenced by Illumina platforms,
-	WGS,
-	MOCK communities with known strains

We found 7 datasets with these parameters and then we added the MOSAIC data that have been already proceed. The accession number of the selected datasets are: 

1. SRR6392275 
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6392275

2. SRR172902
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR172902

3. ERR3413876
https://www.ebi.ac.uk/ena/browser/view/ERR3413876

4. SRR2822457 
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR2822457

5. ERR3597814
https://www.ebi.ac.uk/ena/browser/view/ERR3597814

6. SRR172903
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR172903

7. SRR5275893
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5275893



Each dataset is contained is a folder called Data1, Data2 and so on. The data are in FASTQ format, so we performed quality control and trimming, using a home-made algorithm. 

1.	python /scripts/trimming.py/ -i Data4 -o Data4 -t 10 -v -P illumina

Then we ran GAIA using this command: 

2.	/scripts/gaia2/main_server.sh --pairs pairs.json --input /Trimmed_files --output Gaia_output --threads 10 --mode wgs --index /media/sequentia/synology_office/andreu/Gaia_db_202004/WGS+WTS/Prokaryotes/db.prokaryotes.1.fasta.gz

The pairs.json files have been created for each sequence with:

3.	echo ‘[{“forward”: SRR2822457_1.fastq.gz”, “reverse”:” SRR2822457_2. fastq”}}’ | jq “.” > pairs.json

Then we ran GAIA other three times. 
Why? Because we have million of prokaryotes’ sequences and it is not possible to archive them into one unique file. So, we have four databases and we have to use all of them with our sequences, to ensure higher coverage. So, we obtained from this step four folders: Gaia_output. Gaia_output2, Gaia_output3 and Gaia_output4, in each folder there are two important files: one is the parsed1.txt and the other is the strain.count.txt. 
Now we have to merge the results of parsed1.txt obtained from the 4 outputs. We first created a folder called merged_analysis and then we ran: 

4.	ls Gaia_output3/ | grep SRR | while read line ; do mkdir merged_analysis/$line ; done 

This command allows us to have a copy of the directory with the name of the dataset inside the merged analysis. Important: the folder is empty. Now we’ll run a function to obtain the merged parsed file. 

5.	ls Gaia_output3/ | grep SRR | while read line ; do cat Gaia_output1/$line/parsed1.txt Gaia_output2/$line/parsed1.txt Gaia_output3/$line/parsed1.txt Gaia_output4/$line/parsed1.txt | sort -k1,1V -S 5G | python /scripts/gaia2/merge.py > merged_analysis/$line/parsed1.txt ; done

Now we have to run GAIA for the last time, specifying as output merged_analysis folder and –start 3. Since this moment we’ll work on parsed file. What do it contain? There are 9 columns: 

ERR3413876.2	98	100	439375	0		439375_CP000758.1	1785239	1784948
ERR3413876.2	98	100	529	   0		     529_CP008820.1	     1246290	1245999
ERR3413876.3	98	100	1006551	0		1006551_CP003218.1	1943889	1943657
ERR3413876.3	98	100	1134687	0		1134687_CP008841.1	1632196	1631964
ERR3413876.3	98	100	1134687	0		1134687_CP022348.1	4097086	4096854
ERR3413876.3 98  100	1134687  0	    1134687_CP023185.1  2053546	2053871

Column1: name of the read 
Column2: identity %
Columns3: coverage %
Columns4: TAX ID 
Columns5: code indicating the multiple mapped, uniquely mapped and uniquely + unique region mapped (0,1,2).
Column6: empty 
Column7: TAX ID_locus (chromosome)
Columns8: bp in forward and reverse strand  

 We have to obtain now a BED file with strain ID and genome length.
 So, first we use: 

6.	grep ‘CP\|AE\|CY’ parsed1.txt > parsed1_unique.txt

We need only these patterns because we are looking for the whole genome length. In the databases db.prokaryotes.1/2/3/4.fasta.gz.fai – the format fai is a table –  there are a lot of sequences related to each strains, for examples related to a single protein and so on, so we have to focus our research only on the CP – AE – CY (from NCBI we know that they refer to the whole genome sequencing of those strains). 


7.	awk ‘($5 != 0){print;}’ parsed1_unique.txt > parsed1_final.txt


Other issue: we need those features that have a value in the column 5 different from zero. Why? This code indicates the multi-mapped (0), uniquely mapped (1 – 2 means that are uniquely mapped and that are mapping on a specific and unique region). Zero indicates multi-mapped reads, that are the sequences that map more than one time on the genome, so we have to ignore the features with this code. We are interest in those that have code = 1 or code = 2. Now, we have to take into account one thing: the features with code = 2 are strains’ candidates because if a strain is present inside the sample it has to be in a unique region and uniquely mapped.  So, why we decide to retrieve also the ones with code = 1 if then we’ll use only those with code = 2? We want to know if the coverage is uniform across the genome, it means if we have two strains A and B that belong to the same species, could be that they are different only for the presence of one gene in one of the two strains. So we can know which of the two strains is present in our sample considering the ones that has reads mapping on the unique region (the gene that the other strain hasn’t); and so we can conclude that the strain A is present and not the B in our sample. But, to confirm that the strain A is present in the sample it should have a uniform y constant coverage across the genome (if this not happen, it is a false positive). 

8.	awk '($5==2)' parsed1_final.txt | cut -f7 | sort -u | grep -w -F -f - db.prokaryotes.1.fasta.gz.fai > match1.txt
9.	awk '($5==2)' parsed1_final.txt | cut -f7 | sort -u | grep -w -F -f - db.prokaryotes.2.fasta.gz.fai > match2.txt
10.	awk '($5==2)' parsed1_final.txt | cut -f7 | sort -u | grep -w -F -f - db.prokaryotes.3.fasta.gz.fai > match3.txt
11.	awk '($5==2)' parsed1_final.txt | cut -f7 | sort -u | grep -w -F -f - db.prokaryotes.4.fasta.gz.fai > match4.txt

We used these four functions to obtain the genome length of the candidates’ strains comparing the parsed file with the databases in fai format. So we used awk to specify that we want to match only those features that have the code = 2 of the parsed file, then we cut the column 7 that contains the tax id + the locus (CP, AE or CY), because we only are interested in these tags. We used sort -u to eliminate duplicates and order the results and finally the function grep -w -F -f to intersect our list with the databases lists. 

12.	 for i in match_*.txt; do cat $i >> merged.txt; done 

Then we used this command to merge the results of the four steps below and then we use this: 

13.	awk ‘BEGIN{FS="\t"}{print $1"\t"$2"\t"}’ merged.txt > db.sizes

The output should be like this: 

223283_AE016853.1 6397126
314225_CP000157.1 3052398
292415_CP000116.1 2909809
290400_CP000264.1 4317977
159087_CP000089.1 4501104
323098_CP000115.1 3402093
259536_CP000082.1 2650701
269799_CP000148.1 3997420
203122_CP000282.1 5057531
316056_CP000301.1 5513844

We need this output to use bedtools makewindows with the chromosome and the genome length (it is not a BED file). Now we need to create a FASTA file with the sequences of the candidates’ strains. To do this we can use the Seqtk toolkit: 

14.	seqtk subseq /media/sequentia/synology_office/andreu/Gaia_db_202004/WGS+WTS/Prokaryotes/db.prokaryotes.1.fasta.gz strains.txt > out_1.fa

Another option could be using cat: 

15.	for i in $(cat strain.txt) ; do samtools faidx db.prokaryotes.1.fasta.gz $i; done  >> out_1.fa 

(we have to do the same for the 4 databases and then merge the outputs). 

The same function has been used for the other databases and the final output file is the result of the merged previous outputs. The strains.txt is the first line of the merged.txt file (containing the tax id and the chromosome of interest). We did the merge using this command: 
 
16.	awk ‘BEGIN{RS=">"; FS="\n"; ORS=""} (FNR==1){next} { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }  !(seen[seq]++){ print ">" $0 }' out1.fa out2.fa out3.fa > outfinal.fasta

To merge the output in FASTA format from the previous steps we have to specify with awk that we don’t want replicates, because this could affect our results. So: 
1.	Awk RS defines a line. Awk reads line by line by default. In our case we fix it into the symbol of the FASTA ID, 
2.	Awk FS is any single character or regular expression which you want to use as a input field separator, 
3.	Awk ORS is an Output equivalent of RS. Each record in the output will be printed with this delimiter.
4.	Awk FNR will give your number of records for each input file.
5.	The ‘g’ in gsub() stands for “global,” which means replace everywhere.
6.	The following command: !(seen[seq]++) remove the duplicates lines. 
7.	Then we specify our output

17.	bedtools makewindows -g db.sizes -w 1000 | nucBed -fi output_final.fa -bed stdin | awk 'BEGIN{FS="\t"}{print $1"\t"$2"\t"$3"\t"$5"\t"$10}' | grep -v "#" > db.windows
The output should be like this: 

            223283_AE016853.1	0	1000	0.571000	0
223283_AE016853.1	1000	2000	0.557000	0
223283_AE016853.1	2000	3000	0.565000	0
223283_AE016853.1	3000	4000	0.569000	0
223283_AE016853.1	4000	5000	0.536000	0
223283_AE016853.1	5000	6000	0.567000	0
223283_AE016853.1	6000	7000	0.553000	0
223283_AE016853.1	7000	8000	0.531000	0
223283_AE016853.1	8000	9000	0.541000	0
         223283_AE016853.1	     9000	10000	0.516000	0


The aim of this step is to create windows of 1000 bp of the genome (makewindows – bedtools), then we use nucBed to retrieve different information about those windows such as: 

-	%AT content
-	%GC content
-	Number of As observed
-	Number of Cs observed
-	Number of Gs observed
-	Number of Ts observed
-	Number of Ns observed
-	Number of other bases observed
-	The length of the explored sequence/interval.
-	The sequence extracted from the FASTA file. (optional, if -seq is used)
-	The number of times a user defined pattern was observed. (optional, if -pattern is used.)
Then we use awk to print only the columns of interest: chr, start and end position, GC % content and undefined nucleotides. Finally, with grep -v we exclude the # symbols that can lie in our data. 
Finally, to construct the final file:

18.	awk ‘BEGIN{FS="\t”}{if ($2!=”NA” && NF==9 && $2>=99 && $3>=99) {if ($8>$9) {start=$9 ; end=$8} else {start=$8 ; end=$9} ; print $7”\t”start”\t”end}}’ parsed1_final.txt | sortBed -i stdin | bgzip > trustful.alignments.bed.gz 

The file trustful.alignments.bed.gz is formed by the last three columns of parsed1_final.txt and it is a BED file, containing chromosome name, start and end position. 

19.	 bedtools intersect -a db.windows -b trustful.alignments.bed.gz -c | bgzip > trustful.alignments.windows.bed.gz

Then we do an intersect with bedtools of the previous file with the db.windows and we obtain the trustful.alignments.windows.bed.gz. This file would be like this: 

223283_AE016853.1	0	1000	0.571000	0	0
223283_AE016853.1	1000	2000	0.557000	0	0
223283_AE016853.1	2000	3000	0.565000	0	0
223283_AE016853.1	3000	4000	0.569000	0	0
223283_AE016853.1	4000	5000	0.536000	0	0
223283_AE016853.1	5000	6000	0.567000	0	0
223283_AE016853.1	6000	7000	0.553000	0	0

To create the final table: 

26.	nReads=$(cut -f1 parsed1_final.txt| uniq | wc -l) ; zcat trustful.alignments.windows.bed.gz | awk -v totalReads=$nReads ‘BEGIN{FS=”\t”}; {split($1,a,”_”) genomeLen[a[1]]+=$3-$2 ; mappedReads[a[1]]+=$6 ; if ($6>0) {binsCovered[a[1]]++} else {binsNotCovered[a[1]]++ ; GCcontentNotCovered[a[1]]+=$4 ; NcontentNotCovered[a[1]]+=$5/($3-$2)}}END{for (x in genomeLen) {print x”\t”genomeLen[x]”\t”binsCovered[x]/(binsCovered[x]+binsNotCovered[x])”\t”GccontentNotCovered[x]/+binsNotCovered[x]”\t”NcontentNotCovered[x]/binsNotCovered[x]”\t”mappedReads[x]”\t”totalReads}}’ | sort -k1,1V  > verdict_prova.txt 

Using this code we obtained a table in which there are 7 columns, the last the "Class" it is not present. Now in order to add the last columns indicating the presence or absence of the strain we used these codes: 

27.	 awk ‘NR==1 {print} NR>1 {printf(“%s\t%s\n”, $0, “0”) }’ verdict_prova.txt > verdict_2.txt

With the code above we added a last column in which there are all zero. Then we detected which tax ids are present in our dataframe, using grep -w (names of tax ids of all the strains of the paper). The ones that are present in our list have to be the 1 in the last column instead of the 0, so we use the following code to change it: 

28.	awk ‘BEGIN{FS=OFS=”\t”} $1~/^(tax ids of the strains that are present in the dataset)$/ {$8=”1”}1’ verdict_2.txt > verdict.txt

T
The solution of our analysis: 

29.	SRR6392275: 23/25 21 
30.	SRR172902: 7/22 7 
31.	ERR3413876: 21/30 8 
32.	SRR2822457: 2/20 2 
33.	ERR3597814: 0/10 0
34.	SRR172903: 2/20 2
35.	SRR5275893: 2/20  2

In the last dataset could be that becoming from an industry the tax ids are not present in NCBI? 

We have to take into account that in the paper there are the scientific names of the strains and we had to convert them into tax id, because our strains are annotated with the tax id. Is it possible that there are some errors in the conversion? There are some strains with a scientific name very long and complicated and in NCBI are simpler, so we don’t know if are exactly the same. However, in the first dataset we have obtained a significant number of detected strains 23/25. In the other cases, no. this have to be improved using machine learning method, training the algorithm to recognize the strains and correctly classify them. 

From this first cycle we have obtained the first verdict_final.txt that contains 1440 lines in total. However, this number is too low to perform machine learning approach, so we have decided to increase the final output using a down – sampling. What is it? In this case we use the parsed1_final.txt that contains 10M reads. We can choose randomly 5M and create another file parsed1_5M_final.txt and perform the analysis considering this file. Why we do this? Because we want to increase our samples but also because in this way, we can investigate how the final output changes considering a different number of reads, aka the coverage. Using this approach, we can answer to the following question: 
“Could we predict strains with low coverage, or we have to have always high coverage? What is the reads’ or bins’ total number that we need?” 
First of all, we created a copy of the whole folder SRR6392275_1, for example, called SRR6392275_1_2. 

# RESAMPLING 

36.	cp -r sampleA sampleA_2

Then we used the following command to do the down – sampling of the parsed1.txt: 

37.	cut -f1 parsed1.txt | sort -u | sort -R | head -n 5000000 | grep -w -F -f - parsed1.txt > parsed1_5M.txt

Since now we have to perform all the previous step, so select the “CP-AE-CY” and only the column 5 = 2 for the parsed1.txt and so on, until the creation of another verdict.txt. 


To eliminate heavy files that we are not interested in we can use: 

38.	find Gaia_output* -name "alignments.bam" | while read line ; do rm $line ; done
39.	find Gaia_output* -name "parsed2.txt" | while read line ; do rm $line ; done

           and to compress other files:

40.	find Gaia_output* -name "parsed1.txt" | while read line ; do gzip $line ; done







