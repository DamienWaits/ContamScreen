#This program screens a given assembly or assembly for contamination based on two 
#reference databases supplied by the user as fasta files, parameter 1 and parameter 2. 
#Parameter 1 should include sequences that are similar to those of the targeted assembly. 
#Parameter 2 should include sequences of known or suspected contaminants. These files 
#are merged and used to generate a BLAST database which is run on the target assembly. 
#The best hit for each contig in the assembly is chosen based on e-value score comparison. 
#E-values are compared based on order of magnitude and the sensitivity of this comparison 
#is decided by the user. Parameter 3 should be an integer and corresponds to how much 
#better a “contam” hit must be than the best “good” hit to be considered contaminated. 
#Example: If parameter 3 is input as 5 and the best “good” hit for Contig1 is 1e-10 and the 
#best “contam” hit is 1e-16, Contig1 is returned as “contaminated”. If the best “contam” 
#hit were 1e-13, Contig1 would be returned as “suspect” for manual evaluation. Parameter 
#4 is similar except it is how much better a “good” hit must be than a “contam” hit. Three 
#files are generated as output: TAXON-good.fasta, TAXON-contam.fasta, TAXON-
#suspect.fasta. BLAST output files can be found in the misc_intermediate_files directory.

#This program was developed by Damien S Waits and Kevin M Kocot. Version July 3rd, #2015

#!/bin/bash
#Creates a directory where intermediate output files are stored.
mkdir "misc_intermediate_files"

#Appends “good” to headers in the good file and writes “contam” #to headers in the contaminatant file.
sed -i 's/^>/>good\ /' $1 > good.fasta
sed -i 's/^>/>contam\ /' $2 > contaminated.fasta

#Writes both good.fasta and contaminated.fasta to an all.fasta and generates a blast database since e-values are not 
#informative across multiple blast searches.
cat good.fasta contaminated.fasta > all.fasta
makeblastdb -in all.fasta -dbtype nucl -title ALL -out ALL

#Performs the following actions on all files that have the suffix .fa until “done” is read.
for FILENAME in *.fa
do

#Removes line breaks from taxon.fa using nentferner.pl distributed with HaMStR
nentferner.pl -in=$FILENAME -out=$FILENAME".nent" 
rm $FILENAME
mv $FILENAME".nent" $FILENAME

#Removes unneeded header text that will make blast output unnecessarily long.
sed -i '/^>/ s/ .\+//g' $FILENAME
sed -i '/^>/ s/;.\+//g' $FILENAME

#Deletes sequences shorter than 100 amino acids in length. 
grep -B 1 "[^>].\{100,\}" $FILENAME > $FILENAME".tmp" 
#Deletes blank lines.
sed -i '/^$/d' $FILENAME".tmp"
#Deletes original file.
rm -rf $FILENAME
#Restores original file.
mv $FILENAME".tmp" $FILENAME

#Sets variable "taxon" equal to the part of the file name before the extension.
taxon=`echo $FILENAME | cut --delimiter=. --fields=1` 

#Performs blast search on assembly using the ALL database that includes both good and contam sequences,
#and formats output to a table. 
blastn -db ALL -query $FILENAME -num_descriptions 10 -num_alignments 10 -num_threads 6 > $taxon"_vs_all.txt" 
blast2table2.pl -format 10 -evalue 0.0001 $taxon"_vs_all.txt" > $taxon"_vs_all.table"

#Formats e-values to allow for comparison.
sed 's/[1-9]e-0*//' $taxon"_vs_all.table" > temp.table
sed -i 's/0.0/0/' temp.table

#Loops through all lines in the table format of the blast output.
while read line;
do

#Returns the sixth column of the table.
thisLine=`echo $line | awk '{print $6}'`

#The following nested if statements compares the e-values of hits to the same contig in the targeted assembly.
#Contaminant or true sequences are determined based on the e-value of the hit. If the hit is within a 
#threshold determined by the end user, it is instead output to a file of suspect sequences.
if [ "$lastLine" == "" ]
then
    lastLine=$thisLine
    currentValue=`echo $line | awk '{print $2}'`
    currentState=`echo $line | awk '{print $12}'`
    suspect=false
    if [ "$currentValue" -eq "0" ]
    then
        perfect=$currentState
        perfectLine=$thisLine
    fi
elif [ "$thisLine" == "$lastLine" ]
then
    newValue=`echo $line | awk '{print $2}'`
    newState=`echo $line | awk '{print $12}'`
    if [ "$newValue" -eq "0" ]
    then
        perfect=$newState
        perfectLine=$thisLine
    fi
    if [ "$currentState" == "contam" ] && [ "$newState" == "good" ]
    then
        if [ "$newValue" -gt "$((currentValue+$3))" ]
        then
            currentState=`echo $line | awk '{print $12}'`
            currentValue=$newValue
            suspect=false
        elif [ "$newValue" -gt "$((currentValue-$4))" ]
        then
            suspect=true
        fi
    elif [ "$currentState" == "good" ] && [ "$newState" == "contam" ]
    then    
        if [ "$newValue" -gt "$((currentValue+$4))" ]
        then
            currentState=`echo $line | awk '{print $12}'`
            currentValue=$newValue
            suspect=false
        elif [ "$newValue" -gt "$((currentValue-$3))" ]
        then    
            suspect=true
        fi
    else
        if [ "$newValue" -gt "$currentValue" ]
        then
            currentValue=$newValue
            currentState=`echo $line | awk '{print $12}'`
        fi
    fi
else
    if [ "$perfect" != "" ]
    then
        echo "$perfectLine" >> $perfect"_headers.txt"
    elif [ "$suspect" ]
    then
        echo "$lastLine" >> suspect_headers.txt
    elif [ "$currentState" == "good" ]
    then
        echo "$lastLine" >> good_headers.txt
    elif [ "$currentState" == "contam" ]    
    then    
        echo "$lastLine" >> contam_headers.txt
    fi
    lastLine=$thisLine
    currentValue=`echo $line | awk '{print $2}'`
    currentState=`echo $line | awk '{print $12}'`
    suspect=false
fi
    
done <  $taxon"_vs_all.table"

#Writes out the name of each contig to either good or contam headers files depending on the 
#comparisons from the above if statement.
if [ "$perfect" != "" ]
then        
        echo "$perfectLine" >> $perfect"_headers.txt"
elif  $suspect
then
        echo "$lastLine" >> suspect_headers.txt
elif [ $currentState == "good" ]
then
        echo "$lastLine" >> good_headers.txt
elif [ $currentState == "contam" ]
then
        echo "$lastLine" >> contam_headers.txt
fi

#Removes redundant sequence headers.
sort good_headers.txt | uniq > good_headers_sorted.txt
sort suspect_headers.txt | uniq > suspect_headers_sorted.txt
sort contam_headers.txt | uniq > contam_headers_sorted.txt

#Extracts the sequences that belong to the headers in the above files.
select_contigs.pl -n good_headers_sorted.txt $FILENAME
$taxon"-good.fasta"
select_contigs.pl -n contam_headers_sorted.txt $FILENAME $taxon"-contam.fasta"
select_contigs.pl -n suspect_headers_sorted.txt $FILENAME $taxon"-suspect.fasta"

#Clean-up of garbage files.
rm -f good_headers.txt
rm -f contam_headers.txt
rm -f suspect_headers.txt
rm -f temp.table
mv $taxon"_vs_all.txt" ./misc_intermediate_files/
mv $taxon"_vs_all.table" ./misc_intermediate_files/
done
