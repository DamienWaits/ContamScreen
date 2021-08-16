#Script used to generate statistics on the results of a ContamScreen run.

#!/bin/bash
for file in *.fa
do
        name=`expr match "$file" '.*\(^[a-Z]*_[a-Z]*\)'`
        echo $name
        total=`grep -c ">" $file`
        good=`grep -c ">" $name*good.fasta`
        contam=`grep -c ">" $name*contam.fasta`
        suspect=`grep -c ">" $name*suspect.fasta`
        hits=`awk "BEGIN {print ( $good + $contam + $suspect ) * 100 / $total}"`
        echo -e "$hits \b% had hits"
        goodHits=`awk "BEGIN {print $good * 100 / $total}"`
        contamHits=`awk "BEGIN {print $contam * 100 / $total}"`
        suspectHits=`awk "BEGIN {print $suspect * 100 / $total}"`
        echo -e "$goodHits \b% good"
        echo -e "$contamHits \b% contam"
        echo -e "$suspectHits \b% suspect"
done
