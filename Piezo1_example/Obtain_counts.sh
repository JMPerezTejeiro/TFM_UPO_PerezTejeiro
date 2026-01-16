#!/bin/bash

# merge all the files obtained with FeatureCounts
# this script msut be executed in the same directory where the featureCounts files for each sample of the dataset are located.

# First, we eliminate the output of this script, just in case we have executed it before, avoid a possible overlap
rm counts_m.txt

# Cut columns with Ensmbl id and counts
for num in *.txt;
do cut -f 1,7 $num > $num.cut;
done

# Rename columns
for num in *.cut;
do
  name=${num/\.txt\.cut/}
  # eliminate header
  tail -n +3 $num > $num.middle
  # add name of the sample
  echo -e "ID\t$name" | cat - $num.middle > $num.name

  rm $num.middle
done

# Eliminate .cut files
rm *.cut

# Join all files using awk
awk 'NF > 0 { a[$1] = a[$1] "\t" $2 } END { for (i in a) { print i a[i]; } }' *cut.name > counts_m.txt

# Put the header as first line:
# save the header as a variable
header=$(grep 'ID' counts_m.txt)
# eliminate the header
sed -i '/'"$header"'/d' counts_m.txt
# we add from the beggning
sed -i '1i '"$header"'' counts_m.txt

# eliminate cut.name files
rm *cut.name

echo "done :)"
