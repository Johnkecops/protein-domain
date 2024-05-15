

#!/usr/bin/bash

# This script performs several operations on tab-separated data files.
# It assumes that the input files are tab-separated and that the script is run in the same directory as the input files.

# Step 1: Extract 'spes id' from 255.ass.tab and save it to spes id.tab
grep 'spes id\t' 255.ass.tab > spes_id.tab

# Step 2: Extract the first and sixth columns from spes_id.tab and save it to spes_id.list
awk '{print $1, $6}' spes_id.tab > spes_id.list

# Step 3: Run the connect.awk script on spes_id.list and save the output to spes_id.net
./connect.awk (jabarkan perintah scriptnya) > spes_id.net

# Step 4: Sort the second and first columns of spes_id.net, remove duplicates, and save the count of each unique pair to total_species_per_domain.list
awk '$1>$3{next}{$2,$1,$3}' spes_id.net > spes_id.txt

# Step 5: Move bcc.tab to mammals_ncbi.tab and grep the mammals.tab file to create a list of ncbi_id
mv bcc.tab mammals_ncbi.tab
grep -f mammals.tab ncbi_id > mammals_ncbi.tab

# Step 6: Extract the first and second columns from mammals_ncbi.tab and save it to sp_dom.list
awk '{print $1,$2,$6}' mammals_ncbi.tab > sp_dom.list

# Step 7: Sort the first and second columns of sp_dom.list, count the occurrences of each unique pair, and print the count and original pair
awk 'BEGIN {OFS=FS="\t"}
{ if($3 > $1) {printf("%s %s %s", $3, $2, $1)}
else {printf("%s %s %s", $1, $2, $3)}
for(i=4; i<=NF; i++) {printf(" %s", $i)}
printf("\n") }' dummy_test_duplicates_file.txt | awk ' { arr[$1 $2 $3]++; if(m[$1 $2 $3]=="")
{m[$1 $2 $3]=$0} } END { for(i in arr) {print arr[i], m[i]} }'

