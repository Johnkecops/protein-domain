#!/bin/bash

# This script reads three columns from a file and uses them to run a command on each line.
# The file is assumed to have three columns separated by whitespace.

# Declare three empty arrays to store the data from the file.
declare -a array_a
declare -a array_b
declare -a array_c

# Use awk to read the first, second, and third columns from the input file, respectively, and store them in variables vara, varb, and varc.
vara=$(awk '{ print $1}' input_file.tab)
varb=$(awk '{print $2}' input_file.tab)
varc=$(awk '{print $3}' input_file.tab)

# Initialize a counter variable i to 0.
i=0

# Use a for loop to iterate over each element in the array vara.
# For each element, assign it to the corresponding index in array_a, and increment the counter variable i.
for num in "${vara[@]}"; do
   array_a[$i]=$num
   i=$i+1
done

# Reset the counter variable i to 0.
i=0

# Use another for loop to iterate over each element in the array varb.
# For each element, assign it to the corresponding index in array_b, and increment the counter variable i.
for num in "${varb[@]}"; do
   array_b[$i]=$num
   i=$i+1
done

# Reset the counter variable i to 0.
i=0

# Use another for loop to iterate over each element in the array varc.
# For each element, assign it to the corresponding index in array_c, and increment the counter variable i.
for num in "${varc[@]}"; do
   array_c[$i]=$num
   i=$i+1
done

# Use another for loop to iterate over each element in array_a.
# For each element, run the fastacmd command with the specified parameters, using the corresponding elements from array_b and array_c.
for ((i=0; i<${#array_a[@]}; i++)) do
  fastacmd -d /scr/genomes/Metazoan-Animals/Homo_sapiens/UCSC/hg18/${array_a[$i]}.fa -p F -L ${array_b[$i]},${array_c[$i]} -s ${array_a[$i]} -S 1
done

# #!/bin/bash

# # This script reads three columns from a file and uses them to run a command on each line.
# # The file is assumed to have three columns separated by whitespace.

# # Declare three empty arrays to store the data from the file.
# declare -a array_a
# declare -a array_b
# declare -a array_c

# # Use awk to read the first, second, and third columns from the input file, respectively, and store them in variables vara, varb, and varc.
# vara=$(awk '{ print $1}' /u/praktikum/Download/tt/annofilter3-temp1.tab)
# varb=$(awk '{print $2}' /u/praktikum/Download/tt/annofilter3-temp1.tab)
# varc=$(awk '{print $3}' /u/praktikum/Download/tt/annofilter3-temp1.tab)

# # Initialize a counter variable i to 0.
# i=0

# # Use a for loop to iterate over each element in the array vara.
# # For each element, assign it to the corresponding index in array_a, and increment the counter variable i.
# for num in "${vara[@]}"; do
#    array_a[$i]=$num
#    i=$i+1
# done

# # Reset the counter variable i to 0.
# i=0

# # Use another for loop to iterate over each element in the array varb.
# # For each element, assign it to the corresponding index in array_b, and increment the counter variable i.
# for num in "${varb[@]}"; do
#    array_b[$i]=$num
#    i=$i+1
# done

# # Reset the counter variable i to 0.
# i=0

# # Use another for loop to iterate over each element in the array varc.
# # For each element, assign it to the corresponding index in array_c, and increment the counter variable i.
# for num in "${varc[@]}"; do
#    array_c[$i]=$num
#    i=$i+1
# done

# # Use another for loop to iterate over each element in array_a.
# # For each element, run the fastacmd command with the specified parameters, using the corresponding elements from array_b and array_c.
# for ((i=0; i<${#array_a[@]}; i++)) do
#   fastacmd -d /scr/genomes/Metazoan-Animals/Homo_sapiens/UCSC/hg18/${array_a[$i]}.fa -p F -L ${array_b[$i]},${array_c[$i]} -s ${array_a[$i]} -S 1
# done