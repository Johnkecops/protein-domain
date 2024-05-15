# for p = $1, b=$4, e=$5, s=$3, l= $0

# do

# if (ei >_[e]) && (si<_[s]) && (bi<=[e]) && (p==[pi])}{e = ei, s = si, l = li}
#  else 
# if (ei >_[e]) && (si<_[s]) && (bi<=[e]) && (p==[pi])}{e = ei, s = si, l = li}
# print("#",l) 

# else 
# if {(ei <_[e]) && (bi>=[e]) && (p!==[pi])}{e = ei, s = si, l = li}
# else 
# if (ei >_[e]) && (si<_[s]) && (bi<=[e]) && (p==[pi])}{e = ei, s = si, l = li}
# print("#",l) 
# else if (p ==” “)
# end

# ‘p=p+pi ‘
# ‘b=b+bi’
# ‘e=e+ei’
# ‘s=s+si’
# ‘l=l+li’

# ------------------------------------------------------
# p=p1, b=b1; e=e1; s=s1; l=$0
# # p = protein , b = begin , e = end, s= e value exponent, l = line
# if (bi <= e && p == pi) {
# if (ei > e) {e = ei} 			# take the larger of e, ei
# if (si < s) { s = si; l = li} 		# take the smallest exponent, and keep the best hit
# }
# else { #new cluster
# print "#", l; 		# print best hit from previous cluster
# p=pi; b=bi. e=ei, s=si, l= $0;
# }

# awk '{if($11<_[$5" "$6] || _[$5" "$6]=="")_[$5" "$6]=$11} END{for(i in _)print i" "_[i]}'

#   awk '{if($4<_[$1] || _[$1]=="")_[$1]=$4} END{for(i in _)print i" "_[i]}' anig-overlap-all-sorted.list > min-nmbr
#  awk '{if($5>_[$1] || _[$1]=="")_[$1]=$5} END{for(i in _)print i" "_[i]}' anig-overlap-all-sorted.list > max-nmbr
# join min-nmbr max-nmbr > region-limit
#----
#!/bin/bash

# Set initial values for protein, begin, end, e-value exponent, and line content
p=0
b=0
e=0
s=0
l=""

# Iterate over each line in the input file
while IFS= read -r line; do

    # Extract values from the current line
    p_val=$(echo $line | cut -d' ' -f1)
    b_val=$(echo $line | cut -d' ' -f4)
    e_val=$(echo $line | cut -d' ' -f5)
    s_val=$(echo $line | cut -d' ' -f6)
    l_val=$(echo $line | cut -d' ' -f10)

    # Check if the current protein is the same as the previous one
    if [ "$p_val" -eq "$p" ]; then

        # If the current begin is less than or equal to the current end,
        # and the current protein is the same as the previous one,
        # update the end, e-value exponent, and line content if necessary
        if [ "$b_val" -le "$e_val" ] && [ "$p_val" -eq "$p" ]; then

            # If the current end is less than the current ei, update the end with the current ei
            if [ "$e_val" -lt "$ei" ]; then
                e=$ei
            fi

            # If the current e-value exponent is less than the current si, update the e-value exponent and line content
            if [ "$s_val" -lt "$si" ]; then
                s=$s_val
                l=$l_val
            fi

            # Print the updated line content if it's different from the previous one
            if [ "$l_val" -ne "$l" ]; then
                print "#$l_val"
                l=$l_val
            fi

        # If the current protein is different from the previous one,
        # print the best hit from the previous cluster and reset the variables
        else
            print "#$l"
            p=$p_val
            b=$b_val
            e=$e_val
            s=$s_val
            l=$l_val
        fi

    # If the current protein is different from the previous one,
    # print the best hit from the previous cluster and reset the variables
    else
        print "#$l"
        p=$p_val
        b=$b_val
        e=$e_val
        s=$s_val
        l=$l_val
    fi

done < input_file.txt

# Generate a summary file with minimum and maximum values for certain columns
awk '{if($11<_[$5" "$6] || _[$5" "$6]=="")_[$5" "$6]=$11} END{for(i in _)print i" "_[i]}' input_file.txt > min_nmbr
awk '{if($4<_[$1] || _[$1]=="")_[$1]=$4} END{for(i in _)print i" "_[i]}' input_file.txt > min_nmbr
awk '{if($5>_[$1] || _[$1]=="")_[$1]=$5} END{for(i in _)print i" "_[i]}' input_file.txt > max_nmbr
join min_nmbr max_nmbr > region_limit