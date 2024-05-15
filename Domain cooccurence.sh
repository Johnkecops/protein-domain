#  awk '{if($4<_[$1] || _[$1]=="")_[$1]=$4} END{for(i in _)print i" "_[i]}' anig-overlap-all-sorted.list

# awk '{{p=$1,b=$4,e=$5, s=$3, l=$0} {if (ei >_[e]) && (si<_[s]) && (bi<=[e]) && (p==[pi])}{e = ei, s = si, l = li} else {{p=$1,b=$4,e=$5, s=$3, l=$0} {if (ei >_[e]) {e =ei}}{if(si<s){s=si, l=li}}{print("#",l)}' 
# cont..
# else if l=empty then proceed to next protein


# awk '{{p=$1,b=$4,e=$5, s=$3, l=$0} {if (ei >_[e]) && (si<_[s]) && (bi<=[e]) && (p==[pi])}{e = ei, s = si, l = li} else {{p=$1,b=$4,e=$5, s=$3, l=$0} {if (ei >_[e]) {e =ei}}{if(si<s){s=si, l=li}}{print("#",l)}' 
# cont..
# else if l=empty then proceed to next protein
#----

# Part 1: Process the input file and update the hash table (_) with the values from the 4th column ($4)
awk '{if($4<_[$1] || _[$1]=="")_[$1]=$4} END{for(i in _)print i" "_[i]}' anig-overlap-all-sorted.list

# Part 2: Process the input file and update some variables based on certain conditions
# This part assumes that the input file has the following format:
# protein_id start_position end_position score
# and that the variables [pi], [ei], [si], and [bi] are defined and initialized elsewhere in the script.

# Define the input file
BEGIN {
    FS = "\t"  # Set the field separator to tab
    file = "anig-overlap-all-sorted.list"
}

# Process each line of the input file
{
    # Set the values of p, b, e, s, and l based on the current line
    p = $1
    b = $4
    e = $5
    s = $3
    l = $0

    # Check if the end position (e) is greater than the existing value for the same protein (p)
    if (e > _[p]) {
        # Update the end position (e) with the current value
        _[p] = e

        # Check if the start position (s) is less than the existing value for the same protein (p)
        if (s < _[p, "start"]) {
            # Update the start position (s) with the current value
            _[p, "start"] = s
        }

        # Check if the score (l) is not empty
        if (l != "") {
            # Update the score (l) with the current value
            _[p, "score"] = l
        }
    }
}

# Print the keys (protein_ids) and their corresponding values (start_positions, end_positions, and scores) in the hash table
END {
    for (i in _) {
        if (_[i, "start"] != "") {
            print i "\t" _[i, "start"] "\t" _[i, "end"] "\t" _[i, "score"]
        }
    }
}