#a perl script to remove the redundacy in the genscan prediction FASTA files. The command is:
#BEGIN{$/=">";$"=";"}/(.*?)\n(.+?)\s*>?$/s && push
#@{$h{$2}},$1;END{for(keys%h){print ">@{$h{$_}}\n$_\n"}}

#<redundant-genscan-prediction-files> > <non-redundant-genscanpred-files>
#This script will create a hash with all the sequences as we read them, and check for duplicates by seeing whether that hash element already exists.
#This will remove redundant entries and concatenate their description lines.

#!/usr/bin/perl
#use strict;
#use warnings;

# Check if the input file is provided as an argument
#if (@ARGV != 1) {
#    die "Usage: $0 <input_file>\n";
#}

#my $input_file = $ARGV[0];

# Open the input file for reading
#open(my $in_fh, '<', $input_file) or die "Could not open '$input_file' for reading: $!";

# Initialize a hash to store unique sequences
#my %unique_sequences;

# Initialize a hash to store the count of each sequence
#my %sequence_count;

# Initialize a hash to store the concatenated description lines for each sequence
#my %concatenated_descriptions;

# Read the input file line by line
#while (my $line = <$in_fh>) {
    # Check if the line starts with ">"
 #   if ($line =~ /^>/) {
        # Extract the sequence name (header)
  #      my $header = $line;

        # Remove the ">" character from the header
#        $header =~ s/^>//;

        # Check if the sequence has already been seen
#        if (exists $unique_sequences{$header}) {
            # Increment the count of the sequence
#            $sequence_count{$header}++;

            # Concatenate the description lines for the sequence
#            push @{$concatenated_descriptions{$header}}, $line;
 #       } else {
  #          # Add the sequence to the hash of unique sequences
   #         $unique_sequences{$header} = $line;

            # Initialize the count and description lines for the sequence
#            $sequence_count{$header} = 1;
 #           $concatenated_descriptions{$header} = [$line];
   #     }
  #  }
#}

# Close the input file
#close($in_fh);

# Print the non-redundant sequences and their concatenated description lines
#for my $header (keys %unique_sequences) {
 #   if ($sequence_count{$header} > 1) {
  #      print ">$header\n";
   #     print join("\n", @{$concatenated_descriptions{$header}}), "\n";
    #} else {
     #   print ">$header\n";
      #  print $unique_sequences{$header}, "\n";
    #}#
#}

#!/usr/bin/perl
use strict;
use warnings;

# Check if the input file is provided as an argument
if (@ARGV != 1) {
    die "Usage: $0 <input_file>\n";
}

my $input_file = $ARGV[0];

# Open the input file for reading
open(my $in_fh, '<', $input_file) or die "Could not open '$input_file' for reading: $!";

# Initialize a hash to store unique sequences
my %unique_sequences;

# Initialize a hash to store the count of each sequence
my %sequence_count;

# Initialize a hash to store the concatenated description lines for each sequence
my %concatenated_descriptions;

# Read the input file line by line
while (my $line = <$in_fh>) {
    # Check if the line starts with ">"
    if ($line =~ /^>/) {
        # Extract the sequence name (header)
        my $header = $line;

        # Remove the ">" character from the header
        $header =~ s/^>//;

        # Check if the sequence has already been seen
        if (exists $unique_sequences{$header}) {
            # Increment the count of the sequence
            $sequence_count{$header}++;

            # Concatenate the description lines for the sequence
            push @{$concatenated_descriptions{$header}}, $line;
        } else {
            # Add the sequence to the hash of unique sequences
            $unique_sequences{$header} = $line;

            # Initialize the count and description lines for the sequence
            $sequence_count{$header} = 1;
            $concatenated_descriptions{$header} = [$line];
        }
    }
}

# Close the input file
close($in_fh);

# Print the non-redundant sequences and their concatenated description lines
for my $header (keys %unique_sequences) {
    if ($sequence_count{$header} > 1) {
        print ">$header\n";
        print join("\n", @{$concatenated_descriptions{$header}}), "\n";
    } else {
        print ">$header\n";
        print $unique_sequences{$header}, "\n";
    }
}