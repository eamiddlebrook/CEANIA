# CEANIA
Calculate Exact Average Nucleotide Identity from Alignments

### Read in an aligned multi-fasta and send all pairwise ANI values to STDout.
Usage:
```
python ceania.py ~/${path_to_multifasta_file} ${threads}
```

### Can digest an alignment file with an arbetrary number of alignments and length.
#### Works for everything from:
Single gene alignemnts for 10 organisms
#### To:
Whole genome or concatenated transcript alignments for hundreds of bacteria

### Known Issues:
##### python multithreading seems to hang if ${threads} is greater than the number of physical cores (MacOS)...
