# Task-3-Peptide-mass-analyzer
In this task the program reads in a standard Fasta format file and converts the sequence information into
masses. It considers them to all be +1 charged ions.
To calculate the mass of a peptide, it adds up the masses of its constituent amino
acids. The program can cope with monoisotopic and average masses based on user's selection. 

#Input:
File conanting peptide sequence, identifier and peptide number in Fasta format.

#Output:
Your output should be something of the form-
Prot_name   peptide   mass-to-charge   z   p   sequence
ABC1_YEAST   1         1234.5678       1   0 ALTAMCMNVWEITYHK
ABC1_YEAST   2         894.3019        1   0 SSLPMAR
ABC1_YEAST   3         1031.6791       1   0 MNASFFAAQL

Where “peptide” is the number of the peptide in the protein, z is the charge of the ion (+1 by default),
and p is the number of missed cleavages (relies on task 2 to pass this information on. If
its not available, just assume its 0).

All of the options should be provided on the command line
