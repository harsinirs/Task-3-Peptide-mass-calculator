#
# Task 3: A mass analyser
#converter from sequence to mass (mass/charge), which accepts peptides/proteins in fasta format
#Harsini Raaja Sulochana 11179122
#


import sys
from tabulate import tabulate
def fastaread(filename):
    ##This function takes one argument, the name of an external file. It opens the Fasta file for reading and reads all the lines in one go, returns list of "lines".
    try:
        with open( filename, "r") as fastafile:
            lines = fastafile.readlines();
    except FileNotFoundError:                                       #in case input file refrenced cannot be found
        sys.exit("Sorry, the file does not exist")                  #exits the program after displaying message as it cannot proceed further without file
    return(lines)

def outputfile(lines, z, mtype):
    #This function returns necessary output of sequence, mass, charge and missed cleavages.
    #It accepts three arguments required for calculation- input details, user specified charge and mass type
    output=[["Prot_name", "Peptide","mass-to-charge", "z","p","Sequence"]] #output table headers

    for line in lines:
        if line.startswith('>'):                                     # It expects sequence names, peptide number, missed cleavages to be in the line starting with '>'
            seq_name= line[1:].rstrip('\n').split()[0]               #extracts and stores the sequence name or ID
            mc= line.split("missed=",1)[-1].split()[0]               #extracts and stores the number of missed cleavages
            pno =line.split("peptide ",1)[-1].split()[0]             #extracts and stores peptide number
            continue                                                 #to continue with iterationn before storing values
        elif not line.startswith('>'):                               #It expects peptide sequence to be in the line following the line starting with '>'
            seq = line.split()[0]                                    #extracts and stores peptide sequence

#calls pep2mass function to obtain mass of peptide sequence and stores in output list
        output.append([seq_name, pno, pep2mass(seq,mtype,z), z ,mc, seq]) #output table values
    return(output)

def pep2mass(seq, t, z1):
#This function is used to convert a sequence string, contained in 'seq' to a mass based on type of mass user wants
# It will use one of two possible look-up dictionaries, which is determined by the type m for monoisotropic and a for average mass
#z1 contains charge specified as int value

    mono = {'A' :71.0371, 'C' :103.0092, 'D' :115.0269, 'E' :129.0426, 'F' :147.0684, 'G' :57.0215, 'H' :137.0589, 'I' :113.0841, 'K' :128.0950, 'L' :113.0841, 'M' :131.0405, 'N' :114.0429, 'P' :97.0528, 'Q' :128.0586, 'R' :156.1011, 'S' :87.0320, 'T' :101.0477, 'V' :99.0684, 'W' :186.0793, 'Y' :163.0633, '*' :0.0}

    aver = {'A' :71.08, 'C' :103.14, 'D' :115.09, 'E' :129.12, 'F' :147.18, 'G' :57.05, 'H' :137.14, 'I' :113.16, 'K' :128.17, 'L' :113.16, 'M' :131.19, 'N' :114.10, 'P' :97.12, 'Q' :128.13, 'R' :156.19, 'S' :87.08, 'T' :101.10, 'V' :99.13, 'W' :186.21, 'Y' :163.18, '*' :0.0}
    waterM = 18.0106
    waterA =18.0153
    proton_mass=1.007
    z2=z1*proton_mass                                   #to make mass to charge values of peptides we need to add proton to peptide based on the charge user assigns peptide
    mass=0
    try:
        for i in range(len(seq)):
            if (t =='m'):
                mass = mass + mono[seq[i]]
                pepmass= mass + waterM                  #mass of water is added to peptide equivalent to the terminating groups: H at N terminus and OH at C terminus
                m2c=(pepmass+z2)/z1                     #mass to charge (m+z/z)
            elif (t == 'a'):
                mass = mass + aver[seq[i]]
                pepmass= mass + waterA
                m2c=(pepmass+z2)/z1
    except KeyError:                                     #in case of incorrect or invalid character in sequence not matching with dictionary keys
        sys.exit("Error! Invalid character in Sequence, check sequence in file")
    except ZeroDivisionError:                            #in case charge is 0, it won't be able to divide also 0 charge is not useful for mass spectrometry
        sys.exit("Error! Charge cannot be equal to 0")   #exits the program displaying error message to allow user to correct error
    return round(m2c, 4)                                 #rounds upto 4 decimal places for uniformity

#
# main program
#

if __name__ == "__main__":
    import argparse                                     #argparse aids in user-friendly command line options and also issues errors in case of invalid inputs
    parser = argparse.ArgumentParser(description='convert peptide sequences to masses, read from a Fasta file using mass type specified')
    parser.add_argument("fileName", type = str, help = 'name of fasta file to be read')                            # a positional argument for input file
    group = parser.add_mutually_exclusive_group()       #to ensure either average or monoisotropic mass is used for calculation, they are made mutually exclusive
    group.add_argument('-a',"--average_mass", action='store_true')
    group.add_argument('-m',"--monoisotropic_mass", action='store_true')
    parser.add_argument("charge", type = str, help = 'Specify charge of peptide')                                   # a positional argument for charge user assigns
    parser.add_argument("outfileName", type = str, help = 'name of file output to be written into')                 # a positional argument for output file
    args = parser.parse_args()

    inputfile = args.fileName                           #stores input file name
    c=int(args.charge)                                  #stores charge asigned to peptides
    if args.average_mass:
        masstype='a'
    elif args.monoisotropic_mass:
        masstype='m'
                                                        #stores mass type specified by user

    ifile=[]
    ifile = fastaread(inputfile)                        #stores lines in file as list after reading file

    out=[]
    out=outputfile(ifile,c,masstype)                    #stores output table rows as list returning from function
    outputfile = args.outfileName                       #stores outfile file name

    with open( outputfile, "w") as ofile:
            ofile.write(tabulate(out));                 #tabulates output for clarity and writes into output file before closing it
    ofile.close()
else:
    print("run as module\n")
