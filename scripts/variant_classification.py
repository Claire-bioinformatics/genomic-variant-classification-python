#import statements
import argparse
import logging
import vcf
import os
import gffutils
import math
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import matplotlib.pyplot as plt

#global variables 
count_quality_fail = 0
count_non_coding = 0
count_synonymous = 0
count_non_synonymous = 0

#functions
def cds_coordinate(feature, record, gffdb): #finding the position of the SNP with relation to the CDS
    strand_check = False #forward
    if feature.strand == "-":
        strand_check = True #reverse
    cds_coordinate = 0 
    for child in gffdb.children(feature,order_by = 'start', reverse = strand_check, featuretype="CDS"): #once we found the mRNA where the SNP is, we go and find all the CDS/child in that mRNA in the right order 
        #feature is parent in our main, order_by is to start from the start of mRNA, reverse is for checking the direction, and feature type is for only checking CDS(we don't want exons or whole genes)
        if strand_check == False: #if it is forward strand...
            if record.POS >= child.start and record.POS <= child.end: #you find the CDS that has the SNP this way
                cds_coordinate += record.POS - child.start +1  #if we have found the cds we take the length uptil the SNP
                break
            else: #if the cds with the SNP has not been found yet then we add the whole cds to the count
                cds_coordinate += child.end - child.start +1
        else: #this is reverse strand...
            if record.POS >= child.start and record.POS <= child.end: #you find the CDS that has the SNP this way but from the back
                cds_coordinate += child.end - record.POS +1  #if we have found the cds we take the length uptil the SNP but from the back
                break
            else: #if the cds with the SNP has not been found yet then we add the whole cds to the count
                cds_coordinate += child.end - child.start +1
    return cds_coordinate

def protein_coordinate(cds_coordinate): #calculate the protein coordinate from necleotide/CDS coordinate
    return math.ceil(cds_coordinate/3)                    
                
def synonymous_check(cds_coordinate, protein_coordinate, fastafile, record, feature):
    strand_check = False #i am assuming that we are always on the forward strand
    if feature.strand == "-": #but if there is "-" in GFF file, i am switching to reverse strand
        strand_check = True
    full_sequence = ""
    
    for child in gffdb.children(parent,order_by = "start", reverse = strand_check, featuretype="CDS"): #iterate all the cds of one transcript
        if strand_check == False :
            full_sequence = full_sequence+ child.sequence(fastafile) #get the sequence of cds by adding the iterables
            nucleotide_SNP = str(record.ALT[0]) #give the mutated SNP in a variable
        else: 
            full_sequence = full_sequence+ child.sequence(fastafile, use_strand=True) #reverse strand as our use_strand is complementing
            nucleotide_SNP = str((Seq(str(record.ALT[0]))).complement())#because we want the complement of the base we are storing the sequence as a Seq object before complementing and finally turning it back into a string for further use #also get SNP in complement base
    mutated_sequence = MutableSeq(full_sequence) #putting sequence into an object of the mutableseq class so that we can use function such as translate 
    try:
        mutated_sequence[cds_coordinate - 1] = nucleotide_SNP #mutating the nucleotide sequence from previous line using the SNP
    except IndexError as e:
            logger.error(f"CDS coordinate {cds_coordinate} is out of range: {e}")
    full_sequence = Seq(full_sequence) #putting sequence into an object of the Seq class so that we can use function such as translate 
    original_protein = full_sequence.translate() #translate the nucleotide sequence into protein sequence
    mutated_protein = mutated_sequence.translate()
    if original_protein == mutated_protein: #if SNP has not changed the protein it is synonymous
        check = "Synonymous"
        ref_aa = original_protein[protein_coordinate-1]
        alt_aa = "NA"
    else: #if SNP has caused a change it is non-synonymous
        check = "Non-Synonymous" 
        ref_aa = original_protein[protein_coordinate-1]
        alt_aa = mutated_protein[protein_coordinate-1]

    return check,ref_aa,alt_aa #returning check:synonymous or non-synonymous, ref amino acid, alt amino acid
    

#argparser for making a command line interface

#creates our argparse handler
parser = argparse.ArgumentParser(description="To classify SNP data using GFF, VCF and FASTA files.") 

#creating arguments for CLI
parser.add_argument('--vcffile', required=True, help="This stores the name of the VCF file")
parser.add_argument('--gfffile', required=True,help="This stores the name of the GFF file")
parser.add_argument('--fastafile',required=True, help="This stores the name of the FASTA file")
args = parser.parse_args() #runs the parser

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.FileHandler('3053994_log_file.txt')
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)


#checking if all the files are opening
list_of_files = [args.vcffile,args.gfffile,args.fastafile]#making this variable so that i don't have to write 'try and except' 3 times for 3 different files
flag = False #setting my error flag to false 
try:
    for file in list_of_files: #trying to open and find the files
        with open(file, 'r') as f:
            pass
except FileNotFoundError: #if the file isn't found we set the error flag to true
    logger.error(f"The file {file} was not found")
    flag = True

if flag == True: #if error was found we exit the program
    logger.error("Please check the filenames given at command line - code could not run because files weren't found")
    raise SystemExit(1)
    
logger.info(f"file names given at the command line were: {args.vcffile}, {args.gfffile}, {args.fastafile}\n")

# open the vcf file
vcfReader = vcf.Reader(filename=args.vcffile)

# make or connect to gff database
gff_db = args.gfffile.replace('.gff', '.db') # generate the db filename by replacing the '.gff' extension in the GFF file with '.db'. This will be the name of the GFF db.
if not os.path.isfile(gff_db): # check if the db file already exists.
    logger.info(f'Creating GFF database now {gff_db}...\n')
    gffdb = gffutils.create_db(args.gfffile, dbfn=gff_db, force=True, keep_order=True) # create a new GFF db from the input GFF file.#keep_order = true , keeps the column order of the db the same as gff file
else:
    logger.info(f'Connecting to existing database now {gff_db}...\n') #if the db already exists...
    try:
        gffdb = gffutils.FeatureDB(gff_db, keep_order=True) # connect to the existing GFF db file.
    except ValueError: # if an error occurs while connecting to the db:
        logger.error(f'DB {gff_db} could not be read please try again \n')
        raise SystemExit(1) # Exit the program with a status code of 1, indicating an error

output_table = open("3053994_table.tsv", "w") # open a new file called "3053994_table.tsv" in write mode to store the output table.
output_table.write(f"Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n") # write the header row of the output table, defining the columns: 

try:
    for record in vcfReader: #iterating through VCF file
        if record.QUAL > 20: #finding the records that passed QC
            if record.REF != record.ALT: #finding the records that have SNP change
                for feature in gffdb.region(seqid=record.CHROM, start= record.POS, featuretype='CDS'): #finding all the CDSs that have SNP change
                    if record.POS >= feature.start and record.POS <= feature.end: #if the SNP change is in CDS region
                        for parent in gffdb.parents(feature): #then we want to find all the mRNA of that CDS
                            transcript_id = parent.id
                            cds_coordinate_value = cds_coordinate(parent,record,gffdb) #we want to find the position of SNP in relation to CDS #go see comments of function cds_coordinate
                            protein_coordinate_value = protein_coordinate(cds_coordinate_value) #dividing by 3 to convert from nucleotide to protein length #see comments of function protein_coordinate
                            check,ref_aa,alt_aa = synonymous_check(cds_coordinate_value,protein_coordinate_value,args.fastafile,record,parent) #supplying cds coordinate (SNP position in relation to CDS), protein coordinate (SNP postion in relation to protein), fasta file that has the neucleotide sequences, record which is a row of VCF file, parent which is one transcript that has SNP in GFF file, to do the synonymous check #see synonymous check function for comments
                            
                            if check == "Synonymous": #count synonymous
                                count_synonymous+=1
                            else:
                                count_non_synonymous+=1 #if not synonymous then count it as non-synonmymous
                            break # Exit the parent loop since the SNP has been processed for the current CDS.
                    else: # if the SNP is not in a coding region:
                        count_non_coding += 1
                        check = "Non Coding" #its classification
                        transcript_id = "NA"
                        protein_coordinate_value = "NA"
                        ref_aa = "NA"
                        alt_aa = "NA"
                        break # exit the loop since the SNP has been classified as non-coding.
            output_table.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{check}\t{transcript_id}\t{protein_coordinate_value}\t{ref_aa}\t{alt_aa}\n") ## Write the SNP information to the output table, including its classification and details.
        else: # if the variant's quality score is less than or equal to 20:
            count_quality_fail += 1
except Exception as e: # handle any unexpected exceptions during the processing loop.
    logger.warning(f"An unexpected exception has occured and this is the code: {e}") # log a warning with the exception details.


logger.info(f'number of SNP records that failed the quality check = {count_quality_fail}') # log the total number of SNP records that failed the quality check.

#To generate the barplot
variables= ["Non-coding", "Synonymous", "Non-synonymous"]
values= [count_non_coding, count_synonymous, count_non_synonymous]

plt.bar(variables, values)
plt.title('Counts of what SNP with QUAL>20 are classified as')
plt.xlabel('Types')
plt.ylabel('Counts')
plt.savefig("3053994_barplot.png")

logging.shutdown() # gracefully shut down the logging system
output_table.close() #close the tsv file

logger.info(f"the output files are: 3053994_table.tsv, 3053994_log_file.txt, 3053994_barplot.png. They are located at {os.getcwd()}\n")
# log an informational message specifying the names of the output files generated by the script
# and their location in the current working directory (obtained using `os.getcwd()`).
