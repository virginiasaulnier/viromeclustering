import logging
import os

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
#this is log output, to help me track where the time issues were, and to make sure all loops were working
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logging.debug('This is a log message.')

#input file
file="/Users/Virginiasaulnier/Downloads/v1S2.consensus"

#sort the contigs by size and put the filtered contigs in a separate file
largecontigs= (rec for rec in SeqIO.parse(file,"fasta") if len(rec)>2000)
SeqIO.write(largecontigs, "largecontigs.fasta", "fasta")

#creating 3 arrays to store contig ids of interest
bacteriacontigs = []
def add(i, s):
    size = len(bacteriacontigs)
    if i >= size:
        bacteriacontigs.extend([None]*(i-size+1))
        bacteriacontigs[i] = s

viralcontigs = []
def add(i, s):
    size = len(viralcontigs)
    if i >= size:
        viralcontigs.extend([None]*(i-size+1))
        viralcontigs[i] = s


clustercontigs = []
def add(i, s):
    size = len(viralcontigs)
    if i >= size:
        viralcontigs.extend([None]*(i-size+1))
        viralcontigs[i] = s

#parse contigs from filtered(by size) selection and blast against all bacteria with an e-value less than .0001
for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    result_handle = NCBIWWW.qblast("blastn","nt", seq_record.seq, entrez_query="bacteria[organism]", expect=.0001)
    blast_records = NCBIXML.parse(result_handle)
    logging.debug('This is a log message after parse blast_records.')

#for each blast result(record) of each contig, if the result has at least one alignment
    #  add the sequence record id to an array called bacteriacontigs (only once, then leave loop)
    for record in blast_records:
        i=0
        for alignment in record.alignments:
            while(i == 0):
                bacteriacontigs.append(seq_record.id)
                i = i+1
        logging.debug('This is a log message appending bacteria seq_records.')
        print(seq_record)
print(bacteriacontigs)


#again go through the contigs of appropriate size, if contig also has a matching bacterial alignment(id in bacteriacontigs array)
#blast against all viruses with a percent identity of 1-30% and an e-value of < .0001
#percent identity of 1-30% means if we are looking at a sequence length 1000, we only want matches of 10-300 nucelotides
for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in bacteriacontigs:
         logging.debug('This is a log message in virus blast.')
         result_handle_virus = NCBIWWW.qblast("blastn", "nt", seq_record.seq, entrez_query="viruses[organism] AND perc_ident[1:30]", expect=.0001)
         blast_records_virus = NCBIXML.parse(result_handle_virus)

         # for each blast result(record) of each contig, if the result has at least one alignment
         #  add the sequence record id to an array called viralcontigs (only once, then leave loop)
         for record in blast_records_virus:
             i=0
             for alignment in record.alignments:
                 while(i==0):
                    viralcontigs.append(seq_record.id)
                    i = i +1

print(viralcontigs)


#repeat same procedure as above, but this time using blastx against viral protein database
for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in bacteriacontigs:
        logging.debug('This is a log message in viral protein blast.')
        result_handle_viralprotein=NCBIWWW.qblast("blastx", "pdb", seq_record.seq, entrez_query="viruses[organism] AND perc_ident[1:30]", expect=.0001)
        blast_records_viralprotein=NCBIXML.parse(result_handle_viralprotein)

        #add sequence ID to viralcontig array
        for record in blast_records_viralprotein:
            i=0
            for alignment in record.alignments:
                while(i==0):
                    viralcontigs.append(seq_record.id)
                    i = i +1

print(viralcontigs)

#for all the contigs that had results in both the bacteria and viral blasts, remove those contigs.
clusterids=set(bacteriacontigs) - set(viralcontigs)

#clusterids should now contain only contigs IDs that had alignments in bacteria, and nothing in virus.
print(clusterids)

#open the filtered by size contigs, if the contig ID is listed in clusterids array,
#add the entire sequence record to the array
for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in clusterids:
        clustercontigs.append(seq_record)

#write these records to a file called clustercontigs.fasta
print list(clustercontigs)
SeqIO.write(clustercontigs, "clustercontigs.fasta", "fasta")

#using usearch cluster contigs, outputs are tree.phy, and clusters.txt
retvalue= os.popen("/users/virginiasaulnier/Downloads/usearch9.1.13_i86osx32 -cluster_agg clustercontigs.fasta -treeout tree.phy -clusterout clusters.txt -id 0.80 -linkage min").readlines()
print retvalue