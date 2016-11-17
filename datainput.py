from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logging.debug('This is a log message.')
file="/Users/Virginiasaulnier/Downloads/v1S2tiny.consensus"

largecontigs= (rec for rec in SeqIO.parse(file,"fasta") if len(rec)>2000)
SeqIO.write(largecontigs, "largecontigs.fasta", "fasta")


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

for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    result_handle = NCBIWWW.qblast("blastn","nt", seq_record.seq, entrez_query="bacteria[organism]", expect=.0001)
    blast_records = NCBIXML.parse(result_handle)
    logging.debug('This is a log message after parse blast_records.')

    for record in blast_records:
        i=0
        for alignment in record.alignments:
            while(i == 0):
                bacteriacontigs.append(seq_record.id)
                i = i+1
        logging.debug('This is a log message appending bacteria seq_records.')
        print(seq_record)
print(bacteriacontigs)

#print(bacteriacontigs)

for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in bacteriacontigs:
         logging.debug('This is a log message in virus blast.')
         result_handle_virus = NCBIWWW.qblast("blastn", "nt", seq_record.seq, entrez_query="viruses[organism] AND perc_ident[1:30]", expect=.0001)
         blast_records_virus = NCBIXML.parse(result_handle_virus)

         for record in blast_records_virus:
             i=0
             for alignment in record.alignments:
                 while(i==0):
                    viralcontigs.append(seq_record.id)
                    i = i +1

print(viralcontigs)



for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in bacteriacontigs:
        logging.debug('This is a log message in viral protein blast.')
        result_handle_viralprotein=NCBIWWW.qblast("blastx", "pdb", seq_record.seq, entrez_query="viruses[organism] AND perc_ident[1:30]", expect=.0001)
        blast_records_viralprotein=NCBIXML.parse(result_handle_viralprotein)

        for record in blast_records_viralprotein:
            i=0
            for alignment in record.alignments:
                while(i==0):
                    viralcontigs.append(seq_record.id)
                    i = i +1

print(viralcontigs)

clusterids=set(bacteriacontigs) - set(viralcontigs)

print(clusterids)

for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    if seq_record.id in bacteriacontigs:
        while(i==0):
            clustercontigs.append(seq_record)
            i=i+1

print(clustercontigs)
