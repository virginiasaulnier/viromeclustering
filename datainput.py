from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
#somethingchanged
file="/Users/Virginiasaulnier/Downloads/v1S2short.consensus"

largecontigs= (rec for rec in SeqIO.parse(file,"fasta") if len(rec)>2000)
SeqIO.write(largecontigs, "largecontigs.fasta", "fasta")
for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


x = 1

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


for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    #File = open("/Users/Virginiasaulnier/Downloads/results" + str(x) + ".txt", "w")
    #fasta_string = open(i+"/Users/Virginiasaulnier/Downloads/largecontigs.fasta").read() #or make the names fasta1.fasta and just do open(i).read
    result_handle = NCBIWWW.qblast("blastn","nt", seq_record.seq, entrez_query="bacteria[organism]", expect=.0001)
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        bacteriacontigs.append(seq_record.id)
    print('in seq record')


#print(bacteriacontigs)

for record in blast_records:
    if seq_record.id in bacteriacontigs:
         result_handle_virus = NCBIWWW.qblast("blastn", "nt", seq_record.seq, entrez_query="viruses[organism]" and "1:30[perc_ident]",
                                           expect=.0001)
         blast_records_virus = NCBIXML.parse(result_handle_virus)

         for record in blast_records_virus:
            viralcontigs.append(seq_record.id)
         print(viralcontigs)


#
# for record in blast_records:
#     if seq_record.id in bacteriacontigs:
#         result_handle_viralprotein=NCBIWWW.qblast("blastx", "pdb", seq_record.seq, entrez_query="viruses[organism]"and "1:30[perc_ident]",
#                                            expect=.0001)
#         blast_records_viralprotein = NCBIXML.parse(result_handle_viralprotein)
#
#         for record in blast_records_viralprotein:
#             viralcontigs.append(seq_record.id)
#         print(viralcontigs)

