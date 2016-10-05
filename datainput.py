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

# for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
#     result_handle = NCBIWWW.qblast("blastn","Bacteria and Archaea", seq_record)
#     blast_result = open("my_blast.xml", "w")
#     blast_result.write(result_handle.read())
#     blast_result.close()
#     result_handle.close()

    #
    # result_handle = NCBIWWW.qblast("blastn", "Bacteria and Archaea", seq_record)

# blast_records_bacteria = NCBIXML.parse(result_handle)
# blast_record = NCBIXML.read(result_handle)
# E_VALUE_THRESH = 0.0001

# for alignment in blast_record.alignments:
#     for hsp in alignment.hsps:
#         if hsp.expect > E_VALUE_THRESH:
#             print('****Alignment****')
#             print('sequence:', alignment.title)
#             print('length:', alignment.length)
#             print('e value:', hsp.expect)
#
# file_string = ""
x = 1

bacteriacontigs = []
def add(i, s):
    size = len(bacteriacontigs)
    if i >= size:
        bacteriacontigs.extend([None]*(i-size+1))
        bacteriacontigs[i] = s


for seq_record in SeqIO.parse("largecontigs.fasta", "fasta"):
    #File = open("/Users/Virginiasaulnier/Downloads/results" + str(x) + ".txt", "w")
    #fasta_string = open(i+"/Users/Virginiasaulnier/Downloads/largecontigs.fasta").read() #or make the names fasta1.fasta and just do open(i).read
    result_handle = NCBIWWW.qblast("blastn","nt", seq_record.seq, entrez_query="bacteria[organism]")
    print('in seq record')
    #blast_results = result_handle.read()
    #blast_records = list(blast_results)

    # Next, we save this string in a file:
    # save_file = open("/Users/virginiasaulnier/downloads/my_blast_bacteria"+ str(x) + ".xml", "w")
    # x += 1
    # save_file.write(str(blast_records))
    # save_file.close()

    blast_records = NCBIXML.parse(result_handle)
#     print('went through result parse')
#
#     #or blast_record = NCBIXML.read(result_handle) if you only have one seq in file

    E_VALUE_THRESH = 0.0001
    for blast_record in blast_records:
         print('went into for blast_record')
         for alignment in blast_record.alignments:
            print('in alignment loop')
            for hsp in alignment.hsps:
                 print('in hsp loop')
                 if hsp.expect < E_VALUE_THRESH:
                     print('in expect loop')
                     # seq_record.id += "alignment:",alignment.title+"\n"
#                     # seq_record.id += "e-value:",hsp.expect+"\n"
                     print('****Alignment****')
                     print(seq_record.id)
                     print('sequence:', alignment.title)
                     print('length:', alignment.length)
                     print('e value:', hsp.expect)
                     # save_file = open("/Users/virginiasaulnier/downloads/my_blast_bacteria"+ ".txt", "w")
                     # save_file.write(str(blast_record))
                     bacteriacontigs.insert(x, seq_record.id)
                     #save_file.close()
                     x += 1
#                     File.write(file_string)
#
#
# blast_records = NCBIXML.parse(result_handle)
#save_file.close()



#                 # Get a BLAST record
# blast_record = blast_records.next()
# for alignment in blast_record.alignments:
#     print('in alignment loop 2nd')
#     for hsp in alignment.hsps:
#         print('in hsp loop 2nd')
#         if hsp.expect < E_VALUE_THRESH:
#             print('in expect loop 2nd')
#                 # seq_record.id += "alignment:",alignment.title+"\n"
#                 # seq_record.id += "e-value:",hsp.expect+"\n"
#             print('****Alignment****')
#             print(seq_record.id)
#             print('sequence:', alignment.title)
#             print('length:', alignment.length)
#             print('e value:', hsp.expect)
#
#             x += 1
#             File.write(file_string)
#
#
# #bacteria_results = (rec for rec in SeqIO.parse(file, "fasta") if len(rec) > 2000)
#
# # blast_records_bacteria = NCBIXML.parse(result_handle)
# # blast_record = NCBIXML.read(result_handle)
# # E_VALUE_THRESH = 0.0001
# #
# # for alignment in blast_record.alignments:
# #     for hsp in alignment.hsps:
# #         if hsp.expect > E_VALUE_THRESH:
# #             print('****Alignment****')
# #             print('sequence:', alignment.title)
# #             print('length:', alignment.length)
# #             print('e value:', hsp.expect)
