import subprocess
import subprocess
import numpy as np
from Bio import SeqIO, AlignIO
import re

def most_frequent_string(lst):
    d = {}
    for i in lst:
        if i in d:
            d[i] += 1
        else:
            d[i] = 1
    return max(d, key=d.get)

def get_position(enst):
    obj = subprocess.Popen(
    "grep $'"+enst+"'\t ./downloads/Ensembl_Homo_sapiens.GRCh38.110.bed",
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
    )
    enst_pos = obj.stdout.read() 
    return enst_pos.decode('utf-8').strip()

def get_seq(enst):
    obj = subprocess.Popen(
    'grep '+enst+' ./downloads/Homo_sapiens.GRCh38.cds.all.tsv',
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
    )
    seq = obj.stdout.read() 
    return seq.decode('utf-8').strip()

def to_onehot(seq):
    atgc_dict = {'A':0, 'T':1, 'G':2, 'C':3}
    M = np.zeros((max_len, 4))
    for i, s in enumerate(seq):
        M[i, atgc_dict[s.upper()]] = 1
    return M

def onehot_to_seq(M):
    atgc_dict = {0:'A', 1:'T', 2:'G', 3:'C'}
    seq = ""
    for i in range(M.shape[0]):
        seq += atgc_dict[np.argmax(M[i,:])]
    return seq


def aln_clustal():
    obj = subprocess.Popen(
    'clustalo -i ./downloads/tmp_for_aligning.fa > ./downloads/alignment.fa',
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
    ).wait()

    alignment = AlignIO.read("alignment.fa", "fasta")

def getConsensus(alignment):
    consensus = ""
    for i in range(len(alignment[0])):
        base_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
        for record in alignment:
            if i < len(record):
                base = record[i]
                base_counts[base] += 1
            else:
                base_counts['-'] += 1
        if max(base_counts, key=base_counts.get) == '-':
            del base_counts['-']
            consensus += max(base_counts, key=base_counts.get)    
        consensus += max(base_counts, key=base_counts.get)
    
    return consensus


def multiAlign(sequenced_seqs, ref_seq=None, ref_seq_top5=None):
    new_sequenced_seqs, flags = [], []
    if ref_seq == None:
        ref_seq = sequenced_seqs[0]
    int_flag = np.zeros((len(sequenced_seqs), len(ref_seq)))

    for idx,seq in enumerate(sequenced_seqs):
        if ref_seq_top5 != None:
            start = seq.find(ref_seq_top5)
        else:
            start = 0
        seq = seq[start:]
        new_seq = ""
        new_flag = ""
        len_seq = len(seq)
        len_ref_seq = len(ref_seq)
        i, j = 0, 0
        while(i<len_seq and j<len_ref_seq):
            if ref_seq[j] == seq[i]:
                new_seq += seq[i]
                new_flag += "="

            elif (i+1<len_seq and j+2<len_ref_seq) and seq[i] == ref_seq[j+1] and seq[i+1] == ref_seq[j+2]:
                new_seq += '-'
                new_seq += seq[i]
                new_seq += seq[i+1]
                new_flag += "==="
                i += 1
                j += 2
            elif (i+1<len_seq and j+3<len_ref_seq) and seq[i] == ref_seq[j+2] and seq[i+1] == ref_seq[j+3]:
                new_seq += '--'
                new_seq += seq[i]
                new_seq += seq[i+1]
                new_flag += "===="
                i += 1
                j += 3
            elif (i+1<len_seq and j+4<len_ref_seq) and  seq[i] == ref_seq[j+3] and seq[i+1] == ref_seq[j+4]:
                new_seq += '---'
                new_seq += seq[i]
                new_seq += seq[i+1]
                new_flag += "====="
                i += 1
                j += 4
            else:
                new_seq += seq[i]
                new_flag += "X"
                int_flag[idx, j] = 1
                
            i += 1
            j += 1
            
        new_sequenced_seqs.append(new_seq)
        flags.append(new_flag)
    return new_sequenced_seqs, flags, int_flag

def getEnst(geneName1):
    command = f"""grep -E 'transcript.*gene_name "{geneName1}"' ./downloads/Homo_sapiens.GRCh38.110.gtf | awk '$3 == "exon" {{split($0, a, "; "); for (i in a) if (match(a[i], /^exon_number "[0-9]+"/)) {{split(a[i], b, " "); if (b[2] > max) {{max=b[2]; line=$0}}}}}} END {{print line}}'
    """
    print(command)
    proc = subprocess.Popen(
        command, 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True)
    outinfo, errinfo = proc.communicate()
    return outinfo


in_file = "./downloads/example.reads.fa"
ratio_cutoff = 0.5
out_file = "./downloads/example.reads.fates.tsv"

name_to_seq = {}
name_to_enst = {}
id_1, id_2, enst_1, name_1, enst_2, name_2, concat_name, concat_name_r = "", "", "", "", "", "", "", ""
full_id_list = []
NMD_or_PROTEIN = {}
idx = 0


ENST_MODE = True
'''
ENST_MODE,
>ENST00000541412,TMEM14C|ENST00000481240,TMEM14B|6|1|5fd55169-79d9-4755-b60d-ebe7d51bec93
TACGCTGATATTGCTGGGGAGAGACTGGCTGCTGTG...
otherwise,
>Novel3_177_perfect_19584_F_0_2803_0,ZBED1,DHRSX
ACCTCCAACCTGTCCTACCACCTGGAGAAGAACCACCCCGAGGAATTCTGCGAGTTCGTCAAGAGCAACACGGAGCAGATGCGTGAAGCCTTCGCCACCGCCTTCTCCAAGCTGAAGCCCGAGTCGTCCCAGCAGCCCGGGCAGGACGCGCTGGCCGTCAAGGCCGGCCACGGCTACGACAGCAAGAAGCAGCAGGAGCTGACGGCCGCCGTGCTGGGCCTCATCTGCGAGGGGCTGTACCCAGCCTCCATCGTGGACGAGCCCACCTTCAAGGTGCTGCTGAAGACGGCCGACCCCCGGTATGAGCTGCCCAGCCGGAAGTACATCTCTACCAAGGCCATCCCTGAGAAGTACGGGGCCGTCCGGGAGGTGATCCTGAAGGAGCTGGCCGAGGCCACCTGGTGTGGCATCTCCACCGACATGTGGAGGAGTGAGAATCAGAACCGCGCCTACGTCACGCTGGCCGCCCACTTCCTGGGCCTGGGCGCCCCCAACTGCCTGTCCATGGGCTCCCGCTGCCTGAAGACCTTCGAGGTGCCCGAAGAGAACACGGCGGAGACCATCACGCGAGTGCTCTATGAGGTCTTCATCGAGTGGGGCATCAGCGCCAAGGTCTTCGGGGCCACCACCAACTATGGCAAGGACATCGTGAAGGCGTGCTCCCTGCTGGACGTCGCAGTGCACATGCCCTGCCTGGGCCACACCTTCAATGCCGGCATCCAGCAGGCCTTCCAGCTCCCGAAGCTGGGGGCGCTGCTGTCGCGCTGCCGCAAACTGGTGGAGTACTTCCAGCAGTCTGCCGTGGCCATGTACATGCTCTATGAGAAGCAGAAGCAGCAGAACGTGGCCCACTGCATGCTGGTGAGCAACCGCGTCTCCTGGTGGGGGAGCACGCTGGCCATGCTGCAGCGCCTCAAGGAGCAGCAGTTCGTCATCGCCGGGGTCTTGGTGGAGGACAGCAACAACCACCACCTCATGCTGGAGGCCAGCGAGTGGGCCACCATCGAGGGGCTGGTGGAGCTCCTGCAGCCCTTCAAGCAGGTGGCCGAGATGCTGTCGGCCTCCAGGTACCCCACCATCAGCATGGTGAAGCCGCTGCTGCACATGCTCCTGAACACCACGCTCAACATCAAGGAGACCGACTCCAAGGAGCTCAGCATGGCCAAGGAGGTCATCGCCAAGGAGCTTTCCAAGACCTACCAGGAGACGCCCGAGATCGACATGTTTCTCAACGTGGCCACCTTCCTGGACCCCCGCTACAAGAGGCTGCCCTTCCTCTCCGCCTTCGAGCGGCAGCAGGTGGAGAATCGCGTGGTGGAAGAGGCCAAGGGCCTGCTGGACAAGGTCAAAGACGGCGGCTACCGGCCGGCTGAGGACAAGATCTTCCCGGTGCCCGAGGAGCCTCCCGTCAAGAAGCTCATGCGGACATCCACGCCGCCGCCCGCCAGCGTCATCAACAACATGCTGGCCGAGATCTTCTGCCAGACAGGCGGCGTGGAGGACCAGGAAGAGTGGCATGCCCAGGTGGTGGAGGAGCTGAGCAACTTCAAGTCCCAGAAGGTGCTTGGCCTCAACGAAGACCCCCTCAAGTGGTGGTCAGACCGCCTGGCCCTCTTCCCCCTGCTGCCCAAGGTGCTGCAGAAGTACTGGTGCGTGACGGCCACGCGCGTCGCCCCTGAGCGTCTCTTCGGATCCGCCGCCAACGTGGTCAGCGCCAAGAGGAACCGGCTGGCTCCCGCGCACGTGGACGAGCAGGTGTTTCTGTATGAGAACGCCCGGAGTGGGGCAGAGGCGGAACCCGAGGACCAGGACGAGGGGGAGTGGGGCCTGGACCAGGAGCAGGTGTTCTCCTTGGGGGATGGCGTCAGCGGCGGTTTCTTTGGCATTAGGGACAGCAGCTTCCTGTAGATGTCGCCATTGTCTGCGGCGCGGGCGGCCCTGCGGGTCTACGCGGTAGGCGCCGCGGTGATCCTGGCGCAGCTGCTGCGGCGCTGCCGCGGGGGCTTCCTGGAGCCAGTTTTCCCCCCACGACCTGACCGTGTCGCTATAGTGACGGGAGGGACAGATGGCATTGGCTATTCTACAGCGAAGCATCTGGCGAGACTTGGCATGCATGTTATCATAGCTGGAAATAATGACAGCAAAGCCAAACAAGTTGTAAGCAAAATAAAAGAAGAAACCTTGAACGACAAAGTGGAATTTTTATACTGTGACTTGGCTTCCATGACTTCCATCCGGCAGTTTGTGCAGAAGTTCAAGATGAAGAAGATTCCTCTCCATGTCCTGATCAACAATGCTGGGGTGATGATGGTCCCTCAGAGGAAAACCAGAGATGGATTCGAAGAACATTTCGGCCTGAACTACCTAGGGCACTTCCTGCTGACCAACCTTCTCTTGGATACGCTGAAAGAGTCTGGGTCCCCTGGCCACAGTGCGAGGGTGGTCACCGTCTCCTCTGCCACCCATTACGTCGCTGAGCTGAACATGGATGACCTTCAGAGCAGTGCCTGCTACTCACCCCACGCAGCCTACGCCCAGAGCAAGCTGGCCCTTGTCCTGTTCACCTACCACCTCCAGCGGCTGCTGGCGGCTGAGGGAAGCCACGTGACCGCCAACGTGGTGGACCCCGGGGTGGTCAACACGGACGTCTACAAGCACGTGTTCTGGGCCACCCGTCTGGCGAAGAAGCTTCTCGGCTGGTTGCTTTTCAAGACCCCCGATGAAGGAGCGTGGACTTCCATCTACGCAGCAGTCACCCCAGAGCTGGAAGGAGTTGGTGGCCATTACCTATACAACGAGAAAG
'''
with open(in_file) as f:
    for q in f:
        if '>' in q:
            full_id_list.append(q.rstrip()[1:])
            if ENST_MODE:
                id_1 = q.rstrip()[1:].split("|")[0]
                id_2 = q.rstrip()[1:].split("|")[1]
                enst_1, name_1 = id_1.split(',')[0].split('.')[0], id_1.split(',')[1]
                enst_2, name_2 = id_2.split(',')[0].split('.')[0], id_2.split(',')[1]
                concat_name = name_1 + '_' + name_2
                concat_name_r = name_1 + '_' + name_2
                if concat_name in name_to_enst.keys():
                    name_to_enst[concat_name].append((enst_1, enst_2))
                elif concat_name_r in name_to_enst.keys():
                    name_to_enst[concat_name_r].append((enst_2, enst_1))
                else:
                    name_to_enst[concat_name] = [(enst_1, enst_2)]

            else:
                name_1 = q.rstrip().split(',')[1]
                name_2 = q.rstrip().split(',')[2]
                concat_name = name_1 + '_' + name_2
                concat_name_r = name_1 + '_' + name_2
                if concat_name not in name_to_enst.keys() and concat_name_r not in name_to_enst.keys():
                    outinfo = getEnst(name_1).decode('gbk')
                    matches = re.findall(r'transcript_id "(.*?)"', outinfo)
                    enst1 = matches[0]
                    outinfo = getEnst(name_2).decode('gbk')
                    matches = re.findall(r'transcript_id "(.*?)"', outinfo)
                    enst2 = matches[0]
                    name_to_enst[concat_name] = [(enst_1, enst_2)]
                else:
                    pass
        else:
            seq = q.rstrip()
            if concat_name in name_to_seq.keys():
                name_to_seq[concat_name].append(seq)
            elif concat_name_r in name_to_seq.keys():
                name_to_seq[concat_name_r].append(seq)
            else:
                name_to_seq[concat_name] = [seq]

record_idx = 0


for name, ensts_1_and_2 in name_to_enst.items():

    ensts_1_and_2 = np.array(ensts_1_and_2)
    sequenced_seqs = name_to_seq[name]
    seqs_num = len(sequenced_seqs)
    
    enst1 = most_frequent_string(ensts_1_and_2[:,0])
    enst2 = most_frequent_string(ensts_1_and_2[:,1])


    enst1_pos = get_position(enst1).split('\t')
    enst2_pos = get_position(enst2).split('\t')
    enst1_start, enst1_end, strand = enst1_pos[1], enst1_pos[2], enst1_pos[5]
    enst2_start, enst2_end = enst2_pos[1], enst2_pos[2]


    if strand == '-1' or strand == '-':
        if enst1_end > enst2_end:
            enst_upstream = enst1 
            enst_downstream = enst2
        else:
            enst_upstream = enst2 
            enst_downstream = enst1
    elif strand == '1' or strand == '+':
        if enst1_start < enst2_start:
            enst_upstream = enst1 
            enst_downstream = enst2
        else:
            enst_upstream = enst2 
            enst_downstream = enst1


    if enst_upstream == enst1:
        start = int(enst2_pos[1])
        thick_start = int(enst2_pos[6])
        sizes = [int(i) for i in enst2_pos[-2].split(',')[:-1]]
        shift = [int(i) for i in enst2_pos[-1].split(',')[:-1]]
        juncs = [start+i+j for i,j in zip(shift, sizes)]
        thick_end = int(enst2_pos[7])
        flag = True
        for idx,i in enumerate(juncs):
            if flag and i > thick_start:
                rel = int(np.sum(sizes[:idx])) + thick_start - (start+shift[idx])
                flag = False
            if i > thick_end:
                exon_exon_junc = np.sum(sizes[:idx]) - rel
                break
    else:
        start = int(enst1_pos[1])
        thick_start = int(enst1_pos[6])
        sizes = [int(i) for i in enst1_pos[-2].split(',')[:-1]]
        shift = [int(i) for i in enst1_pos[-1].split(',')[:-1]]
        juncs = [start+i+j for i,j in zip(shift, sizes)]
        thick_end = int(enst1_pos[7])
        flag = True
        for idx,i in enumerate(juncs):
            if flag and i > thick_start:
                rel = int(np.sum(sizes[:idx])) + thick_start - (start+shift[idx])
                flag = False
            if i > thick_end:
                exon_exon_junc = int(np.sum(sizes[:idx])) - rel
                break
    

    ref_seq = get_seq(enst_upstream).split('\t')[1]
    ref_seq_top5 = ref_seq[:5]
    ref_seq_downstream = get_seq(enst_downstream).split('\t')[1]
    try:
        ref_seq_downstream_top5 = ref_seq_downstream[exon_exon_junc:exon_exon_junc+5]
    except:
        print(exon_exon_junc)
        pass

    cut_sequenced_seqs = []
    for query_seq in sequenced_seqs:
        last_p = query_seq.find(ref_seq_downstream_top5)
        if last_p == -1:
            cut_sequenced_seqs.append(query_seq)
        else:
            cut_sequenced_seqs.append(query_seq[:last_p])


    ref_len = len(ref_seq)
    new_sequenced_seqs,flags, int_flag = multiAlign(cut_sequenced_seqs, ref_seq, ref_seq_top5)
    

    try:
        new_sequenced_seqs.append(ref_seq)
        
        sum_of_flag = np.sum(int_flag, axis=0)
        inserted_exon_start = 0
        for i in range(sum_of_flag.shape[0]):
            if sum_of_flag[i] > int_flag.shape[0]*0.6 and sum_of_flag[i+1] > int_flag.shape[0]*0.6:
                inserted_exon_start = i - i % 3
                break
        new_new_sequenced_seqs = []
        for seq, new_seq in zip(cut_sequenced_seqs, new_sequenced_seqs):
            actual_start = inserted_exon_start - new_seq[:inserted_exon_start].count('-')
            new_new_sequenced_seqs.append(seq[actual_start:])
        

        STOP_CODONS = ['TAG','TAA','TGA']
        truncted_seqs,_,_ = multiAlign(new_new_sequenced_seqs)
        
        corrected_seq = getConsensus(truncted_seqs)
        truncated_length = 55
        for i in range(0, len(corrected_seq), 3):
            if corrected_seq[i:i+3] in STOP_CODONS:
                if len(corrected_seq)-i > truncated_length:
                    NMD_or_PROTEIN[record_idx] = ("NMD", corrected_seq, len(corrected_seq)-i)
                else:
                    NMD_or_PROTEIN[record_idx] = ("PROTEIN", corrected_seq, len(corrected_seq)-i)
                # record_idx += 1
                # break
        record_idx += 1
    except:
        pass
    
with open(out_file, 'w') as w:
    for id, value in NMD_or_PROTEIN.items():
        prediction,seq,_ = value
        id1, id2= full_id_list[id].split('|')[0], full_id_list[id].split('|')[1]
        id1, id2 = id1.split(',')[1], id2.split(',')[1]
        print('>'+id1+','+id2+','+prediction+'\n'+seq, file=w)
        