import matplotlib.pyplot as plt
import subprocess
import matplotlib
import re

def grep_bed_or_fasta(enst=None, mode=True, strand=None):
    if mode:
        enst = enst.split('.')[0]
        proc = subprocess.Popen(
        f"grep {enst} ./downloads/Ensembl_Homo_sapiens.GRCh38.110.bed", 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True)
        outinfo, errinfo = proc.communicate()
        return outinfo
    else:
        cmd = "./downloads/bin/bedtools getfasta -fi ./downloads/hg38.fa -bed ./downloads/tmp.bed > ./downloads/seq1.fa" if strand=='+' \
            else "./downloads/bin/bedtools getfasta -fi ./downloads/hg38.fa -bed ./downloads/tmp.bed -s > ./downloads/seq1.fa"
        proc = subprocess.Popen(
        cmd, 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True)
        outinfo, errinfo = proc.communicate()
        return None

def getEnst(geneName1):
    command = f"""grep -E 'transcript.*gene_name "{geneName1}"' ./downloads/Homo_sapiens.GRCh38.110.gtf | awk '$3 == "exon" {{split($0, a, "; "); for (i in a) if (match(a[i], /^exon_number "[0-9]+"/)) {{split(a[i], b, " "); if (b[2] > max) {{max=b[2]; line=$0}}}}}} END {{print line}}'
    """
    proc = subprocess.Popen(
        command, 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True)
    outinfo, errinfo = proc.communicate()
    return outinfo


def to_str(ori_list):
    return [str(i) for i in ori_list]


def Plot(start, end, chr_id):

    start, end = start, end
    SIZE_SCALE = abs(end - start) // 7
    HALF_SIZE_SCALE = SIZE_SCALE // 2

    fig = plt.figure(figsize=(16, 3))
    ax = fig.add_subplot(111, aspect='equal')
    plt.xlim([start, end])
    plt.ylim([0, SIZE_SCALE])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    xticks_start = (start//SIZE_SCALE + 1) * SIZE_SCALE
    xticks_list = [i for i in range(xticks_start, end, SIZE_SCALE)]
    ax.set_xticks(xticks_list)
    tmp_list = [i//1000 for i in xticks_list]
    ax.set_xticklabels([str(i//1000)+','+str(i%1000) for i in tmp_list], fontsize=10)




    def plot_tracks(chromStart, chromEnd, blockSizes, blockShift, thickStart, thickEnd, name,SIZE_SCALE, color, mode):
        mid_line = SIZE_SCALE//4 if mode == "ref" else SIZE_SCALE*3//4
        plt.plot([chromStart, chromEnd], [mid_line, mid_line], linestyle='-', color='black', alpha=0.7, lw=1, zorder=1)
        mid_point = (chromStart + chromEnd) // 2
        ARROW_SIZE = SIZE_SCALE//40
        plt.plot([mid_point-ARROW_SIZE, mid_point+ARROW_SIZE], [mid_line-ARROW_SIZE*3, mid_line], linestyle='-', color='black', alpha=0.8, lw=0.8)
        plt.plot([mid_point-ARROW_SIZE, mid_point+ARROW_SIZE], [mid_line+ARROW_SIZE*3, mid_line], linestyle='-', color='black', alpha=0.8, lw=0.8)
        plt.text(chromEnd+ARROW_SIZE, mid_line-ARROW_SIZE, name, fontsize=10 , color='black')

        def add_pathes(thickStart, thickEnd, chromStart, chromEnd, blockSizes, blockShift, SIZE_SCALE):
            
            for start, size in zip(blockShift, blockSizes):
                if start + size <= thickStart-chromStart or start >= thickEnd-chromStart:
                    rect = matplotlib.patches.Rectangle((chromStart + start, mid_line-SIZE_SCALE//20), size, SIZE_SCALE//10, color=color, alpha=1)
                    ax.add_patch(rect)
                elif start < thickStart-chromStart and start + size > thickStart-chromStart:
                    truncated_blockSize = (thickStart-chromStart) - start
                    rect = matplotlib.patches.Rectangle((chromStart + start, mid_line-SIZE_SCALE//20), truncated_blockSize, SIZE_SCALE//10, color=color, alpha=1)
                    ax.add_patch(rect)
                    new_start = chromStart + start + truncated_blockSize
                    new_size = size - truncated_blockSize
                    rect = matplotlib.patches.Rectangle((new_start, mid_line-SIZE_SCALE*3//20), new_size, SIZE_SCALE*3//10, color=color, alpha=1)
                    ax.add_patch(rect)
                elif start < thickEnd-chromStart and start + size > thickEnd-chromStart:
                    truncated_blockSize = (start + size) - (thickEnd-chromStart)
                    rect = matplotlib.patches.Rectangle((chromStart + start, mid_line-SIZE_SCALE*3//20), size - truncated_blockSize, SIZE_SCALE*3//10, color=color, alpha=1)
                    ax.add_patch(rect)
                    new_start = chromStart + start + size - truncated_blockSize
                    rect = matplotlib.patches.Rectangle((new_start, mid_line-SIZE_SCALE//20), truncated_blockSize, SIZE_SCALE//10, color=color, alpha=1)
                    ax.add_patch(rect)
                else:
                    rect = matplotlib.patches.Rectangle((chromStart + start, mid_line-SIZE_SCALE*3//20), size, SIZE_SCALE*3//10, color=color, alpha=1)
                    ax.add_patch(rect)

        add_pathes(thickStart, thickEnd, chromStart, chromEnd, blockSizes, blockShift, SIZE_SCALE)



    def plot_from_bed(file, color, mode,SIZE_SCALE):
        idx = 1
        reference_color = {1:"#482ff7", 2:"#086972"}
        old_start, old_end = 0, 0
        with open(file) as f:
            for q in f:
                q = q.rstrip().split('\t')
                chromStart, chromEnd = int(q[1]), int(q[2])
                name = q[3]
                blockSizes = [int(i) for i in q[10].split(',')[:-1]]
                blockShift = [int(i) for i in q[11].split(',')[:-1]]
                thichStart, thickEnd = int(q[6]), int(q[7])
                if mode == 'query':
                    plot_tracks(chromStart, chromEnd, blockSizes, blockShift, thichStart, thickEnd, name, SIZE_SCALE, color=color, mode=mode)
                else:
                    if chromStart < old_end and old_end < chromEnd:
                        idx += 1
                    elif old_start > chromStart and old_start < chromEnd:
                        idx += 1
                    plot_tracks(chromStart, chromEnd, blockSizes, blockShift, thichStart, thickEnd, name, SIZE_SCALE, color=reference_color[idx], mode=mode)
                    old_start, old_end = chromStart, chromEnd

    plot_from_bed("./downloads/reference.bed", color="#482ff7", mode="ref", SIZE_SCALE=SIZE_SCALE)
    plot_from_bed("./downloads/query.bed", color="#ff304f", mode="query", SIZE_SCALE=SIZE_SCALE)

    ax.set_xlabel('Chromosome '+chr_id, fontsize=15, fontweight='bold')
    ax.xaxis.set_label_coords(0.5, -0.2)

    return fig

def get_fig(enst1, enst2,info, seq, idx_fusion, primaryAligning):
    # try:
    if True:

        outinfo_1 = grep_bed_or_fasta(enst1).decode('gbk')
        outinfo_2 = grep_bed_or_fasta(enst2).decode('gbk')
        with open("./downloads/reference.bed", 'w') as w:
            print(outinfo_1,end="", file=w)
            print(outinfo_2,end="", file=w)

        outinfo_1 = outinfo_1[:].split('\n')[0].split("\t")
        outinfo_2 = outinfo_2[:].split('\n')[0].split("\t")
        chr_id = outinfo_1[0]
        strand = outinfo_1[5]
        positions = [int(outinfo_1[1]), int(outinfo_1[2]), int(outinfo_2[1]), int(outinfo_2[2])]
        start = min(positions) - 500
        end = max(positions) + 500


        with open("./downloads/tmp.bed", 'w') as f:
            print("\t".join(["chr"+chr_id, str(start), str(end), "region_tmp", "0", strand]), file=f)

        grep_bed_or_fasta(mode=False, strand=strand)

        with open("./downloads/seq2.fa", 'w') as f:
            print('>'+info+'\n'+seq, end="",file=f)

        align_line = 6 if primaryAligning else 7
        proc = subprocess.Popen(
        f"./downloads/bin/blat -t=dna -q=dna ./downloads/seq1.fa ./downloads/seq2.fa ./downloads/blat_out > /dev/null ; head -{align_line} ./downloads/blat_out | tail -1", 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True)
        outinfo, errinfo = proc.communicate()

        align_res = outinfo.decode('gbk').split('\t')
        '''
        match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
        ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        939	4	0	0	25	48	15	7824	+	ENST00000378842.1,GALT|ENST00000555981.5,IL11RA	3067	792	1783	chr9:34646362-34658983	12621	2458	11225	33	29,6,16,22,6,43,16,19,38,15,27,11,17,46,69,12,15,19,46,17,16,53,55,18,5,22,59,27,36,19,31,24,89,	792,824,831,848,870,877,920,938,958,997,1013,1043,1055,1073,1119,1189,1201,1219,1239,1291,1308,1325,1379,1435,1455,1462,1489,1550,1578,1614,1634,1670,1694,	2458,2488,2494,2510,2635,2641,2685,2701,3048,3086,3101,3128,3139,3156,8855,8924,8937,9239,9258,10381,10399,10415,10468,10523,10541,10673,10698,10757,10937,11058,11077,11111,11136,

        '''
        blockSizes = [int(i) for i in align_res[-3].split(',')[:-1]] # 29,6,16,22,6,43,16,19,38,15,27,11,17,46,69,12,15,19,46,17,16,53,55,18,5,22,59,27,36,19,31,24,89,	
        Ref_pos = [int(i) for i in align_res[-1].split(',')[:-1]] # 2458,2488,2494,2510,2635,2641,2685,2701,3048,3086,3101,3128,3139,3156,8855,8924,8937,9239,9258,10381,10399,10415,10468,10523,10541,10673,10698,10757,10937,11058,11077,11111,11136,
        Fatsq_pos = [int(i) for i in align_res[-2].split(',')[:-1]] # 792,824,831,848,870,877,920,938,958,997,1013,1043,1055,1073,1119,1189,1201,1219,1239,1291,1308,1325,1379,1435,1455,1462,1489,1550,1578,1614,1634,1670,1694,
        
        Ref_start, Ref_end, Fastq_start, Fastq_end = Ref_pos[0], Ref_pos[-1], Fatsq_pos[0], Fatsq_pos[-1]
        last_exonSize = blockSizes[-1]
        chromStart = start + Ref_pos[0] - Fatsq_pos[0] # 34646362 + 2458 - 792
        exonStart  = start + Ref_pos[0] # chromStart + 792
        exonEnd = start + Ref_pos[-1] + blockSizes[-1] # 34646362 + 11136 + 89
        chromEnd = exonEnd + len(seq) - last_exonSize - Fastq_end # exonEnd + 3067 - 1783
        blockSizes = [Fastq_start] + blockSizes[:] + [len(seq) - Fastq_end - last_exonSize]
        Ref_pos = [0] + [i - Ref_start + Fastq_start for i in Ref_pos[:]]
        Ref_end = Ref_pos[-1]
        Ref_pos += [Ref_end + last_exonSize]
    

        blockSizes = ",".join(to_str(blockSizes[:]))+","
        Ref_pos = ",".join(to_str(Ref_pos[:]))+","
        query_bed = to_str([chr_id, chromStart, chromEnd, "readthrough"+str(idx_fusion), 0, strand, exonStart, exonEnd, 0, len(Fatsq_pos)+2, blockSizes, Ref_pos])

        with open("./downloads/query.bed", 'w') as w:
            print("\t".join(query_bed), end="", file=w)


        final_start = start-1000 if start<chromStart else chromStart - 1000
        final_end = end+1000 if end>chromEnd else chromEnd + 1000

        fig = Plot(final_start, final_end, chr_id)
        return fig




# input using pasteboard
user_input = '''>ENST00000378842.1,GALT|ENST00000555981.5,IL11RA
AGATTTTTCCAGCGGATCCCCCGGTGGCCTCATGTCGCGCAGTGGAACCGATCCTCAGCAACGCCAGCAGGCGTCAGAGGCGGACGCCGCAGCAGCAACCTTCCGGGCAAACGCCCACCAGCATATCCGCTACAACCCGCCTGCAGGATGAGTGGGTGCTGGTGTCAGCTCACCGCATGAAGCGGCCCTGGCAGGGTCAAGTGGAGCCCCAGCTTCTGAAGACAGTGCCACCGCCATGACCCTCTCAACCCTCTGTGTCCTGGGGCCATCCGAGCCAACGGAGAGGTGAATCCCCAGTACGATAGCACCTTCCTGTTTGACAACGACTTCCCCAGCTCTGCAGCCTGATGCCCCCAGTCCAGGACCCAGTGATCATCCGCCTTTTCCAAGCAAAGTCTGCTCGAGGAGTCTGTAAGGTCATTGTGCTTCCACCCCTGGTCGGATGTAACGCTGCCACTCATGTCGCGTCCGCTGAGATCCGGGCTGTTTGTTGATGCATGGGCCTCAGTCACAGAGGAGCTGGGTGCCCAGTACCCTTGGTGCAGATCTTTGAAAACAAAGGTGCCATGATGGGCTGTTCTAACCCCCACCCCCACTGGCCAGGTATGGGCCAGCAGTTTCCTGCCAGATATTGCCGCCAGCGTGAGGAGCGATCTCAGCAGGCCATAAGTGAGTCAGCATGGAGAGCCCCTGCTAATGGAGTACAGCCGCCAGGAGCTACTCAGGAAGGAACGTCTGCCGCCAACCAGTGAGCACTGGTTAGTACTGGTGTCCCCTTCTGGGGCAACATGCGCCCTACCAGACACTGCTGCTGCCCCGTCCGCGCATGTGGCGGCGGCTACCTGAGGCTGACCCCTGCTGAGCGTGATGATCTAGCCCTCCATCATGAAGAAGCTCTTGACCAAGTATGACAACCTCTTGAGACGTCCTTTCCCTCTACTCCATGGGCTGGCATGGCGGCTCCCACAGGATCAGAGGCTGGGGCCAACTGGAACCCATTGGCAGCTGCACGGCTCATTACTACCCTCCGCTCCTGCGCTTTTCTGCCACTGTCTCGGAAATTCATGGTTGGTCTACGAAATGCTTGCTCAGGCTCAGAGGGACCTCACCCCTGAGCAGATGAGCAGCAGCTGCTCAGGGCTGAGCAGGGTCCTGGTGGCCGTGGCTACAGCCCTGGTGTCTGCCTCCCTCCCCCTGCCCCAGGCCTGGGGCCCCCACCCAGGGGTCCAGTATGGGCAAGCCAGGCAGGTCCGTGAAGCTGTGTTGTCCTGGAGTGACTGCCGGGGACCCACAGTGTCCTGGTTTCGGATGGGGAGCCAAAGCTGGCTCCAGGGACCTGACTCTGGCTAGGGGCATGAACTGGTCCTGGCCCAGGCAGGACAGCACTGATGAGGGCACCTACATCTGCCAGACCCTGGATGGTGCACTTGGGGGGCACAGTGACCCTGCAGCTACGGGCTACCCCTCCAGCCCGCCCTGTTGTCCTCCTTGCCAAGCAGCCGACTATGAGAACTTCTCTTGCACTTGGAGTCCCAGCCAGATCAGCGGCGTTTACCCACCCGCTACCTCACCTCCTAACAGGAAGAAGACAGTCCTAGGAGCTGATAGCCAGAGGAGGAGTCCATCCACAGGGGCCCTGGCCATGCCCACAGGATCCCCTAGGGGGCTGCCCCGCTGTGTTGTCCACGGGGCTGGTTCTGGAGCCAGTACCGGATTAATGTGACTGAGGTGAACCCACTGGGTGCCAGCACACGCCTGCTGGATGTGAGCTTGCAGAGCATCTTGCGCCCTGACCCACACCCAGGGCCTGCCGGGTGAGAGTCCAGTACCAGGTTACCCCCGACGCCTGCGAGCCAGCTGGACATACCCTGCCTCCTGGCCCGTGCGAGCCCCCACTTCCTGCTCAAGTTCCGTTTGCAGTACCCGTCCGGCGCACCATCCAGCCTGGTCCACGGTGGAGCCAGCTGGACTAGGAGGAGGTGATCACAGATGCTGTGGCTGGGCTGCCCATGCTGTACGAGTCAGTGCCCGGGACTTTCTAGATGCTGGCACTGGAGCACCTGGAGCCCGGAGAGGCCGGGGAACTCCGAGCACTCGGGACCATACCAAAGGAGATACCAGCATGGGGCCTAGCTACACACGCCAGCCAGAGGTGGAAGCCTCAGGTGGACAGCCCTGCTCCTCCAAGGCCCTCCCTCCAACCACACCCTCGGCTACTTGATCACAGGGACTCTGTGGAGCAGGTAGCTGTGCTTGGCGTCTTTGGGAATCCTTTCTTTCCTGGGACTGGTGGCTGGGGCCCTGGCACTGCGGGCTCTGGTAAGTGACTGCCATTGGTCCCTCAGCCTCTGATCCTCACACATGCTCTGATGCCCATAGACCACATTCATCTCCACCCTTCATGACTGCCCGCTGAACCTGTCTGATTTGCTGGAACTACCTCCCCATACCTCCATCCCCCATGCCCCACTTGATTTTAACTGATTCCTCTCCTGACCCTTTACTAATAAACCCTTTGGCGGGGAGACTGAGATAACCCACATTGTTGGAGAGACAGCTGCCTTTCTATGCCCCAGGCTGAGGCTGAGACGGGGGGGAAGGATGGATCCCCAAAGCCTGGGTTCTTGGCCTCAGTGATTCCAGTGGACAGGCGTCCAGAGCTCCAAAACCTGTAGAGGACCCCAGGAGGGCTCGTTTCGGCAGATTCCACCTATAATTCTGTCTTGCTGGTGTGGATAGAAACCAGGCAGGACAGTAGTATCCCTATGGTTGGATCTCAGCTGGAAGTTCTGTTTGGAGCCCATTTCTTGTGAGACCCTGTATTTCAAATTTGCAGCTGAAAGGTGCTTGTAACTCTGATTTCACCCCAGAGTTGGAGTTCTGCTCAAGGAAACGTGTGTAATGTGTACTCTGTGTCCATGTGTGACCATGTGTCTGTGAGGCAGGGACATGTATTCTCTGCATGCATGTATGTAGGTGCCTGGGGAGTGTGTGTGGGTCCTTGGCTCTTGGCCTTTCCCCTTGCAGGGGTTGTTGCAGGTGTGAATAAAGAGAATAAGGAAGTTCT
>ENST00000397832.7,IQCJ|ENST00000445224.6,SCHIP1
CTCTTTGAGGCTTTCTGATGCCTAGAGTGAGCTCAACTCAAAGGTTTCCAGCCTCACACTCGCCTCACATTCCCCCACAGTCACATTGCGCTCTGTGATTCTGAGGAATACAGTGTGCCAGCATCCGATCCAGTCTCCTTTCACCTGCAGGTGTTCCAGAAACTTCAAAATGCGTCTGGAAGAACTGAAAAGATTGCAGAATCCTCTAGAACAAGTTAATGATGGAAAATATTCATTTGAAAACATTCAGCGAGCATGGCGAGAGTACCTGCAGCGGCAGGAGCCCCTGGGGAAGAGGAGCCCGTCCCCACCCTCTGTCTCCTCAGAGAAGCTGAGCAGCTCTGTCAGCATGAACACCTTCTCCGACAGCAGCACACCCGATTACCGAGAGGATGGGATGGATCTAGGCAGTGACGCCGGCAGCAGCAGCAGCAGCAGCCGCGCCAGTTCACAGTCCAACTCCACCAAAGTGACCCCTTGCTCCGAGTGCAAATCTTCATCGTCGCCGGGGGGCAGCCTGGACTTGGTGTCTGCCCTGGAGGACTATGAGGAGCCCTTCCCGGTCTACCAGAAGAAGGTGATTGATGAGTGGGCGCCGGAGGAGGACGGGGAGGAGGAGGAAGAGGAGGACGAGCGCGACCAGCGAGGGTACCGGGATGACCGCTCTCCGGCCCGGGAACCGGGGGACGTAAGCGCCAGGACCCGCAGCGGCGGCGGCGGGGGCAGGAGCGCCACCACCGCCATGCCGCCCCCGGTGCCCAACGGCAACCTCCACCAGCACGACCCCCAGGACCTCAGGCACAATGGCAACGTGGTGGTGGCTGGCCGGCCGAGCTGTTCCCGGGGCCCCCGCCGGGCGATCCAAAAGCCCCAGCCGGCTGGGGGCCGGCGCAGTGGCCGCGGCCCGGCGGCTGGGGGGCTCTGCCTTCAGCCCCCAGACGGCGGGACGTGCGTCCCCGAAGAGCCCCCGGTGCCACCTATGGATTGGGAGGCGCTGGAGAAGCATCTGGCCGGGCTGCAGTTCCGGGAGCAGGAGGTACGGAACCAGGGCCAGGCGAGGACCAACTCCACCTCCGCACAGAAAAATGAGAGAGAGTCTATCAGACAGAAGTTGGCACTTGGAAGCTTCTTTGATGATGGCCCAGGAATTTATACCAGCTGTAGCAAAAGTGGGAAGCCAAGCCTTTCCTCCCGACTGCAGAGTGGGATGAACTTGCAGATATGCTTTGTCAACGACAGTGGCAGTGATAAGGACAGTGATGCTGATGACAGTAAGACTGAAACCAGCTTGGACACCCCCTTGTCTCCCATGAGCAAACAGAGTTCTTCCTATTCTGATAGAGACACTACTGAAGAGGAGTCTGAATCCTTGGATGACATGGACTTCCTTACAAGGCAAAAGAAATTGCAAGCTGAAGCCAAAATGGCCCTTGCCATGGCCAAACCAATGGCCAAAATGCAAGTAGAAGTGGAGAAACAGAACAGGAAAAAGTCTCCCGTCGCTGATCTTCTGCCACACATGCCTCATATAAGTGAATGCTTGATGAAAAGAAGTTTAAAACCCACCGACCTGAGAGACATGACTATTGGGCAGCTACAAGTGATAGTCAATGATCTCCATTCCCAGATAGAAAGCTTGAATGAAGAGTTGGTCCAGCTGCTTCTCATCCGAGATGAGCTGCACACAGAGCAGGATGCCATGCTGGTGGACATTGAAGACTTGACCAGACATGCTGAAAGTCAGCAGAAGCACATGGCAGAGAAAATGCCTGCAAAGTGA
'''
# user_input = '''>tmpname,IQCJ,SCHIP1
# CTCTTTGAGGCTTTCTGATGCCTAGAGTGAGCTCAACTCAAAGGTTTCCAGCCTCACACTCGCCTCACATTCCCCCACAGTCACATTGCGCTCTGTGATTCTGAGGAATACAGTGTGCCAGCATCCGATCCAGTCTCCTTTCACCTGCAGGTGTTCCAGAAACTTCAAAATGCGTCTGGAAGAACTGAAAAGATTGCAGAATCCTCTAGAACAAGTTAATGATGGAAAATATTCATTTGAAAACATTCAGCGAGCATGGCGAGAGTACCTGCAGCGGCAGGAGCCCCTGGGGAAGAGGAGCCCGTCCCCACCCTCTGTCTCCTCAGAGAAGCTGAGCAGCTCTGTCAGCATGAACACCTTCTCCGACAGCAGCACACCCGATTACCGAGAGGATGGGATGGATCTAGGCAGTGACGCCGGCAGCAGCAGCAGCAGCAGCCGCGCCAGTTCACAGTCCAACTCCACCAAAGTGACCCCTTGCTCCGAGTGCAAATCTTCATCGTCGCCGGGGGGCAGCCTGGACTTGGTGTCTGCCCTGGAGGACTATGAGGAGCCCTTCCCGGTCTACCAGAAGAAGGTGATTGATGAGTGGGCGCCGGAGGAGGACGGGGAGGAGGAGGAAGAGGAGGACGAGCGCGACCAGCGAGGGTACCGGGATGACCGCTCTCCGGCCCGGGAACCGGGGGACGTAAGCGCCAGGACCCGCAGCGGCGGCGGCGGGGGCAGGAGCGCCACCACCGCCATGCCGCCCCCGGTGCCCAACGGCAACCTCCACCAGCACGACCCCCAGGACCTCAGGCACAATGGCAACGTGGTGGTGGCTGGCCGGCCGAGCTGTTCCCGGGGCCCCCGCCGGGCGATCCAAAAGCCCCAGCCGGCTGGGGGCCGGCGCAGTGGCCGCGGCCCGGCGGCTGGGGGGCTCTGCCTTCAGCCCCCAGACGGCGGGACGTGCGTCCCCGAAGAGCCCCCGGTGCCACCTATGGATTGGGAGGCGCTGGAGAAGCATCTGGCCGGGCTGCAGTTCCGGGAGCAGGAGGTACGGAACCAGGGCCAGGCGAGGACCAACTCCACCTCCGCACAGAAAAATGAGAGAGAGTCTATCAGACAGAAGTTGGCACTTGGAAGCTTCTTTGATGATGGCCCAGGAATTTATACCAGCTGTAGCAAAAGTGGGAAGCCAAGCCTTTCCTCCCGACTGCAGAGTGGGATGAACTTGCAGATATGCTTTGTCAACGACAGTGGCAGTGATAAGGACAGTGATGCTGATGACAGTAAGACTGAAACCAGCTTGGACACCCCCTTGTCTCCCATGAGCAAACAGAGTTCTTCCTATTCTGATAGAGACACTACTGAAGAGGAGTCTGAATCCTTGGATGACATGGACTTCCTTACAAGGCAAAAGAAATTGCAAGCTGAAGCCAAAATGGCCCTTGCCATGGCCAAACCAATGGCCAAAATGCAAGTAGAAGTGGAGAAACAGAACAGGAAAAAGTCTCCCGTCGCTGATCTTCTGCCACACATGCCTCATATAAGTGAATGCTTGATGAAAAGAAGTTTAAAACCCACCGACCTGAGAGACATGACTATTGGGCAGCTACAAGTGATAGTCAATGATCTCCATTCCCAGATAGAAAGCTTGAATGAAGAGTTGGTCCAGCTGCTTCTCATCCGAGATGAGCTGCACACAGAGCAGGATGCCATGCTGGTGGACATTGAAGACTTGACCAGACATGCTGAAAGTCAGCAGAAGCACATGGCAGAGAAAATGCCTGCAAAGTGA
# '''

user_input = user_input.split('>')

# or load from file
user_input = []
idx = -1
with open("./downloads/example.fa") as f:
    for q in f:
        if '>' == q[0]:
            idx += 1
            user_input.append(q.rstrip())
        else:
            user_input[idx] += ('\n' + q.rstrip())
        


# PLOT

primaryAligning = True

idx = 0
error_html = ""
for record in user_input:
    if record != '':
        info, seq = "", ""
        idx += 1
        i = record.rstrip().split('\n')
        for k in i:
            if ',' in k :
                info += k
            else:
                seq += k

        if 'ENST' in info:
            enst1, enst2 = info.split('|')[0].split(',')[0].split('.')[0], info.split('|')[1].split(',')[0].split('.')[0]
        else:
            geneName1, geneName2 = info.split(',')[1], info.split(',')[2]
            outinfo = getEnst(geneName1).decode('gbk')
            matches = re.findall(r'transcript_id "(.*?)"', outinfo)
            enst1 = matches[0]
            outinfo = getEnst(geneName2).decode('gbk')
            matches = re.findall(r'transcript_id "(.*?)"', outinfo)
            enst2 = matches[0]
        fig = get_fig(enst1, enst2, info, seq, idx, primaryAligning)

        '''
        select secondary alignment
        '''
        primaryAligning = False 

        fig.show()
        

