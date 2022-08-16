import re
import numpy as np

fasta = open("/data/Database/hg19/ucsc.hg19.fasta", 'r')
chr_this = ""

chr_comp = re.compile(pattern=">chr")
gcgc_comp = re.compile(pattern="GCGC")
mapping_short_seq = "GCGC"
line_this = ''
chr_pos = 1
out_bed = "out.bed"

bed = open(out_bed, "w")


# 转换为数值
def to_numeric(_string):
    trans_dict = {"A": 1, "T": -1, "G": 2, "C": -2, "a": 1, "t": -1, "c": 2, "g": -2, "N": 256, "n": 256}
    _string = _string.replace("\n", '')
    _string = list(_string)
    for _i in np.arange(0, len(_string)):
        _string[_i] = trans_dict[_string[_i]]
    return _string


def mapping_pos(short_seq, long_seq, cut_off):
    out_pos = []
    short_seq = to_numeric(short_seq)
    long_seq = to_numeric(long_seq)
    for _long_pos in np.arange(start_pos, start_pos + len(long_seq) - len(short_seq)):
        mapping_score = 0
        for _short_pos in np.arange(0, len(short_seq)):
            mapping_score += abs(long_seq[_long_pos + _short_pos - start_pos] - short_seq[_short_pos])
            if mapping_score > cut_off:
                break
        if mapping_score<= cut_off:
            out_pos.append(_long_pos)
    return out_pos


# 逐行读取，节约内存
# 双行拼接，防止接头处出现匹配序列


line_number = 1

while True:
    line_pre = line_this
    line_this = fasta.readline()
    chr_pos+=len(line_pre)
    #print("length line pre:{0}\nLength line this :{1}\nLine this={2}\nLine pre={3}".
    #format(len(line_pre),len(line_this),line_this,line_pre))
    #input()
    line_number += 1
    if not line_this:
        break
    elif chr_comp.match(line_this):
        chr_this = line_this.replace('>', "")
        chr_this = chr_this.replace("\n", "")
        chr_pos = 0
        line_this = ""
        print("New chr")
        continue
    elif not line_pre:
        continue
    else:
        line_running = line_pre+line_this
        line_running = line_running.replace("\n", "")
        pos_list = []
        pos_list = mapping_pos(short_seq=mapping_short_seq,long_seq=line_running, cut_off=0)
        if not pos_list:
            #print("poslist empty")
            #print(line_running)
            chr_pos += len(line_pre)
            continue
        else:
            for i in np.arange(0, len(pos_list)):
                start_pos = chr_pos + pos_list[i]-76
                end_pos = chr_pos + len(mapping_short_seq) + pos_list[i]-76
                #print("{0}\t{1}\t{2}\n".format(chr_this, start_pos, end_pos), file=bed, end="")
                print("{0}\t{1}\t{2}\n".format(chr_this, start_pos, end_pos), end="\t")
                print("Line Number:{0}\nChr name:{1}\nStart pos:{2}\nEnd pos:{3}\n".format(line_number,chr_this,start_pos,end_pos))
                #print(line_running[pos_list[i]:pos_list[i]+4])
                input()
            chr_pos += len(line_pre)
    if not chr_pos % 10:
        print(chr_pos)


fasta.close()
bed.close()
