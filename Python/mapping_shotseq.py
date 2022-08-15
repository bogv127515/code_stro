import numba
import numpy as np
import argparse
from Bio.SeqIO import parse
import sys
from numba import jit
from datetime import datetime
import os

'''
Authorized date 25/07/2022

程序功能：
检索用户给出的基因组序列文件中的一段短序列，设置特定阈值来对位置进行匹配，输出bed文件

参数列表：
-f --fasta_file ：待检索的fasta序列文件
-s --short_sequence：在基因组中待匹配的短序列 
-o --out_prefix：输出bed文件头的名称
-m --max_mismatch：匹配时允许最大错配碱基数量，默认为0

'''

__author__ = "liujiaxing"
__mail__ = "jiaxingliu@chosenmed.com"


# 通过窗口滑动得到最小错配的方法得到待匹配序列位置信息
@numba.njit()
def mapping_pos_mismatch(reads_seq, ref_seq, max_mismatch):
    out_pos = []
    for _ref_pos in np.arange(0, len(ref_seq) - len(reads_seq)):
        mismatch_counts = 0
        for _short_pos in range(0, len(reads_seq)):
            if ref_seq[_ref_pos + _short_pos] != reads_seq[_short_pos]:
                mismatch_counts += 1
                if mismatch_counts > max_mismatch:
                    break
                else:
                    continue
        if mismatch_counts <= max_mismatch:
            out_pos.append(_ref_pos)
    return out_pos

#判断参数是否正确，如果不正确raise一个错误FLAG
def judge(fasta,short_seq,out_pre,mm):
    date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ERROR_FLAG = False
    if not os.path.exists(fasta):
        print("{t}\ttest,FASTA file {fa} seems not exist or unreadable\nPlease check it!".
        format(t=date_now,fa=fasta))
        ERROR_FLAG = True
    elif os.path.exists(out_pre):
        print("{d}\ttest ,out put file {out} seems exist\nTry remove or rename?".
        format(d=date_now,out=out_pre))
        ERROR_FLAG = True
    else:
        for i in short_seq:
            if i not in ["A","G","T","C"]:
                ERROR_FLAG = True
    return ERROR_FLAG


def main():
    #参数部分
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='Author:\t{0}\nMail:\t{1}'.format(__author__, __mail__))
    parser.add_argument("-f", "--fasta_file", help="Input fasta file name", required=True)
    parser.add_argument("-s", "--short_sequence", type=str,help="Short sequences which will mapped to fasta", required=True)
    parser.add_argument("-o", "--out_file",
                        help="Out bed file name prefix,example: -o examp get out put file\n\texamp.bed", required=True)
    parser.add_argument("-m", "--max_mismatch", default=0, help="Max mismatch number[defult=0]", required=False)
    args = parser.parse_args()
    
    #判断参数是否正确
    if judge(fasta=args.fasta_file,short_seq=args.short_sequence,
    out_pre=args.out_file + ".bed",mm=args.max_mismatch):
        sys.exit()
    
    fasta = open(args.fasta_file, "r")
    bed_name = args.out_file + ".bed"
    bed = open(bed_name, 'w')
    shot_seq = args.short_sequence

    #使用BioPython读取整个染色体序列信息，生成迭代器，放到内存里
    for record in parse(fasta, "fasta"):
        #提取染色体序列信息，变成字符串类型（字符串几乎没有长度限制）
        seq = str(record.seq)
        seq = seq.upper()
        #匹配函数返回匹配位置列表
        pos_list = mapping_pos_mismatch(ref_seq=seq, reads_seq=shot_seq, max_mismatch=args.max_mismatch)
        #获取染色体名称
        chrom = record.id
        #实时汇报进度
        print(chrom)
        #写入bed文件
        for i in range(0, len(pos_list)):
            print("{0}\t{1}\t{2}\n".format(chrom, pos_list[i]+1, pos_list[i] + len(shot_seq)+1),
                  file=bed, end="")
    #关闭读写文件流
    fasta.close()
    bed.close()


if __name__ == "__main__":
    main()

