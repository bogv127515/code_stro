import argparse
import pandas as pd
import os


def trans_count(ref_seq,alt_seq,ref_list,alt_list,reverse=False):

    """
    Parameters:
        ref_seq - 基因组中应有的位点信息
        alt_seq - 肿瘤组织中突变的位点信息
        ref_list - 检测出全部突变位点在基因组中应有的序列
        alt_list - 全部的突变位点序列
        reverse - 是否进行反向 ref<->alt 匹配

    Returns:
        返回输入突变类型在输入队列中的频度信息

    Raises:
        KeyError - raises an exception 
    
    
    """
    line_count = len(ref_list)-1
    atten_count = 0
    while line_count!=0:
        if ref_list[line_count]==ref_seq and alt_list[line_count]==alt_seq:
            atten_count+=1
        if reverse:
            if ref_list[line_count]==alt_seq and alt_list[line_count]==ref_seq:
                atten_count+=1
        line_count-=1
    return atten_count



def main(sample):
    #找到VCF或是统计表格的路径
    #Python中可以使用*作为通配符匹配路径，操作会更为灵活
    vcf_name = "{0}/{1}/{1}/Somatic/MergeAnno/{1}.MergeAnno.xls".format(args.WorkingDir,sample)
    vcf_df = pd.read_csv(vcf_name,header=0,sep="\t")

    #提取参考列和突变列信息，方便统计类型
    ref = vcf_df["Ref"]
    alt = vcf_df["Alt"]

    #调用突变类型频度统计函数，统计位点类型
    gc_trans = trans_count(ref_seq="C",alt_seq="G",ref_list=ref,alt_list=alt,reverse=True)
    gt_plus_ca_trans = trans_count(ref_seq="G",alt_seq="T",ref_list=ref,alt_list=alt)+trans_count(ref_seq="C",alt_seq="A",ref_list=ref,alt_list=alt)
    
    #打开写出结果表格
    with open("out.xls",'a') as outer:
        #统计所有位点类型信息
        #此处不用加一，只有用到索引数字的时候才进行加减一操作
        total=len(vcf_df["Chr"])
        #统计G<->C突变事件的占比
        gc_percentage = (gc_trans/total)*100
        gc_percentage = str(round(gc_percentage,2))

        #统计G->T及C->A事件发生的百分比
        gt_plus_ca_percentage = (gt_plus_ca_trans/total)*100
        gt_plus_ca_percentage = str(round(gt_plus_ca_percentage,2))

        #按表格形式写出到结果表格中
        #如果有很多列的话可以尝试一下创建一个表头字典，用字典写入处理过后的每行或每个样本的数据
        outer.write("{0}\t{1}\t{2}%\t{3}%\n".format(sample,str(total),gc_percentage,gt_plus_ca_percentage))
    #正常情况下应该关闭文件输出流，但使用with open的时候直接结束with语句块就行，不用额外进行文件流的关闭

    #最后汇报交过完成情况
    print("Sample {0} done!".format(sample))


if __name__ =="__main__":

    #只使用样本列表文件以及工作目录参数，方便在任何位置（例如工作目录下）用其他方式调用该脚本
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s','--sampleList',help='sample name',required=True)
    parser.add_argument('-wd','--WorkingDir',help='Celcu Directory ',required=True)
    args=parser.parse_args()

    #输入表头，可以随时更改
    with open("out.xls","a") as outer:
        outer.write("样本名称\t突变位点数\tC>G + G>C 占比\tG>T + C>A 占比\n")
        outer.close()

    #读入样本表格
    sample_list = pd.read_csv(args.sampleList,header=0,sep="\t")
    #拆分样品信息输入主流程处理函数中逐一处理
    for sn in sample_list["order_number"]:
        main(sample=sn)
    print("All done. Cheers!")
