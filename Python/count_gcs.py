#!/software/python3/Python-v3.7.0/bin/python3.7

import argparse
import re

#统计有可能是假阳性的突变位点G>C C>G
def main():
        parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-s','--sample',help='sample name',required=True)
#       parser.add_argument('-o','--outfile',help='outfile_name',required=True)
        parser.add_argument('-md','--mainDir',help='maindir',required=True)
        args=parser.parse_args()
        global headep
        global site_dict

        vcf_filename = "{0}/Analysis/{1}/{1}/Somatic/MergeAnno/{1}.MergeAnno.xls".format(args.mainDir,args.sample)
        with open(vcf_filename,'r') as vcf:
#               int ca,cg,ct=0
#               int ga,gt,gc=0
#               int ta,tg,tc=0
#               int ag,at,ac=0
                total=0
                cel=0
                for line in vcf:
                        if line.startswith("Chr"):
                                headep = re.split(pattern="\t",string=line)
                        else:
                                site_info = re.split(pattern="\t",string=line)
                                site_dict = dict(zip(headep,site_info))
                        if site_dict["Ref"]=="G" and site_dict["Alt"]=="C":
                                cel=cel+1
                        if site_dict.get["Ref"]=="C" and site_dict["Alt"]=="G":
                                cel=cel+1
                        total=total+1
                percentage = (cel/total)*100
                percentage = "{0}%".format(str(percentage))
                with open("out.xls",'a') as outer:
                        outer.write("{0}\t{1}\t{2}\t{3}\n".format(args.sample,str(total),percentage,str(cel)))
                outer.close()
        vcf.close()


if __name__ == '__main__':
        headep = {}
        site_dict = {}
        main()
