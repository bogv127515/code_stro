import os
import re
import sys

''' lib path '''
B_DIR   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(B_DIR,'conf'))
import settings
sys.path.append(settings.LIB_DIR)
from deal_ref_alt import deal_ref_alt
from amino_acid_trans import get_transcript,three2clinc_c,three2clinc_p,three2one

''' get nm '''
def deal_info(ann,trans_dic):
	aa		= []
	aa_new		= []
	trans_info 	= ann.split(';')[0].split(',') 
	can_flag 	= False
	func		= trans_info[0].split('|')[1]
	gene		= trans_info[0].split('|')[3]
#	if 'synonymous_variant' in info:
		
	for info in trans_info :
		info_list 	= info.split('|')
		sub_func	= info_list[1]
		if sub_func == 'intergenic_region':continue
		(sub_gene,nm,exon,code,pro) = (info_list[3],info_list[6],info_list[8],info_list[9],info_list[10])
		#c.6666G>A
		if sub_gene != gene : continue
		if sub_func != func : continue
		if code.startswith('n.') : continue  

		code_new	= three2clinc_c(code)
		pro_new		= three2clinc_p(pro)
		
		nm = nm.split(':')[-1]
		exon1 = 'exon' + exon.split('/')[0]  if '/' in exon else ''
		if sub_func == 'intron_variant' : exon1 = exon.replace('exon','intron')
		sub_aa		= ':'.join([gene,nm,exon1,code,pro,exon]).replace('::',':').strip(':')
		sub_aa_new = ':'.join([gene,nm,exon1,code_new,pro_new,exon]).replace('::',':').strip(':')
	
		if sub_aa not in aa : 
			aa.append(sub_aa)
			aa_new.append(sub_aa_new)
		if gene in trans_dic:
			for sub in trans_dic[gene]:
				if sub == nm:
					Canonical_three		= sub_aa
					Canonical_three_new	= sub_aa_new
					can_flag = True
				break
				
	if aa == [] : 
		aa 		= ['.']
		aa_new 	= ['.']
	if can_flag == False: 
		Canonical_three		= aa[0]
		Canonical_three_new 	= aa_new[0]
	Canonical_one		= three2one(Canonical_three_new)	
	return(func,','.join(aa),Canonical_one,Canonical_three,Canonical_three_new)


def snpEff_deal(vcf,trans_dic):
	###20201228 同义突变结果更改。
	with open(vcf,'r') as f:
		pos_dic = {}
		for line in f:
			if line.startswith('#'): continue
			##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
			dat	= line.strip().split('\t')
			(ch,pos,ref,alt,info) = (dat[0],dat[1],dat[3],dat[4],dat[7])
			
			ch_list = ch.split('_')
			if len(ch_list) > 1 :continue

			(res0,res1) = deal_ref_alt(ch,pos,ref,alt)
			res1 	= res1.replace('\t','_')
			print(res1)
			(ch1,start,end,ref1,alt1) = res1.split('_')
			m = re.match(r'.*ANN=(.*)',info)
			synonymous_snpEff_flag = False
			(aa,Canonical_one,Canonical_three) = ('.','.','.')
			if m :
				synonymous_snpEff_flag = True if 'synonymous_variant' in m.group(1) else False 
				##add at 20201228
				(func,aa,Canonical_one,Canonical_three,Canonical_three_new) = deal_info(m.group(1),trans_dic)
			if aa == '.' : (Canonical_one,Canonical_three) = ('.','.')
			pos_dic[res1] = [aa,Canonical_one,Canonical_three,synonymous_snpEff_flag,func]

	return(pos_dic)
	
def snpEff_deal_bk(vcf,trans_dic):
	###20201228 同义突变结果更改。
	with open(vcf,'r') as f:
		pos_dic = {}
		for line in f:
			if line.startswith('#'): continue
			##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
			dat	= line.strip().split('\t')
			(ch,pos,ref,alt,info) = (dat[0],dat[1],dat[3],dat[4],dat[7])
			
			ch_list = ch.split('_')
			if len(ch_list) > 1 :continue

			(res0,res1) = deal_ref_alt(ch,pos,ref,alt)
			res1 	= res1.replace('\t','_')
			print(res1)
			(ch1,start,end,ref1,alt1) = res1.split('_')
			m = re.match(r'.*ANN=(.*)',info)
			print(m.group(1))
			synonymous_snpEff_flag = False
			(aa,Canonical_one,Canonical_three) = ('.','.','.')
			if m :
				synonymous_snpEff_flag = True if 'synonymous_variant' in m.group(1) else False 
				##add at 20201228
				(func,aa,Canonical_one,Canonical_three,Canonical_three_new) = deal_info(m.group(1),trans_dic)
			if aa == '.' : (Canonical_one,Canonical_three) = ('.','.')
			pos_dic[res1] = [aa,Canonical_one,Canonical_three,synonymous_snpEff_flag]

	return(pos_dic)

if __name__ == '__main__':
	snpeff_vcf 	= sys.argv[1]
	f_trans 	= sys.argv[2] if len(sys.argv) > 2 else '/data/Pipeline/liaorui/database/Transcript/all_canonical_trans.txt' 

	trans_dic	= get_transcript(f_trans)
	pos_dic		= snpEff_deal(snpeff_vcf,trans_dic)

	for pos in pos_dic:
		print(pos,'\t','\t'.join(pos_dic[pos]))
