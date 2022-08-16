import os
import re
import sys

''' lib path '''
B_DIR   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(B_DIR,'conf'))
import settings
sys.path.append(settings.LIB_DIR)
from deal_ref_alt import deal_ref_alt
from amino_acid_trans import three2clinc_c,three2clinc_p,three2one

vcf2annovar_genefunction = {"synonymous_variant": "exonic", "missense_variant": "exonic", "inframe_insertion": "exonic",
                            "inframe_deletion": "exonic", "stop_gained": "exonic","frameshift_variant": "exonic",
                            "coding_sequence_variant": "exonic",
                            "upstream_gene_vairant":"upstream","downstream_gene_variant":"downstream",
                            "intron_variant":"intronic","splice_region_variant":"intronic","intergenic_variant":"intergenic",
                            "3_prime_UTR_variant":"UTR3","5_prime_UTR_variant":"UTR5",
                            "splice_donor_variant":"splicing","splice_acceptor_variant":"splicing",
                            "non_coding_transcript_exon_variant":"ncRNA_exonic","non_coding_transcript_intron_variant":"ncRNA_intronic"}




def vep_deal(f_in):
	with open(f_in,'r') as f:
		pos_dic = {}
		for line in f:
			if line.startswith('#') : continue
			# Uploaded_variation Location Allele Gene Feature Feature_type 
			# Consequence cDNA_position CDS_position Protein_position Amino_acids  
			# Codons  Existing_variation Extra 
			dat = line.strip().split('\t')
		
			#part1 get chr start end ref alt

			#chr1_9770691_AGAG/-    	chr1:9770691-9770694
			#chr10_63852322_TT/AG   	chr10:63852322-63852323 AG
			#chr1_9784501_-/ACCGAGGGG   chr1:9784500-9784501    ACCGAGGGG
			result = re.split(r'[_/]',dat[0])
			if len(result) > 4 : 
				print (result)
			(ch,pos,ref,alt) = re.split(r'[_/]',dat[0])[0:4]			
			if ref== '-' : pos = str(int(pos)-1)					
			(res0,res1) = deal_ref_alt(ch,pos,ref,alt)
			res1	= res1.replace('\t','_')
			(ch,start,end,ref,alt) = res1.split('_')		
	
			#part2 get aacchange
			exon_func 	= dat[6]
			info 		= re.split(r'[;=]',dat[-1])
			info_dic 	= {}
			for n in range(0,len(info),2):
				info_dic[info[n]] = info[n+1]

			if 'SYMBOL' in info_dic:
				gene = info_dic['SYMBOL']
			else:
				gene = 'Null'
			HGVS_OFFSET	= ''
			if 'HGVS_OFFSET' in info_dic:
				new_id = ch + ':' + str(int(start) + int(info_dic['HGVS_OFFSET']))
				HGVS_OFFSET = 'HGVS_OFFSET=' + info_dic['HGVS_OFFSET'] + ';' + new_id	

			if 'HGVSc' not in info_dic : 
				(hgvs_one,hgvs_three,hgvs_three_new) = ['.','.','.']
			else:
				pro_three	= info_dic['HGVSp'].split(':')[1] if 'HGVSp' in info_dic else ''
				pro_three_new	= three2clinc_p(pro_three)
				if '%3D' in pro_three_new:
					if exon_func == 'synonymous_variant' : 
						pro_three_new = pro_three_new.replace('%3D','=')
					else: pro_three_new = ''

				pro_one = three2one(pro_three_new)
				(nm,c)	= info_dic['HGVSc'].split(':')

				if 'EXON' in info_dic :
					number = 'exon' + info_dic['EXON'].split('/')[0]
					num_total = info_dic['EXON']
				elif 'INTRON' in  info_dic:
					number = 'intron' + info_dic['INTRON'].split('/')[0]
					num_total = info_dic['INTRON']
				else:
					number = ''	
					num_total = ''

				code 		= nm + ':' + number + ':' + c
				hgvs_three	= (code + ':' + pro_three).strip(':') + ':' + num_total
				hgvs_three_new 	= (code + ':' + pro_three_new).strip(':') + ':' + num_total
				hgvs_one	= (code + ':' + pro_one).strip(':') + ':' + num_total

				hgvs_three 	= hgvs_three.replace('::',':')
				hgvs_three_new 	= hgvs_three_new.replace('::',':')
				hgvs_one	= hgvs_one.replace('::',':')

			pos_dic[res1] = [HGVS_OFFSET,hgvs_one,hgvs_three,hgvs_three_new,gene]
		return(pos_dic)


def vep_deal_test(f_in):
	with open(f_in,'r') as f:
		pos_dic = {}
		for line in f:
			if line.startswith('#') : continue
			# Uploaded_variation Location Allele Gene Feature Feature_type 
			# Consequence cDNA_position CDS_position Protein_position Amino_acids  
			# Codons  Existing_variation Extra 
			dat = line.strip().split('\t')
		
			#part1 get chr start end ref alt

			#chr1_9770691_AGAG/-    	chr1:9770691-9770694
			#chr10_63852322_TT/AG   	chr10:63852322-63852323 AG
			#chr1_9784501_-/ACCGAGGGG   chr1:9784500-9784501    ACCGAGGGG
			result = re.split(r'[_/]',dat[0])
			if len(result) > 4 : 
				print (result)
			(ch,pos,ref,alt) = re.split(r'[_/]',dat[0])[0:4]			
			if ref== '-' : pos = str(int(pos)-1)					
			(res0,res1) = deal_ref_alt(ch,pos,ref,alt)
			res1	= res1.replace('\t','_')
			(ch,start,end,ref,alt) = res1.split('_')		
	
			#part2 get aacchange
			variaint_func 	= dat[6]
			info 		= re.split(r'[;=]',dat[-1])
			info_dic 	= {}
			for n in range(0,len(info),2):
				info_dic[info[n]] = info[n+1]

			if 'SYMBOL' in info_dic:
				gene = info_dic['SYMBOL']
			else:
				gene = 'Null'
			HGVS_OFFSET	= ''
			if 'HGVS_OFFSET' in info_dic:
				new_id = ch + ':' + str(int(start) + int(info_dic['HGVS_OFFSET']))
				HGVS_OFFSET = 'HGVS_OFFSET=' + info_dic['HGVS_OFFSET'] + ';' + new_id	

			if 'HGVSc' not in info_dic : 
				(hgvs_one,hgvs_three,hgvs_three_new) = ['.','.','.']
			else:
				pro_three	= info_dic['HGVSp'].split(':')[1] if 'HGVSp' in info_dic else ''
				pro_three_new	= three2clinc_p(pro_three)
				if '%3D' in pro_three_new:
					if variaint_func == 'synonymous_variant' : 
						pro_three_new = pro_three_new.replace('%3D','=')
					else: pro_three_new = ''

				pro_one = three2one(pro_three_new)
				(nm,c)	= info_dic['HGVSc'].split(':')

				if 'EXON' in info_dic :
					number = 'exon' + info_dic['EXON'].split('/')[0]
					num_total = info_dic['EXON']
				elif 'INTRON' in  info_dic:
					number = 'intron' + info_dic['INTRON'].split('/')[0]
					num_total = info_dic['INTRON']
				else:
					number = ''	
					num_total = ''

				code 		= nm + ':' + number + ':' + c
				hgvs_three	= (code + ':' + pro_three).strip(':') + ':' + num_total
				hgvs_three_new 	= (code + ':' + pro_three_new).strip(':') + ':' + num_total
				hgvs_one	= (code + ':' + pro_one).strip(':') + ':' + num_total

				hgvs_three 	= hgvs_three.replace('::',':')
				hgvs_three_new 	= hgvs_three_new.replace('::',':')
				hgvs_one	= hgvs_one.replace('::',':')

			#part3 get gene_function for recall
			if variaint_func in vcf2annovar_genefunction.keys():
				variaint_func = vcf2annovar_genefunction[variaint_func]
			else:variaint_func = variaint_func

			pos_dic[res1] = [HGVS_OFFSET,hgvs_one,hgvs_three,hgvs_three_new,gene,variaint_func]
		return(pos_dic)

if __name__ == '__main__':
	f_vep 	= sys.argv[1]
	pos_dic = vep_deal(f_vep) 
	#for sub in pos_dic:
	#	print(sub,pos_dic[sub][1],pos_dic[sub][2],pos_dic[sub][3])	
