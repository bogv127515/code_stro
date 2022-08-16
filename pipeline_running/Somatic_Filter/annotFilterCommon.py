import numpy as np
#import matplotlib.pyplot as plt
import scipy.stats as stats
import os
import re
import sys
import json


####################################### lib path  #####################################################
path = '/data/pipe_panel/database/pylib'
sys.path.append(os.path.join(path))
B_DIR   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(B_DIR,'conf'))
import settings
sys.path.append(settings.LIB_DIR)
from deal_ref_alt import deal_ref_alt
from amino_acid_trans import get_transcript,one2clinc_c,one2clinc_p	


#################################### Annot Need Files #################################################
annot_files_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'annot_files')
f_trans 		= os.path.join(annot_files_dir,'hg19_canonical_trans.txt')
f_hgvs			= os.path.join(annot_files_dir,'hg19_variants_hgvs.txt')
f_cosmic_hotspot= os.path.join(annot_files_dir,'hg19_cosmic_hotspots.txt')
f_pathogenic	= os.path.join(annot_files_dir,'hg19_pathogenic_variant.txt')
f_black_27		= os.path.join(annot_files_dir,'hg19_panel27_black.txt')
f_black_599		= os.path.join(annot_files_dir,'hg19_black.txt')
f_correct_hgvs	= os.path.join(annot_files_dir,'f_correct_hgvs.txt')
f_germline_white_hgvs	= os.path.join(annot_files_dir,'hg19_Germline_Somatic_variants_hgvs.txt')
f_ddpcr			= os.path.join(annot_files_dir,'hg19_ddpcr_hgvs.txt')
f_OMIM			= os.path.join(annot_files_dir,'hg19_OMIM.list')
f_OncoKB		= os.path.join(annot_files_dir,'hg19_OncoKB.list')
f_Black_YCZ		= os.path.join(annot_files_dir,'hg19_Black_YCZ_20211009.list')
f_InPanel		= os.path.join(annot_files_dir,'NanOnco_Plus_Panel_v2_62.sort.bed')
f_InExon		= os.path.join(annot_files_dir,'somatic_mut.last.bed')
f_HGNC			= os.path.join(annot_files_dir,'HGNC.genesymbol.xls')
f_tmb			= os.path.join(annot_files_dir,'TMB_mut.last.bed')
f_false 		= {
				'27' 	: os.path.join(annot_files_dir,'27_false_site.txt'),
				'599' 	: os.path.join(annot_files_dir,'599_false_site.txt') }

trans_dic = get_transcript(f_trans)

out_all_title_pair = ['Sample','Chr','Start','End','Ref','Alt','CHR','Gene','Type','Transcript','cHGVS',
				   'pHGVS','VAF','Consequence','Affected_Exon','Transcript_ann','cHGVS_ann','pHGVS_ann',
				   'Affected_Exon_ann','Gene.refGene','Func.refGene',
				   'ExonicFunc.refGene','AAChange.refGene','cytoBand','avsnp150','snp138','DM',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG',
				'cosmic88_coding','cosmic_count','cosmic_hotspots',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'genome_AF','genome_AF_eas','exome_AF','exome_AF_eas',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_all','1000g2015aug_eas',
				'1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur',
				'Annovar_Trans','Vep_Trans','SnpEff_Trans','SnpEff_annot','Annot_Software',
				'Vep_Trans_Three','SnpEff_Trans_Three','HGVS_OFFSET','InterVar_automated',
				'Canonical.transcript','GT','Coverage','Ref_Reads','Alt_Reads','Var',
				'N_GT','N_Coverage','N_Ref_Reads','N_Alt_Reads','N_Var',
				'Filter','FILTER_flag','QC_flag','Gene_flag','LOD_flag','Func_flag','Exon_func_flag',
				'Synonymous_flag','DB_flag','MAF_flag','SNP_flag','Cosmic_flag','White_flag',
				'Black_flag','Annot_flag','Count_27','Count_599','PASS_Count','Adjust_PASS',
				'Type','Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','Black_YCZ',
				'InPanel','InExon','InTMB','HGNC','Trans_Flag','Trans_Start',
				'Trans_Ref','Trans_Alt','Trans_Input','Strain','Trans_C','Trans_P','Trans_Function'
]


title_mrd_all = ['Sample','Chr','Start','End','Ref','Alt','Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','avsnp150','snp138','DM','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding','cosmic_count','cosmic_hotspots','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','genome_AF','genome_AF_eas','exome_AF','exome_AF_eas','esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_all','1000g2015aug_eas','1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','Annovar_Trans','Vep_Trans','SnpEff_Trans','SnpEff_annot','Annot_Software','Vep_Trans_Three','SnpEff_Trans_Three','HGVS_OFFSET','InterVar_automated','AAChange.1','GT','Coverage','Ref_Reads','Alt_Reads','Var','N_GT','N_Coverage','N_Ref_Reads','N_Alt_Reads','N_Var','Filter','FILTER_flag','QC_flag','LOD_flag','Func_flag','Exon_func_flag','Synonymous_flag','DB_flag','MAF_flag','SNP_flag','Cosmic_flag','White_flag','Black_flag','Normal_flag','MAF_RD_flag','SBias_flag','DatabaseSave_flag','Annot_flag','PASS_Count','Adjust_PASS','Type','Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','N_Var','Black_YCZ','MRD_PASS_flag']

#out_all_title_pair = ['Sample','Chr','Start','End','Ref','Alt','CHR','Gene','Type','Transcript','cHGVS',
#				   'pHGVS','VAF','Consequence','Affected_Exon','Transcript_ann','cHGVS_ann','pHGVS_ann',
#				   'Affected_Exon_ann','Gene.refGene','Func.refGene',
#				   'ExonicFunc.refGene','AAChange.refGene','cytoBand','avsnp150','snp138','DM',
#				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG',
#				'cosmic88_coding','cosmic_count','cosmic_hotspots',
#				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
#				'genome_AF','genome_AF_eas','exome_AF','exome_AF_eas',
#				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_all','1000g2015aug_eas',
#				'1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur',
#				'Annovar_Trans','Vep_Trans','SnpEff_Trans','SnpEff_annot','Annot_Software',
#				'Vep_Trans_Three','SnpEff_Trans_Three','HGVS_OFFSET','InterVar_automated',
#				'Canonical.transcript','GT','Coverage','Ref_Reads','Alt_Reads','Var',
#				'Filter','FILTER_flag','QC_flag','Gene_flag','LOD_flag','Func_flag','Exon_func_flag',
#				'Synonymous_flag','DB_flag','MAF_flag','SNP_flag','Cosmic_flag','White_flag',
#				'Black_flag','Annot_flag','Count_27','Count_599','PASS_Count','Adjust_PASS',
#				'Type','Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','Black_YCZ',
#				'InPanel','InExon','InTMB','HGNC','Trans_Flag','Trans_Start',
#				'Trans_Ref','Trans_Alt','Trans_Input','Strain','Trans_C','Trans_P','Trans_Function'
#]

out_all_title_single = ['Sample','Chr','Start','End','Ref','Alt','CHR','Gene','Type','Transcript','cHGVS',
				   'pHGVS','VAF','Consequence','Affected_Exon','Transcript_ann','cHGVS_ann','pHGVS_ann',
				   'Affected_Exon_ann','Gene.refGene','Func.refGene',
				   'ExonicFunc.refGene','AAChange.refGene','cytoBand','avsnp150','snp138','DM',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG',
				'cosmic88_coding','cosmic_count','cosmic_hotspots',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'genome_AF','genome_AF_eas','exome_AF','exome_AF_eas',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_all','1000g2015aug_eas',
				'1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur',
				'Annovar_Trans','Vep_Trans','SnpEff_Trans','SnpEff_annot','Annot_Software',
				'Vep_Trans_Three','SnpEff_Trans_Three','HGVS_OFFSET','InterVar_automated',
				'Canonical.transcript','GT','Coverage','Ref_Reads','Alt_Reads','Var',
				'Filter','FILTER_flag','QC_flag','Gene_flag','LOD_flag','Func_flag','Exon_func_flag',
				'Synonymous_flag','DB_flag','MAF_flag','SNP_flag','Cosmic_flag','White_flag',
				'Black_flag','Annot_flag','Count_27','Count_599','PASS_Count','Adjust_PASS',
				'Type','Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','Black_YCZ',
				'InPanel','InExon','InTMB','HGNC','Trans_Flag','Trans_Start',
				'Trans_Ref','Trans_Alt','Trans_Input','Strain','Trans_C','Trans_P','Trans_Function'
]

out_all_title_germline = ['Sample','CHR','Start','End','Ref','Alt','Chr','Gene','Type','Transcript','cHGVS','pHGVS','VAF','Consequence','Affected_Exon','Transcript_ann','cHGVS_ann','pHGVS_ann','Affected_Exon_ann','Gene.refGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','avsnp150','snp138','DM','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding','cosmic_count','cosmic_hotspots','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','genome_AF','genome_AF_eas','exome_AF','exome_AF_eas','esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_all','1000g2015aug_eas','1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','Annovar_Trans','Vep_Trans','SnpEff_Trans','SnpEff_annot','Annot_Software','Vep_Trans_Three','SnpEff_Trans_Three','HGVS_OFFSET','InterVar_automated','Canonical.transcript','GT','Coverage','Ref_Reads','Alt_Reads','Var','Filter','FILTER_flag','QC_flag','Gene_flag','Func_flag','Exon_func_flag','Synonymous_flag','DB_flag','MAF_flag','SNP_flag','Cosmic_flag','InPanel','InExon','HGNC','Trans_Flag','Trans_Start','Trans_Ref','Trans_Alt','Trans_Input','Strain','Trans_C','Trans_P','Trans_Function'
]

out_sub_title_Germline = out_all_title_germline

out_sub_title_Somatic_pair = ['Chr','Start','End','Ref','Alt',
				'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
				'cytoBand','1000g2015aug_all','avsnp150','snp138',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
				'1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
				'GT','AAChange.1','Ref_Reads','Alt_Reads','Var','Description_OncoKB','Ipen_YCZ','DP4',
				'FATHMM_pred','N_Var','Black_YCZ','QC_flag','LOD_flag','Synonymous_flag',
				'MAF_flag','Black_flag'
]

out_sub_title_Somatic_pair = ['Chr','Start','End','Ref','Alt',
                                'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
                                'cytoBand','1000g2015aug_all','avsnp150','snp138',
                                'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
                                'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                                'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
                                '1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
                                'GT','AAChange.1','Ref_Reads','Alt_Reads','Var','Description_OncoKB','Ipen_YCZ','DP4',
				'FATHMM_pred','N_Var','Black_YCZ'
]

out_sub_title_Somatic_pair_YCZ = ['Chr','Start','End','Ref','Alt',
                                'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
                                'cytoBand','1000g2015aug_all','avsnp150','snp138',
                                'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
                                'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                                'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
                                '1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
                                'GT','AAChange.1','Ref_Reads','Alt_Reads','Var','Description_OncoKB','Ipen_YCZ','DP4',
				'FATHMM_pred','N_Var','Black_YCZ','QC_flag','LOD_flag','Synonymous_flag',
				'MAF_flag','Black_flag'
]

out_sub_title_Somatic_single_YCZ = ['Chr','Start','End','Ref','Alt',
                                'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
                                'cytoBand','1000g2015aug_all','avsnp150','snp138',
                                'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
                                'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
                                'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
                                '1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
                                'GT','AAChange.1','Ref_Reads','Alt_Reads','Var','Description_OncoKB','Ipen_YCZ',
				'DP4','FATHMM_pred','Black_YCZ','QC_flag','LOD_flag','Synonymous_flag','MAF_flag','Black_flag'
]

out_sub_title_Somatic_single = ['Chr','Start','End','Ref','Alt',
				'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
				'cytoBand','1000g2015aug_all','avsnp150','snp138',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
				'1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
				'GT','AAChange.1','Ref_Reads','Alt_Reads','Var','Description_OncoKB','Ipen_YCZ',
				'DP4','FATHMM_pred','Black_YCZ'
]

out_sub_check_title_pair = ['Sample','IDD','ID','Check','Count_27','Count_599',
				'Filter','MuTect2','VarDict','VarScan','PASS_Count',
				'Type','Chr','Start','End','Ref','Alt','Bam',
				'Ref_Reads','Alt_Reads','Var','GT','AAChange.1',
				'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
				'cytoBand','1000g2015aug_all','avsnp150','snp138',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
				'1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
				'Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','N_Var','Black_YCZ'
]

out_sub_check_title_single = ['Sample','IDD','ID','Check','Count_27','Count_599',
				'Filter','MuTect2','VarDict','VarScan','PASS_Count',
				'Type','Chr','Start','End','Ref','Alt','Bam',
				'Ref_Reads','Alt_Reads','Var','GT','AAChange.1',
				'Func.refGene','Gene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene',
				'cytoBand','1000g2015aug_all','avsnp150','snp138',
				'CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG','cosmic88_coding',
				'SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred',
				'esp6500siv2_all','ExAC_ALL','ExAC_EAS','1000g2015aug_eas','1000g2015aug_sas',
				'1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur','InterVar_automated',
				'Description_OncoKB','Ipen_YCZ','DP4','FATHMM_pred','Black_YCZ'
]

germline_clinvar_list = ["Affects","Affects,association","Affects,risk_factor","association,risk_factor","Benign/Likely_benign,risk_factor","Benign,drug_response","Benign,risk_factor","Conflicting_interpretations_of_pathogenicity","Conflicting_interpretations_of_pathogenicity,Affects","Conflicting_interpretations_of_pathogenicity,drug_response","Conflicting_interpretations_of_pathogenicity,other,risk_factor","Conflicting_interpretations_of_pathogenicity,risk_factor","drug_response","drug_response,other","drug_response,protective,risk_factor","drug_response,risk_factor","Likely_benign,drug_response","Likely_pathogenic","Likely_pathogenic,association","Likely_pathogenic,drug_response","Likely_pathogenic,other","Likely_pathogenic,risk_factor","Pathogenic","Pathogenic/Likely_pathogenic","Pathogenic/Likely_pathogenic,drug_response","Pathogenic/Likely_pathogenic,other","Pathogenic/Likely_pathogenic,risk_factor","Pathogenic,Affects","Pathogenic,association","Pathogenic,association,protective","Pathogenic,drug_response","Pathogenic,other","Pathogenic,protective","Pathogenic,risk_factor","protective,risk_factor","risk_factor","Uncertain_significance,drug_response","Uncertain_significance,risk_factor","Uncertain_significance,association"]

### get canonical transcript and correct hgvs format ##################################################
''' get canonical trans from annovar,vep,snpeff '''
def get_canonical_trans_bk(dat,idd,trans_dic):

	canonical_trans = ''
	Canonical_soft	= 'Annovar'
	annot_flag	= True
	
	var_type 	= dat['Type']
	
	annovar_trans 	= dat['Annovar_Trans']
	vep_trans	= dat['Vep_Trans']
	snpeff_trans 	= dat['SnpEff_Trans']
	#offset		= dat['HGVS_OFFSET']
	
	annovar_ls	= annovar_trans.split(':')
	vep_ls		= vep_trans.split(':')
	snpeff_ls	= snpeff_trans.split(':')
	
	#并不一定选择的就是NM系列的
	annovar_nm	= annovar_ls[1] if len(annovar_ls) > 2 else ''  
	(vep_nm,vep_exon)	= vep_ls[1:3] if len(vep_ls) > 3 else ['','']   
	(snpeff_gene,snpeff_nm,snpeff_exon)	= snpeff_ls[0:3] if len(snpeff_ls) > 3 else ['','',''] 

	#判断哪些软件注释出经典转录本的flag，每个flag匹配经典转录本包含关系的唯一状态
	# nm_flag : 0 无匹配 
	# 1：annovar 3：annovar snpeff 5：annovar vep 7：annovar snpeff vep
	# 2：snpeff 6：snpeff vep
	# 4: vep
 
	'''
		Mod date:20220809
		1. 在annovar策略生效的基础上：
			● Flag为0、1、3、5、7，即无经典转录本和annovar包含经典转录本的情况，写入annovar结果
			● Flag为4、6，即含有Vep不含有annovar结果时，考虑Vep临床意义，写入Vep
			● Flag为2的时候，即只有SnpEff有经典转录本，写入SnpEff结果
		2. 在snpeff策略生效的基础上：
			● Flag为0、2、3、6、7，即无经典转录本和SnpEff包含经典转录本条件下，写入SnpEff结果
			● Flag为4、5，即无SnpEff结果同时包含Vep结果时，选择Vep结果
			● Flag为1，即只有Annovar包含经典转录本条件下，选择Annovar结果
		3. 在Vep策略生效的基础上：
			● Flag为0、4、5、6、7，即无经典转录本和Vep包含经典转录本条件下，写入Vep结果
			● Flag为1、3，即无Vep同时存在annovar结果时，优先选择annovar结果
			● Flag为2，即只有SnpEff包含经典转录本条件下，选择SnpEff结果

	'''


	dict_key = trans_dic.keys()
	nm_flag = 0
	if annovar_nm in dict_key: nm_flag = nm_flag + 1
	if snpeff_nm in dict_key: nm_flag = nm_flag + 2
	if vep_nm in dict_key: nm_flag = nm_flag + 4

	if var_type == 'snv' and annovar_trans != '.':
		print("annovar_trans start")
		if nm_flag in [0,1,3,5,7]:
			canonical_trans = annovar_trans           ##snv 注释为annovar的结果，更为准确。
			dat['Gene_c'] = canonical_trans.split(':')[0]
		elif nm_flag in [4,6]:
			canonical_trans = vep_trans
			Canonical_soft	= 'Vep' #转录本筛选方面，VEP考虑更多临床意义。
			dat['Gene_c'] = canonical_trans.split(':')[0]
		else: 
			canonical_trans = snpeff_trans
			Canonical_soft = 'SnpEff'
			dat['Gene_c'] = canonical_trans.split(':')[0]
	elif vep_nm == annovar_nm and not vep_exon.startswith('intron'): 
		print("vep start")
		##非snv的条目采用VEP进行注释，vep和annovar的转录本相同。
		if nm_flag in [0,4,5,6,7]:
			canonical_trans = vep_trans
			Canonical_soft	= 'Vep'
			dat['Gene_c'] = canonical_trans.split(':')[0]
		elif nm_flag in [1,3]:
			canonical_trans = annovar_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
		else:
			canonical_trans = snpeff_trans
			Canonical_soft = 'SnpEff'
			dat['Gene_c'] = canonical_trans.split(':')[0]
	elif snpeff_gene == dat['Gene.refGene'] and snpeff_nm.startswith('NM_') and 'exon' in snpeff_exon:
		print("snpeff start")
		##非snv的条目，且vep和annovar的转录本相同不相同。则采用SnpEff。
		if nm_flag in [0,2,3,6,7]:
			canonical_trans = snpeff_trans
			Canonical_soft = 'SnpEff'
			dat['Gene_c'] = canonical_trans.split(':')[0]
		elif nm_flag in [1,5]:
			canonical_trans = annovar_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
		else:
			canonical_trans = vep_trans
			Canonical_soft	= 'Vep' 
			dat['Gene_c'] = canonical_trans.split(':')[0]
	else:
		if nm_flag in [1,3,5,7]:
			#if dat['AAChange.refGene'] == 'UNKNOWN': 
			canonical_trans = annovar_trans 
			dat['Gene_c'] = canonical_trans.split(':')[0]
		elif nm_flag in [4,6,5]:
			canonical_trans = vep_trans
			Canonical_soft  = 'Vep'
			dat['Gene_c'] = canonical_trans.split(':')[0]
		else:
			canonical_trans = snpeff_trans
			Canonical_soft = 'SnpEff'
			dat['Gene_c'] = canonical_trans.split(':')[0]
		if '-AS1' in dat['Gene_c'] and annovar_trans == '.':
			canonical_trans = snpeff_trans
			Canonical_soft = 'SnpEff'
			dat['Gene_c'] = canonical_trans.split(':')[0]
		if '-AS1' in dat['Gene_c'] and annovar_trans != '.':
			canonical_trans = annovar_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
		#else:
		#	canonical_trans = annovar_trans
		#	dat['Gene_c'] = canonical_trans.split(':')[0]
		##其他情况采用Vep的结果。

	'''add at 20210427 Canonical selected from three software'''
	list_nm = [annovar_nm,vep_nm,snpeff_nm]
	l_trans = [annovar_trans,vep_trans,snpeff_trans]
	list_sf = ['Annovar','Vep','SnpEff']
	n = 0
	Canonical_soft_t = Canonical_soft

	#if dat['Start'] == '106158509': print(Canonical_soft,idd,canonical_trans,list_nm[1],list_nm[2],dat['Gene'])
	if dat['Gene_c'] != '.' and dat['Gene_c'] in dict_key:
		jd = trans_dic[dat['Gene_c']]
		for nm_t in list_nm:
			n = n + 1
			if Canonical_soft == 'Annovar' and nm_t == annovar_nm and list_nm[0] in jd:continue
			if Canonical_soft == 'Vep' and nm_t == vep_nm and list_nm[1] in jd : continue
			if Canonical_soft == 'SnpEff' and nm_t == snpeff_nm and list_nm[2] in jd: continue
			for mm in range(len(jd)):
				if nm_t == jd[mm]:
					if 'intron' in l_trans[n-1] :continue
					if 'p.' not in l_trans[n-1] :continue
					canonical_trans = l_trans[n-1]
					Canonical_soft_t= list_sf[n-1]
					dat['Gene_c'] = canonical_trans.split(':')[0]
				break
			if Canonical_soft != Canonical_soft_t: break

		Canonical_soft = Canonical_soft_t

#	if dat['Gene'] == 'TET2':print(Canonical_soft,idd,canonical_trans)
	if canonical_trans == '.' or annovar_trans == 'UNKNOWN' or 'UTR5' in annovar_trans : annot_flag = False
	dat['Gene.refGene'] = dat['Gene_c'] if  dat['Gene_c'] != '.' else dat['Gene.refGene']
	## 未被注释到。

	return(canonical_trans,annot_flag,Canonical_soft,dat['Gene.refGene'])

''' get_canonical_trans_test from annovar,vep,snpeff'''
def get_canonical_trans(dat,idd,trans_dic):
    
	canonical_trans = ''
	Canonical_soft	= 'Annovar'
	annot_flag	= True
	var_type 	= dat['Type']
	annovar_trans 	= dat['Annovar_Trans']
	vep_trans	= dat['Vep_Trans']
	snpeff_trans 	= dat['SnpEff_Trans']
	dict_key = trans_dic.keys()
	annovar_ls	= annovar_trans.split(':')
	vep_ls		= vep_trans.split(':')
	snpeff_ls	= snpeff_trans.split(':')
	
	#并不一定选择的就是NM系列的
	annovar_nm	= annovar_ls[1] if len(annovar_ls) > 2 else ''  
	(vep_nm,vep_exon)	= vep_ls[1:3] if len(vep_ls) > 3 else ['','']   
	(snpeff_gene,snpeff_nm,snpeff_exon)	= snpeff_ls[0:3] if len(snpeff_ls) > 3 else ['','','']

	canonical_nm = []
	for i in trans_dic.values():
		canonical_nm.extend(i)
	annovar_include,vep_include,snpeff_include=False,False,False
	if annovar_nm.split('.')[0] in canonical_nm: annovar_include = True
	if snpeff_nm.split('.')[0] in canonical_nm: snpeff_include = True
	if vep_nm.split('.')[0] in canonical_nm: vep_include = True

	annovar_flag = var_type == 'snv' and annovar_trans != '.'
	vep_flag = vep_nm == annovar_nm and not vep_exon.startswith('intron')
	snpeff_flag = snpeff_gene == dat['Gene.refGene'] and snpeff_nm.startswith('NM_') and 'exon' in snpeff_exon

	#!在上一版本判断策略的条件下，进行二次判断，在主策略选择未注释结果时自动校正为其他已有经典转录本的结果
	if annovar_flag :
		if annovar_include:
			canonical_source = "annovar"
		elif vep_include:
			canonical_source = "vep"
		elif snpeff_include: 
			canonical_source = "snpeff"
		else:canonical_source = "annovar"
	elif vep_flag:
		if vep_include:
			canonical_source = "vep"
		elif annovar_include:
			canonical_source = "annovar"
		elif snpeff_include:
			canonical_source = "snpeff"
		else: canonical_source = "annovar"
	elif snpeff_flag:
		if snpeff_include:
			canonical_source = "snpeff"
		elif annovar_include:
			canonical_source = "annovar"
		elif vep_include:
			canonical_source = "vep"
		else: canonical_source = "annovar"
	else:
		canonical_source = "annovar"

	if canonical_source == "annovar":
		canonical_trans = annovar_trans
		dat['Gene_c'] = canonical_trans.split(':')[0]
	elif canonical_source == "vep":
		canonical_trans = vep_trans
		Canonical_soft = "Vep"
		dat['Gene_c'] = canonical_trans.split(':')[0]
	elif canonical_source == "snpeff":
		canonical_trans = snpeff_trans
		Canonical_soft = "SnpEff"
		dat['Gene_c'] = canonical_trans.split(':')[0]
	else:canonical_trans = "."

	if '-AS1' in dat['Gene_c'] and annovar_trans == '.':
		canonical_trans = snpeff_trans
		Canonical_soft = 'SnpEff'
		dat['Gene_c'] = canonical_trans.split(':')[0]
	if '-AS1' in dat['Gene_c'] and annovar_trans != '.':
		canonical_trans = annovar_trans
		dat['Gene_c'] = canonical_trans.split(':')[0]
	'''add at 20210427 Canonical selected from three software'''
	list_nm = [annovar_nm,vep_nm,snpeff_nm]
	l_trans = [annovar_trans,vep_trans,snpeff_trans]
	list_sf = ['Annovar','Vep','SnpEff']
	n = 0
	Canonical_soft_t = Canonical_soft
 
	#!如果按照以上优先级选择出未注释结果，同时又其他软件注释出含有NM_转录本的情况，予以第三次校正
 
	if canonical_trans==".":
		if annovar_nm.startswith('NM_'):
			canonical_trans = annovar_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
			Canonical_soft = "Annovar"
		elif vep_nm.startswith('NM_'):
			canonical_trans = vep_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
			Canonical_soft = "Vep"
		elif snpeff_nm.startswith('NM_'):
			canonical_trans = snpeff_trans
			dat['Gene_c'] = canonical_trans.split(':')[0]
			Canonical_soft = "Snpeff"
	
 	#!20220812 新需求：转录本可信性排序以及转录本从众统一策略
	''' 预计实现方法：
 		从众策略：建立当前基因对应转录本字典，获取其中每个基因注释到各个转录本的频度情况，结构{"ABL2":["NM_007314.3",1],"ALK":["NM_004304.4",1]……}
				对于还未写入
    			最终经典转录本判定后，对对应的键值进行检索。如果没有对应转录本信息，入栈给初始计数1。然后在字典方法中写排序方法，将高报出的写在前面
				
   		检查经典转录本是否返回全部转录本信息，还是只有一个？设法使经典转录本方法开始包含全部转录本信息'''


	#if dat['Start'] == '106158509': print(Canonical_soft,idd,canonical_trans,list_nm[1],list_nm[2],dat['Gene'])
	if dat['Gene_c'] != '.' and dat['Gene_c'] in dict_key:
		jd = trans_dic[dat['Gene_c']]
		for nm_t in list_nm:
			n = n + 1
			if Canonical_soft == 'Annovar' and nm_t == annovar_nm and list_nm[0] in jd:continue
			if Canonical_soft == 'Vep' and nm_t == vep_nm and list_nm[1] in jd : continue
			if Canonical_soft == 'SnpEff' and nm_t == snpeff_nm and list_nm[2] in jd: continue
			for mm in range(len(jd)):
				if nm_t == jd[mm]:
					if 'intron' in l_trans[n-1] :continue
					if 'p.' not in l_trans[n-1] :continue
					canonical_trans = l_trans[n-1]
					Canonical_soft_t= list_sf[n-1]
					dat['Gene_c'] = canonical_trans.split(':')[0]
				break
			if Canonical_soft != Canonical_soft_t: break

		Canonical_soft = Canonical_soft_t



	#	if dat['Gene'] == 'TET2':print(Canonical_soft,idd,canonical_trans)
	if canonical_trans == '.' or annovar_trans == 'UNKNOWN' or 'UTR5' in annovar_trans : annot_flag = False
	dat['Gene.refGene'] = dat['Gene_c'] if  dat['Gene_c'] != '.' else dat['Gene.refGene']
	print("idd:"+idd)
	print("Canonical_soft:\t"+Canonical_soft)
	print("canonical_trans:\t"+canonical_trans)
	## 未被注释到。
	return(canonical_trans,annot_flag,Canonical_soft,dat['Gene.refGene'])

''' get Ipen_YCZ datavase dic at 20211202'''
def Ipen_YCZ_dic(f_Ipen_YCZ):
	Ipen_YCZ_gene_dic = {}
	with open(f_Ipen_YCZ,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = dict(zip(title,line.strip().split('\t')))
			Ipen_YCZ_gene_dic[dat['Gene']] =  dat['Ipen_YCZ']
	return(Ipen_YCZ_gene_dic)

''' get OMIM and OncoKB database dic at 20210512 '''
def omim_dic(f_OMIM):
	omim_gene_dic = {}
	with open(f_OMIM,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = dict(zip(title,line.strip().split('\t')))
			omim_gene_dic[dat['gene']] =  dat['phenotypes']
	return(omim_gene_dic)

def beddic(f_bed):
	dic_bed	= {}
	with open(f_bed,'r') as f:
		for line in f:
			if line.startswith('#'):continue
			dat = line.strip().split('\t')
			ch,start,end = dat[0:3]
			for i in range(int(start),int(end)+1):
				key = ch + ':' + str(i)
				dic_bed[key] = 1
	return(dic_bed)

def HGNC(f_HGNC):
	dic_HGNC = {}
	with open(f_HGNC,'r') as f:	
		title = f.readline().strip().split('\t')
		for line in f:
			dat = line.strip().split('\t')
			dic = dict(zip(title,dat))
			dic_HGNC[dic['Approved symbol']] = 1
	return(dic_HGNC)
	
			
def oncokb_dic(f_OncoKB):
	oncokb_gene_dic = {}
	with open(f_OncoKB,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = dict(zip(title,line.strip().split('\t')))
			key = ':'.join([dat['gene'],dat['Alteration']])
			oncokb_gene_dic[key] = dat['Oncogenic|Mutation Effect|Alteration Description']
	return(oncokb_gene_dic)

def Black_YCZ_dic(f_Black_YCZ):
	Black_YCZ_dic = {}
	with open(f_Black_YCZ,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = dict(zip(title,line.strip().split('\t')))
			gene = dat['Gene.refGene']
			if ':p.' in dat['AAChange.1']:
				psite = dat['AAChange.1'].split(':p.')[1]
				key = gene+':'+psite
				Black_YCZ_dic[key] = 'Y'
			else:
				print(dat['AAChange.1'])
	return(Black_YCZ_dic)

''' get ddpcr site '''
def get_ddpcr(f_ddpcr):
	ddpcr_dic = {}
	with open(f_ddpcr,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = dict(zip(title,line.strip().split('\t')))
			ddpcr_dic[dat['idd']] = 1
	return(ddpcr_dic)
	

''' get correct hgvs and white site(which Level is P/LP) '''
def get_hgvs(f_hgvs):
	hgvs_site_dic = {}
	###白名单
	with open(f_hgvs,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			if line.startswith('#'):continue
			dat = dict(zip(title,line.strip().split('\t')))
			hgvs_site_dic[dat['ID']] = [dat['Level'],dat['Canonical_One'],dat['Canonical_Three']]
	return(hgvs_site_dic)

'''20210104 强制修改部分注释'''
def get_correct_site():
	hgvs_correct_site = {}
	with open(f_correct_hgvs,'r') as inf:
		title = inf.readline().strip().split('\t')
		for line in inf:
			dat  = dict(zip(title,line.strip().split('\t')))
			hgvs_correct_site[dat['ID']] = dat['hgvs'] + ';' + dat['hgvs_correct']
	return(hgvs_correct_site)	

def get_white_site(dat,hgvs_site_dic):
	idd = '_'.join([dat['Chr'],dat['Start'],dat['End'],dat['Ref'],dat['Alt']])
	white_flag 	= False
	Canonical_One 	= dat['Canonical.transcript']
	##如果为白名单的点，则均采用白名单的hgvs格式进行输出。
	if idd in hgvs_site_dic :
		dat['Annot_flag'] 	= True
		#dat['Annot_Software']	= 'DB'
		level 			= hgvs_site_dic[idd][0]
		Canonical_One 		= hgvs_site_dic[idd][1]
		Canonical_Three 	= hgvs_site_dic[idd][2]
		if level == 'P':
			print('white site:	',hgvs_site_dic[idd])
			white_flag = True	
	return(white_flag,Canonical_One)

def get_germline_hgvs(f_germline_white_hgvs):
	dic_germline_hgvs = {}
	with open(f_germline_white_hgvs,'r') as in_f:
		title = in_f.readline().strip().split('\t')
		for line in in_f:
			item = line.strip().split('\t')
			dat = dict(zip(title,line))
			dic_germline_hgvs[dat['ID']] = [dat['Level'],dat['Canonical_One'],dat['Canonical_Three']]
	return(dic_germline_hgvs)
			
#######################################################################################################
''' get variant type '''
def get_variant_type(dat):
	(ch,start,end,ref,alt) = (dat['Chr'],dat['Start'],dat['End'],dat['Ref'],dat['Alt'])
	ref_len = len(ref)
	alt_len = len(alt)
	typ = ''
	if (ref_len == 1 and ref != '-' ) and (alt_len == 1 and alt != '-'):    # A C
		typ = 'SNV'
	elif ref == '-' or  alt.startswith(ref):    #-  A   #A  AT  #AT ATAC
		typ = 'Insertion'
	elif alt == '-' or ref.startswith(alt):     #T  -   #TA T   #TTA    TT
		typ = 'Deletion'
	else :
		typ = 'Complex'
	return(typ)


#######################################################################################################
''' filter FILTER somatic : three for two '''
def filter_FILTER_somatic(dat,n):
	fil 		= dat['Filter']
	filter_list	= fil.split('|')
	(mutect_filter,vardict_filter,varscan_filter) = filter_list
	filter_flag	= True

	count = 0
	for sub in filter_list:
		if sub == 'PASS': count+=1
	if count < n : filter_flag = False
	
	# 多等位基因捞回
	if mutect_filter == 'multiallelic' : filter_flag = True
	# 正常人中高比例位点去除
	if 'panel_of_normals' in mutect_filter : filter_flag = False
	
	return(filter_flag,count)


### correcting splicing annotation result. for del and ins, for example 'c.1889+1->GGGGG' and 'c.1889+1GGGGG>-'
def correct_splicing(dat):
	idd = '_'.join([dat['Chr'],dat['Start'],dat['End'],dat['Ref'],dat['Alt']])
	gene,trans,exon,c_site = dat['AAChange.1'].split(':')
	if dat['Type'] == 'del' and '>' in c_site:
		c_site_a = re.match(r'(c.(\d+)(\+|\-)(\d+))([A-Z]+)(>-)',c_site)
		c_site_c = list(c_site_a.groups())[0] + 'del'
		c_site = c_site_c
		AAChange = ':'.join([gene,trans,exon,c_site])
		print('Correct splicing site:',idd,'\t','before:',dat['AAChange.1'],'after:',AAChange)
	elif dat['Type'] == 'ins' and '>' in c_site:
		c_site_a = re.match(r'(c.(\d+)(\+|\-)(\d+))(->)([A-Z]+)',c_site)
		c_site_c = list(c_site_a.groups())[0] + 'ins' + list(c_site_a.groups())[5]
		AAChange = ':'.join([gene,trans,exon,c_site_c])
		print('Correct splicing site:',idd,'\t','before:',dat['AAChange.1'],'after:',AAChange)
	else:
		AAChange = dat['AAChange.1']
	return(AAChange)	

### software stats #####################################################################################
''' get MuTect2 && VarDict && VarScan stat '''
def get_stat_mutect2_vardict_varscan(title,value,ana_type):
	format_dic = {
			'mutect2': ['GT:AD:AF:DP:F1R2:F2R1:SB',
						'GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB',
						'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:SB',
						'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:PGT:PID:PS:SB',
						'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:SB',
						'GT:AD:AF:DP:F1R2:F2R1:FT:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:SB',
						'GT:AD:AF:DP:F1R2:F2R1:FT:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:PGT:PID:PS',
						'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:PGT:PID:PS:SB',
						'GT:AD:AF:DP:F1R2:F2R1:FT:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:PGT:PID:PS:SB',
			],
			'vardict': ['GT:DP:VD:AD:AF:RD:ALD',
						'GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM',
			],
			'varscan': ['GT:GQ:DP:RD:AD:FREQ:DP4',
						'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR',
			]
	}

	#GT  : AD    : AF        :DP    :F1R2:F2R1:SB
	#0/1 : 776,12: 8.703e-03 :788   :397,1 :321,3:233,543,12,0   
	#0/0 : 283,1 : 5.472e-03 :284   :118,0 :137,0:96,187,1,0

	#GT  : DP    :VD  :ALD    :RD         :AD      :AF       :BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM    
	#0/1 : 2424  :37  :11,26  :1181,1206  :2387,37 :0.0153   :2,2:7.1:1:29.5:1:0.01975:2.31403:60:36:0.0151:0:1.2   
	#0/0 : 1222  :0   :0,0    :641,580    :1221,0  :0        :0:0:0:0:0:1:0:0:0:0:0:0

	#GT  : GQ    :DP     :RD      :AD   :FREQ        :DP4    
	#0/1 : .     :318    :207     :39   :15.85%      :153,54,32,7    
	#1/1 :.      :231    :164     :18   :9.89%       :112,52,12,6

	#生成Anno信息对应的字典
	title_list 	= title.split(':')
	value_list 	= value.split(':')
	dat_dic 	= dict(zip(title_list,value_list))

	#查看Mutect2生成的标签中是否有和提供字典相匹配的标签
	if title in format_dic['mutect2']:
		#标签逗号分隔，AD:原始片段，突变片段
		(ref_depth,alt_depth) = dat_dic['AD'].split(',')
		var = '%.2f' %(float(dat_dic['AF'])*100) + '%'
		Ref_F,Ref_R,Alt_F,Alt_R = dat_dic['SB'].split(',')
	elif title in format_dic['vardict']:
		(ref_depth,alt_depth) = dat_dic['AD'].split(',')
		var = '%.2f' %(float(dat_dic['AF'])*100) + '%'
		Ref_F,Ref_R = dat_dic['RD'].split(',')
		Alt_F,Alt_R = dat_dic['ALD'].split(',')
	elif title in format_dic['varscan']:
		ref_depth = dat_dic['RD']
		alt_depth = dat_dic['AD']
		var_tmp = dat_dic['FREQ'].replace('%','') #frequency
		var = '%.2f' %(float(var_tmp)) + '%'
		if ana_type == 'Pair':
			Ref_F,Ref_R,Alt_F,Alt_R = dat_dic['DP4'].split(',')
		else:
			Ref_F,Ref_R,Alt_F,Alt_R = dat_dic['RDF'],dat_dic['RDR'],dat_dic['ADF'],dat_dic['ADR']
	else:
		raise TypeError ('The title : ' + title + ' is error, please check it')
	
	if dat_dic['GT'] == '0|0' : dat_dic['GT'] = '0/0'
	if dat_dic['GT'] == '0|1' : dat_dic['GT'] = '0/1'
	if dat_dic['GT'] == '1|1' : dat_dic['GT'] = '1/1'
	##20201216 添加Fisher Test
	fs_pvalue = stats.fisher_exact([[int(Ref_F),int(Alt_F)],[int(Ref_R),int(Alt_R)]])[1]
	ref_alt_FR = ','.join([Ref_F,Ref_R,Alt_F,Alt_R])	
	return(dat_dic['GT'],dat_dic['DP'],ref_depth,alt_depth,var,fs_pvalue,ref_alt_FR)


''' get GATK HaploTyperCaller stat '''
def get_stat_haplotypercaller(title,value):

	format_dic = [
				'GT:AD:DP:GQ:PL',
				'GT:AD:DP:GQ:PGT:PID:PL:PS',
				'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:PGT:PID:PS:SB',
				'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:OBF:OBP:OBQ:OBQRC:SB',
				'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:PGT:PID:PS:SB',
				'GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC:SB',
	]

	#GT     : AD      :DP    :GQ   :PL  
	#0/1    : 74,281  :402   :99   :7144,0,1756

	title_list 	= title.split(':')
	value_list 	= value.split(':')
	dat_dic 	= dict(zip(title_list,value_list))

	if title in format_dic:
		(ref_depth,alt_depth) = dat_dic['AD'].split(',')
		if dat_dic['DP'] == '0' : 
			var = '0.00%'
		else:
			var = '%.2f' %(float(dat_dic['AF'])*100) + '%' if 'AF' in dat_dic else '%.2f' %(int(alt_depth)/int(dat_dic['DP'])*100) + '%'
	else:
		raise TypeError ('The title : ' + title + ' is error, please check it')

	if dat_dic['GT'] == '0|1' : dat_dic['GT'] = '0/1'
	if dat_dic['GT'] == '1|1' : dat_dic['GT'] = '1/1'
	return(dat_dic['GT'],dat_dic['DP'],ref_depth,alt_depth,var)

#######################################################################################################
''' filter qc somatic single '''
def filter_qc_somatic_single(line_dic,filter_vaf,filter_depth,vaf_spec,ddpcr_site_dic):
	depth = float(line_dic['Coverage'])
	alt_depth = float(line_dic['Alt_Reads'])
	tumor_ratio = line_dic['Var']
	tumor_ratio = float ('%.4f' % (float(tumor_ratio.strip('%'))/100))
	qc_flag = True if int(depth) >= 50 and int(alt_depth) >= int(filter_depth) and tumor_ratio >= filter_vaf else False
#	lod_flag = True if int(depth) >= 50 and tumor_ratio >= vaf_spec and tumor_ratio < filter_vaf else False
	lod_flag = True if int(depth) >= 50 and tumor_ratio < filter_vaf else False
	#add at 20210525 according the RD
	if line_dic['ID'] in ddpcr_site_dic and int(depth)>=1000 and int(alt_depth)>=20:
		if tumor_ratio >= float(vaf_spec) and tumor_ratio < filter_vaf:
			lod_flag = False
			qc_flag  = True	
	return(qc_flag,lod_flag)

''' filter qc somatic pair '''
def filter_qc_somatic_pair(line_dic,filter_vaf,filter_depth,vaf_spec,ddpcr_site_dic):
	tumor_ratio	= line_dic['Var']
	normal_ratio	= line_dic['N_Var']
	tumor_ratio	= float ('%.4f' % (float(tumor_ratio.strip('%'))/100))
	normal_ratio	= float ('%.4f' % (float(normal_ratio.strip('%'))/100))
	depth		= int(line_dic['Coverage'])
	alt_depth	= int(line_dic['Alt_Reads'])
	diff 		= float(tumor_ratio) - float(normal_ratio)
	qc_flag		= False
	lod_flag	= False

	if diff <= 0:qc_flag = False;lod_flag = False
	if int(depth) >= 50 and (tumor_ratio >= filter_vaf) and (int(alt_depth) >= int(filter_depth)) : qc_flag = True
	if int(depth) >= 50 and (tumor_ratio >= float(vaf_spec)) and (tumor_ratio < filter_vaf) : lod_flag = True
	Normal_fliter = float ('%.2f' %(abs((tumor_ratio-normal_ratio)/(tumor_ratio+normal_ratio))))
	if Normal_fliter <= 0.36 : qc_flag = False;lod_flag= False
	#if normal_ratio >= float(filter_vaf) : qc_flag = False;lod_flag= False

	if line_dic['ID'] in ddpcr_site_dic and (int(depth)>=1000) and (int(alt_depth)>=20):
		if tumor_ratio >= float(vaf_spec) and tumor_ratio < filter_vaf:
			lod_flag = False
			qc_flag	 = True
	return(qc_flag,lod_flag)

''' filter qc germline '''
def filter_qc_germline(min_alt_d,depth,alt_depth,ref_depth,ratio):
	ratio = float('%.4f' % (float(ratio.strip('%'))/100))
	'''20201231 遗传提出阈值修改。过滤频率低于0.3，depth(alt+ ref) 低于50 的位点。'''
	qc_flag	 = True if int(alt_depth) >= min_alt_d else False
#	qc_flag2 = True if ratio >= 0.3 and int(alt_depth) + int(ref_depth) >= 50 and int(alt_depth) >= min_alt_d else False
#	qc_flag = True if ratio > 0.1 and int(alt_depth) >= min_alt_d else False
	return(qc_flag)

#######################################################################################################
''' filter gene '''
def filter_gene(dat,f_gene):
	gene_dic = {}
	with open(f_gene,'r') as f:
		for line in f:
			if line.startswith('#') : continue
			(gene,typ) = line.strip().split('\t')
			if 'snv' in typ :
				gene_dic[gene] = 0			
	gene_list = []
	for i in dat['Gene.refGene'].split(';') :
		if i not in gene_list : gene_list.append(i)

	if len(gene_list) == 1 : dat['Gene'] = gene_list[0].split('-AS1')[0]

	gene_flag = False
	for gene in gene_list :
		if gene in gene_dic :
			gene_flag = True
			break
	return(gene_flag,gene_dic)
	

#######################################################################################################
''' filter func '''
def filter_func(dat,synonymous_snpEff_flag):
	#print(synonymous_snpEff_flag)
	func 	    = dat['Func.refGene']
	gene_detial = dat['GeneDetail.refGene']
	
	if func == 'splicing': dat['ExonicFunc.refGene'] = func
	if dat['ExonicFunc.refGene'] == 'unknown' and gene_detial != '.' :
		dat['ExonicFunc.refGene'] = func
	
	exon_func = dat['ExonicFunc.refGene']

	''' 01 func '''
	func_flag = True if func == 'exonic' or func =='splicing' or func == 'exonic;splicing' else False
	
	''' 02 exonic func '''
	exon_flag = False if exon_func == '.' or exon_func == 'unknown' else True 	
	synonymous_flag = True if exon_func == 'synonymous SNV' else False

	''' 03 add 20201228 correct the synonymous_flag by Canonical.transcript and SnpEff_annot'''
	protein	  = dat['Canonical.transcript'].split(':')[-1]
	if 'p.' in protein:
		aa1 = re.search("(p.)([A-Z]*)([0-9]*)([A-Z]*)",protein).group(2)
		aa2 = re.search("(p.)([A-Z]*)([0-9]*)([A-Z]*)",protein).group(4)
		if '=' in protein: synonymous_flag = True
		if aa1 == aa2  :
			if aa1 != '': synonymous_flag = True
	if dat['Annot_Software'] == 'SnpEff' and synonymous_snpEff_flag: 
		synonymous_flag = True
		print('SnpEff synonymous site:',dat['SnpEff_annot'])
	return(func_flag,exon_flag,synonymous_flag)
	

#######################################################################################################
''' filter maf'''
def filter_maf(dat,maf,type_ana):
	kg_list	= ['1000g2015aug_all','1000g2015aug_eas','1000g2015aug_sas','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eur']
	exac_list 	= ['ExAC_ALL','ExAC_EAS']
	esp_list	= ['esp6500siv2_all']
	gnomad_list	= ['genome_AF','genome_AF_eas','exome_AF','exome_AF_eas']
	data_list	= kg_list + exac_list + esp_list
#	data_list 	= kg_list + exac_list + esp_list + gnomad_list
#	data_list 	= kg_list + exac_list + gnomad_list
	maf_flag 	= True
	if type_ana == 'Germline':
		for sub in data_list :
			if dat[sub] == '.' : continue
			if float(dat[sub]) >= maf :
				maf_flag = False
				break
	else:	
		for sub in data_list :
			if dat[sub] == '.' : continue
			if float(dat[sub]) > maf :
				maf_flag = False
				break
	return(maf_flag)
	
	
### database info : Pathogenic(Clinvar,HGMD) #########################################################	
def get_pathogenic():
	pathogenic_dic = {}
	with open(f_pathogenic,'r') as f:
		for line in f:
			if line.startswith('#') : continue
			dat = line.strip().split('\t')
			idd = '_'.join(dat[0:5])
			cantrans = dat[6]
			pathogenic_dic[idd] = [cantrans,dat[8],dat[9]]
	return(pathogenic_dic)


### somatic database info : Clinvar,Cosmic,CommonSnp #################################################
def filter_database_somatic(dat):
	clinvar     = dat['CLNSIG']
	common_snp	= dat['common_snp']
	snp_flag = True if common_snp == '.' else  False
	clinvar_filter_list = ['Benign','Benign/Likely_benign','Likely_benign']
	clinvar_flag = False if clinvar in clinvar_filter_list else True
	return(clinvar_flag,snp_flag)

### germline database info : Clinvar,CommonSnp #########################################################
def filter_database_germline(dat):
	clinvar     = dat['CLNSIG']
#	db_flag = False
	'''
	for item  in germline_clinvar_list:
		if clinvar == item:
			db_flag = True
			break
		else:
			db_flag = False
	'''
	snp_flag = True if dat['common_snp'] == '.' else  False
	clinvar_filter_list = ['Benign','Benign/Likely_benign','Likely_benign']
	db_flag =False if clinvar in clinvar_filter_list else True  ##  增加过滤这三情况
	return(db_flag,snp_flag)

### cosmic occurence count  ############################################################################
def get_cosmic_occurence(dat):
	cosmic_count = 0
	cosmic      = dat['cosmic88_coding']
	#ID=COSM1745359,COSM436522;OCCURENCE=1(urinary_tract),1(skin),1(lung)
	M = re.match(r'.*OCCURENCE=(.*)$',cosmic)
	if M :
		occure_list = M.group(1).split(',')	
		for sub in occure_list:
			sub = sub.strip(')')
			(sub_count,sub_name) = sub.split('(')[0:2]
			cosmic_count += int(sub_count)
	cosmic_count_flag = True if cosmic_count >= 2 else False
	return(cosmic_count_flag,cosmic_count)


### False site (before samples)#########################################################################	
''' get false site count from before sample '''
def get_pos_count(panel):
	pos_dic = {}
	with open(f_false[panel],'r') as f:
		##Gene   ID      Chr     Start   End     Ref     Alt     Count   Max     Min     Ratios  Samples
		for line in f:
			if line.startswith('#') : continue
			dat = line.strip().split('\t')
			pos_dic[dat[1]] = dat[7]
	return(pos_dic)


### Black site ########################################################################################	
''' get black(flase positive) site '''

def get_black_site(black_f,ana_type):
	black_site_dic = {}
	with open(black_f,'r') as f:
		title = f.readline().strip().split('\t')
		for line in f:
			dat = line.strip().split('\t')
			dic_l = dict(zip(title,dat))
			black_site_dic[dic_l['Var']] = dic_l['95%_upper']
	return(black_site_dic)	

def filter_black_site(dat,black_site_dic,project):
	idd = '_'.join([dat['Chr'],dat['Start'],dat['End'],dat['Ref'],dat['Alt']])
	black_flag = False
	ratio = float(dat['Var'].strip('%'))/ 100
	var = str('%.3f' %(ratio))
	if idd in black_site_dic:
		cutoff = float(black_site_dic[idd])
		if ratio < cutoff :	
			print(project+' black site: '+idd+'\t'+dat['Filter']+ ',Var:' + var + ',black_cutoff:' + str(cutoff))
			black_flag = True
	return(black_flag)
	
#######################################################################################################	
def get_cosmic_hotspot():
	hotspot_dic = {}
	with open(f_cosmic_hotspot,'r') as f:
		for line in f :
			#Chr	Pos	Ref	Alt	Gene	Protein	Context	gdna	Info
			if line.startswith('#') : continue
			dat = line.strip().split('\t')
			idd = '_'.join(dat[0:4])
			hotspot_dic[idd] = '\t'.join(dat)
	return(hotspot_dic)

def filter_cosmic_hotspot(dat,hotspot_dic):
	idd = '_'.join([dat['Chr'],dat['Start'],dat['Ref'],dat['Alt']])
	flag_hotspot = False
	if idd in hotspot_dic:
		info = dat['Var'] + '(' + dat['Alt_Reads'] + '/' + dat['Ref_Reads'] + ')' + ' ' + dat['Filter']
		print('The cosmic site : {f}  {i} : {s}'.format(f = dat['Adjust_PASS'], i = info,s = hotspot_dic[idd]))
		flag_hotspot = True
	return(flag_hotspot)


#######################################################################################################	
''' get title '''
def get_title(typ):
	if typ == 'ALL_Single':
		return(out_all_title_single)
	elif typ == 'ALL_Germline':
		return(out_all_title_germline)
	elif typ == 'ALL_Pair':
		return(out_all_title_pair)
	elif typ == 'SUB_Germline':
		return(out_sub_title_Germline)
	elif typ == 'SUB_Somatic_Pair':
		return(out_sub_title_Somatic_pair)
	elif typ == 'SUB_Somatic_Single':
		return(out_sub_title_Somatic_single)
	elif typ == 'SUB_Somatic_Pair_YCZ':
		return(out_sub_title_Somatic_pair_YCZ)
	elif typ == 'SUB_Somatic_Single_YCZ':
		return(out_sub_title_Somatic_single_YCZ)
	elif typ == 'TSO500':
		return(out_all_title_single_ts0500)
	elif typ == 'Check_Pair':
		return(out_sub_check_title_pair)
	elif typ == 'Check_Single':
		return(out_sub_check_title_single)
	else:
		raise TypeError ('%s is wrong, it must be ALL or SUB not others' %typ)	
		

########################################################################################################
''' write all filter site : to clinic'''
def write_all_filter_site(dat):
	all_filter_flag = False
	if dat['QC_flag'] and dat['Annot_flag'] and dat['Gene_flag']:
		if dat['Func_flag'] and dat['Exon_func_flag'] : all_filter_flag = True
		if dat['White_flag'] : all_filter_flag = True

	return(all_filter_flag)
	

''' write pass filter site but save LOD for 18 genes'''
def write_pass_lod_gene_site(dat):
	pass_lod_gene_flag = False

	if dat['LOD_flag'] and dat['Annot_flag'] :
		if dat['FILTER_flag'] and dat['Func_flag'] and dat['Exon_func_flag'] and \
not dat['Synonymous_flag'] and dat['MAF_flag'] :
			pass_lod_gene_flag = True
			if not dat['DB_flag'] :
				print('ClinVar False Site : ',dat['ID'],dat['CLNSIG'])
				pass_lod_gene_flag = False
		if dat['Black_flag'] : pass_lod_gene_flag = False
		if dat['White_flag'] : pass_lod_gene_flag = True
	return(pass_lod_gene_flag)

''' write pass filter site : to clinic (the last and must be true) '''
def write_pass_filter_site(dat):
	pass_filter_flag = False

	if dat['QC_flag'] and dat['Annot_flag'] and dat['Gene_flag']:
		if dat['FILTER_flag'] and dat['Func_flag'] and dat['Exon_func_flag'] and \
not dat['Synonymous_flag'] and dat['MAF_flag'] :
			pass_filter_flag = True
			if not dat['DB_flag'] :
				print('ClinVar False Site : ',dat['ID'],dat['CLNSIG'])
				pass_filter_flag = False

		if dat['Black_flag'] : pass_filter_flag = False
		if dat['White_flag'] : pass_filter_flag = True
	return(pass_filter_flag)	


''' write pass fite : for TMB analysis '''
def write_pass_site(dat):
	pass_flag = False
	
	if dat['White_flag'] or dat['FILTER_flag'] : pass_flag = True
	if dat['Black_flag'] : pass_flag = False
	return(pass_flag)

	
''' write all title out '''
def write_out(dat,title_list):
	data_filter_list = []
	for title in title_list:
		if title in dat:
			data_filter_list.append(str(dat[title]))
		else:	
			raise TypeError ('<%s> is not in level title' %title)
	return('\t'.join(data_filter_list) + '\n')

