import os
import re
import sys
import argparse

sys.path.append(os.path.dirname(__file__))
sys.path.append('/data/Pipeline/PMC/TISSUSE/Database/pylib/')
import annotFilterCommon as aFC
from annovarTrans import annovar_deal
from vepTrans import vep_deal, vep_deal_test
from snpEffTrans import snpEff_deal
from amino_acid_trans import one2three
from TransVar import Trans

''' main '''
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Filter GATK MuTect && VarDict && VarScan Somatic Annot Result')
	parser.add_argument('-s', '--sample',help = 'sample name',required = True)
	parser.add_argument('-a', '--ana_type',help = 'sample model, pair or single',required = True)
	parser.add_argument('-vcf','--vcf',help = 'vcf file',required = True)
	parser.add_argument('-g', '--gene_list',help = 'panel gene list',required = True)
	parser.add_argument('-vep', '--vep_annot',help = 'vep annot file',action='store',default = False)
	parser.add_argument('-eff', '--snpeff_annot',help = 'snpeff annot file',action='store',default = False)
	parser.add_argument('-fv','--filter_vaf',help = 'filter vaf for report,default=0.02',default = 0.02)
	parser.add_argument('-fd','--filter_depth',help = 'filter alt  depth for report,default=8',default = 8)
	parser.add_argument('-fv_s','--vaf_spec',help = 'filter vaf for report,default=0.01',default = 0.01)
	parser.add_argument('-fd_s','--depth_spec',help = 'filter alt depth for report,default=19',default = 19)
	parser.add_argument('-p','--project',help = 'project for 27 or 599 group',default = '27')	
	parser.add_argument('-IPY', '--Ipen_YCZ',help = 'Ipen_YCZ file',required = True)	
	parser.add_argument('-b', '--black',help = 'black file',required = True)

	annot_files_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'annot_files')
	f_hgvs		= os.path.join(annot_files_dir,'hg19_variants_hgvs.txt')
	args		= vars(parser.parse_args())
	sample		= args['sample']
	ana_type	= args['ana_type']
	in_vcf		= args['vcf']
	f_gene		= args['gene_list']
	f_vep_vcf	= args['vep_annot']
	f_eff_vcf	= args['snpeff_annot']
	filter_vaf	= float(args['filter_vaf'])
	filter_dep	= int(args['filter_depth'])
	vaf_spec	= float(args['vaf_spec'])
	depth_spec	= int(args['depth_spec'])
	project		= args['project']
	f_Ipen_YCZ	= args['Ipen_YCZ']
	f_black		= args['black']
	f_HGNC		= aFC.f_HGNC
	f_InPanel	= aFC.f_InPanel
	f_InExon	= aFC.f_InExon

	s_dir 		= os.path.dirname(os.path.dirname(os.path.abspath('.')))
	tumor_name 	= sample + '_T' if ana_type == 'Pair' else sample
	recal_bam 	= os.path.join(s_dir,tumor_name,'Realign',tumor_name + '.recal.bam')
	print('Filter Vaf   :',filter_vaf)
	print('Filter Depth :',filter_dep)
	
	''' 打开待写入的文件 共6个 '''
	#File 01: all variants and all information and flag
	f_all_out		= open(sample + '.Somatic.All.Annot.xls','w')	
	#File 02: all variants except qc_fail and fun_fail variants, give out
	f_all_filter_out	= open(sample + '.Somatic.All.Filter.Annot.xls','w')
	#File 03: all pass variants
	f_pass_out		= open(sample + '.Somatic.PASS.Annot.xls','w')		
	#File 04: all pass variants and p/lp variants ,then filter the qc,func and maf, give out
	f_pass_filter_out	= open(sample + '.Somatic.PASS.Filter.Annot.xls','w')
	#File 05: all pass variants and p/lp variants,then filter the qc,func and maf,but 0.01<=vaf <0.02
	f_pass_lod_gene_out	= open(sample + '.Somatic.lod_gene.Check.xls','w')
	#File 06: add id and sample info for f_pass_filter_out
	f_report_check_out	= open(sample + '.Somatic.Check.xls','w')	

	''' 输出各个表格的表头 '''
	out_all_title 	= aFC.get_title('ALL_' + ana_type)
	out_sub_title 	= aFC.get_title('SUB_Somatic_' + ana_type)
	out_check_title = aFC.get_title('Check_' + ana_type) 

	
	f_all_out.write('\t'.join(out_all_title) + '\n')
	f_all_filter_out.write('\t'.join(out_sub_title) + '\n')
	f_pass_out.write('\t'.join(out_all_title) + '\n')
	f_pass_filter_out.write('\t'.join(out_sub_title) + '\n')
	f_pass_lod_gene_out.write('\t'.join(out_check_title) + '\n')
	f_report_check_out.write('\t'.join(out_check_title) + '\n')

	''' get VEP annot result '''
	vep_dic = vep_deal_test(f_vep_vcf) if f_vep_vcf != False else {}

	''' get snpEff annot result '''
	eff_dic = snpEff_deal(f_eff_vcf,aFC.trans_dic) if f_eff_vcf != False else {}
	
	''' get site count samples before '''
	pos_27_dic 	= aFC.get_pos_count('27')
	pos_599_dic	= aFC.get_pos_count('599')
	
	''' get pathogenic site '''
	dm_dic = aFC.get_pathogenic()

	''' get cosmic hostspot site '''
	cosmic_hotspot_dic = aFC.get_cosmic_hotspot()

	''' get white site '''
	hgvs_site_dic = aFC.get_hgvs(f_hgvs)

	''' get black site '''
	black_site_dic = aFC.get_black_site(f_black,project)

	''' get hgvs_correct site '''
	hgvs_correct_site = aFC.get_correct_site()

	''' get OncoKB database '''
	oncokb_gene_dic = aFC.oncokb_dic(aFC.f_OncoKB)

	''' get ddpcr site '''
	ddpcr_site_dic	= aFC.get_hgvs(aFC.f_ddpcr)

	''' get Black_YCZ database '''
	Black_YCZ_dic = aFC.Black_YCZ_dic(aFC.f_Black_YCZ)
	
	''' get Ipen_YCZ database'''
	Ipen_YCZ_dic = aFC.Ipen_YCZ_dic(f_Ipen_YCZ)

	dic_HGNC = aFC.HGNC(f_HGNC)
	dic_InPanel = aFC.beddic(f_InPanel)
	dic_InExon  = aFC.beddic(f_InExon)
	dic_tmb	    = aFC.beddic(aFC.f_tmb)
	### main process ######################################################################################

	dic_annovar = {
		"synonymous SNV"	    :"Synonymous_substitution",
		"stopgain"		    :"Nonsense_substitution",
		"startloss"		    :"Nonsense_substitution",
		"nonsynonymous SNV"	    :"Missense_substitution",
		"frameshift insertion"	    :"Frameshift_insertion",
		"frameshift deletion"	    :"Frameshift_deletion",
		"nonframeshift deletion"    :"Inframe_deletion",
		"nonframeshift insertion"   :"Inframe_insertion",
		"splicing"		    :"Splice_Site_mutation",
		"stoploss"		    :"Elongation_mutation",
		"frameshift substitution"   :"Truncation_mutation",
		"nonframeshift substitution":"Missense_substitution/Synonymous_substitution",
		"unknown"		    :"Other",
		"."			    :"Other"}

	with open(in_vcf, 'r') as f:
		res 		= ''
		title_tmp   = f.readline().strip().split('\t')
		title	    = []

		for item in title_tmp:
			if 'WithVer' in item:
				title_aa = re.sub('WithVer','',item)
				title.append(title_aa)
			else:
				title.append(item)

		title_len   = title.index('Otherinfo1') + 1
		for line in f:
			dat 		= line.strip().split('\t')
			line_dic 	= dict(zip(title,dat))
			ch_len		= len(line_dic['Chr'].split("_"))
			if ch_len > 1: continue
			idd = '_'.join([line_dic['Chr'],line_dic['Start'],line_dic['End'],line_dic['Ref'],line_dic['Alt']])
			line_dic['ID'] 		= idd
			line_dic['Sample']	= sample
			line_dic['IDD'] 	= '_'.join([sample,line_dic['ID']])
			line_dic['Check'] 	= ''
			line_dic['Bam'] 	= recal_bam
	
			#if line_dic['Ref'] == '0' or line_dic['Alt'] == '0' : continue
			#if len(line_dic['Ref']) > 50 or len(line_dic['Alt']) > 50 : continue
			ch,start,end = line_dic['Chr'],int(line_dic['Start']),int(line_dic['End'])
			mm,nn,jj = 0,0,0
			#查看突变区间是否的每个碱基是否与热点区域重合，并统计重合位点数量
			for i in range(start,end+1):
				key = ch + ':' + str(i)
				if key in dic_InPanel:
					mm = mm + 1
				if key in dic_InExon:
					nn = nn + 1
				if key in dic_tmb:
					jj = jj + 1

			#Raise三个指标的Flag
			if mm == 0 :
				line_dic['InPanel'] = 'No'
			if mm == (end+1-start):
				line_dic['InPanel'] = 'Yes'
			if mm < (end+1-start) and mm != 0:
				line_dic['InPanel'] = 'PartlyIn'
			if nn == 0 :
				line_dic['InExon'] = 'No'
			if nn == (end+1-start):
				line_dic['InExon'] = 'Yes'
			if nn < (end+1-start) and nn != 0:
				line_dic['InExon'] = 'PartlyIn'
			if jj == 0:
				line_dic['InTMB'] = 'No'
			if jj == (end+1-start):
				line_dic['InTMB'] = 'Yes'
			if jj < (end+1-start) and jj != 0:
				line_dic['InTMB'] = 'PartlyIn'

			############################# 01 basic information #########################################
			''' get variant type '''
			''' 1>1 SNP - > n Insert n > - Deletion'''
			line_dic['Type'] = aFC.get_variant_type(line_dic)
	
			''' get FILTER info '''
			#写入三款软件是否PASS的指标
			line_dic['Filter'] = dat[title_len -1 + 9]
			(line_dic['MuTect2'],line_dic['VarDict'],line_dic['VarScan']) = line_dic['Filter'].split('|')
			## just one or more softwares pass is pass, filter panel_of_normals , save the pseudo-multiallelic
			#小于规定个数的软件得到此位点予以去除，返回是否通过以及通过软件个数信息
			(line_dic['FILTER_flag'],line_dic['PASS_Count']) = aFC.filter_FILTER_somatic(line_dic,1) 	
			
			############################# 02 stat #####################################################
			''' get tumor and normal sample stats and qc filter '''
			stats_title = dat[title_len -1 + 11]
			stats_tumor = dat[title_len -1 + 12]
			
			#stats_title是分析类型，例如ALL_Single，依据类型给出字典
			tumor_stat_list = aFC.get_stat_mutect2_vardict_varscan(stats_title,stats_tumor,ana_type)
			#得出基因型，覆盖率百分比，原始、突变读数，BS（不带百分号的百分比值），dp4（匹配到4个位置的reads数量）
			line_dic['GT'],line_dic['Coverage'] = tumor_stat_list[0],tumor_stat_list[1]
			line_dic['Ref_Reads'],line_dic['Alt_Reads'],line_dic['Var'] = tumor_stat_list[2:5]
			line_dic['BS_pvalue'],line_dic['DP4'] = tumor_stat_list[5],tumor_stat_list[6]

			if ana_type == 'Pair':
				stats_normal     = dat[title_len -1 + 13]
				normal_stat_list = aFC.get_stat_mutect2_vardict_varscan(stats_title,stats_normal,ana_type)
				line_dic['N_GT'], line_dic['N_Coverage'] = normal_stat_list[0],normal_stat_list[1]
				line_dic['N_Ref_Reads'],line_dic['N_Alt_Reads'],line_dic['N_Var'] = normal_stat_list[2:5]
				line_dic['N_BS_pvalue'],line_dic['N_DP4'] = normal_stat_list[5],normal_stat_list[6]

				''' filter qc '''
				#LOD：log odds ratio（对数差异比）
				line_dic['QC_flag'],line_dic['LOD_flag'] = aFC.filter_qc_somatic_pair(line_dic,filter_vaf,filter_dep,vaf_spec,ddpcr_site_dic)

			else:
				line_dic['Clinvar_Origin'] 	= 'na'
				line_dic['ChosenMed_Status'] 	= 'na'
				
				''' filter qc '''
				from annotFilterCommon import filter_qc_somatic_single
				line_dic['QC_flag'],line_dic['LOD_flag'] = aFC.filter_qc_somatic_single(line_dic,filter_vaf,filter_dep,vaf_spec,ddpcr_site_dic)

			############################# 03 datbase ####################################################
			''' get dm '''
			#如果pathogenic数据库中存在记录则写入记录，否则写入点
			line_dic['DM'] = dm_dic[line_dic['ID']] if line_dic['ID'] in dm_dic else '.'
			
			''' get database(cosmic,clinvar,commonsnp) info '''
			#与pathogenic类似
			(line_dic['DB_flag'],line_dic['SNP_flag']) = aFC.filter_database_somatic(line_dic)	
				
			''' get cosmic occurence count,cosmic_count >=2 Cosmic_flag is Ture '''
			#与pathogenic类似
			(line_dic['Cosmic_flag'],line_dic['cosmic_count']) = aFC.get_cosmic_occurence(line_dic)
			
			''' add 20190813 black site filter '''
			#是否在黑名单字典中，升起黑名单标志
			line_dic['Black_flag'] 	= aFC.filter_black_site(line_dic,black_site_dic,project)

			''' add 20190829 sample site before '''
			#在之前的成熟Panel列表中检出的个数
			line_dic['Count_27'] 	= pos_27_dic[idd] if idd in pos_27_dic else 0
			line_dic['Count_599'] 	= pos_599_dic[idd] if idd in pos_599_dic else 0
			
			########################## 04 canonical transcript ##########################################
			''' get Annovar canonical transcript '''
			aa_change 	= line_dic['AAChange.refGene']
			#if aa_change == 'UNKNOWN' or re.findall('wholegene',aa_change) : continue
			#idd是包含染色体信息，位置起止信息，ref、alt序列以下划线分割的字符串。之前vcf字典化是以该值作为key的
			#最新的和之前的注释信息
			(line_dic['Annovar_Trans'],line_dic['Annovar_Trans_Before']) = annovar_deal(line_dic,idd,aFC.trans_dic)[0:2]	

			''' add 20190715 get vep annotaion result '''
			#在每一个软件生成的字典中检索当前行字典名称是否存在与对应的vcf注释文件字典中
			if idd in vep_dic:
				(line_dic['HGVS_OFFSET'],line_dic['Vep_Trans'],line_dic['Vep_Trans_Three'])=vep_dic[idd][0:3]
				gene				=  vep_dic[idd][-2]
				line_dic['Vep_Trans']		= gene + ':' + line_dic['Vep_Trans']
				line_dic['Vep_Trans_Three']	= gene + ':' + line_dic['Vep_Trans_Three']
			else:
				(line_dic['Vep_Trans'],line_dic['Vep_Trans_Three']) = ('.','.')
				
			''' add 20190715 get snpeff annotaion result '''
			''' add 20201228 snpeff synonymous_snpEff_flag '''
			if idd in eff_dic:
				(line_dic['SnpEff_annot'],line_dic['SnpEff_Trans'],line_dic['SnpEff_Trans_Three'],synonymous_snpEff_flag)= eff_dic[idd][0:4]
			else:
				(line_dic['SnpEff_annot'],line_dic['SnpEff_Trans'],line_dic['SnpEff_Trans_Three'],synonymous_snpEff_flag)= ('.','.','.','.')

				
			''' add 20190910 makesure canonical transcript according three  '''
			#使用典型转录本校正结果
			(line_dic['Canonical.transcript'],line_dic['Annot_flag'],line_dic['Annot_Software'],line_dic['Gene.refGene']) = aFC.get_canonical_trans(line_dic,idd,aFC.trans_dic)
			
			''' liujiaxing add 20220811 makesure Gene.Func mapping canonical transcript '''
			if line_dic['Annot_Software'] == 'Vep': line_dic['Gene.Func'] = vep_dic[idd][-1]
			if line_dic['Annot_Software'] == 'SnpEff': line_dic['Gene.Func'] = eff_dic[idd][-1]
			
			''' add 20190813 the correct hgvs for some site(updating...) and white site() '''
			(line_dic['White_flag'],line_dic['Canonical.transcript']) =aFC.get_white_site(line_dic,hgvs_site_dic)
			line_dic['AAChange.1'] = line_dic['Canonical.transcript']
			
			################# 05 filter ,correct the Synonymous_flag for some site(add 20201228) ################
			''' filter gene list '''
			#Gene_flag是检测当前dict里面有哪些基因是在提供的-g参数文件中的，在的话返回在的基因列表和True标签
			line_dic['Gene_flag'],gene_dic 	 = aFC.filter_gene(line_dic,f_gene)	
			#if not line_dic['Gene_flag'] : continue
			''' filter func '''
			#
			(line_dic['Func_flag'],line_dic['Exon_func_flag'],line_dic['Synonymous_flag']) = aFC.filter_func(line_dic,synonymous_snpEff_flag)
			
			''' filter maf in population'''
			#各个人种基本突变频率大于0.01则升起MAF的False标签
			line_dic['MAF_flag'] = aFC.filter_maf(line_dic, 0.01,'')

			''' correct the hgvs provided by ycz begun at 20210105'''
			if idd in hgvs_correct_site.keys():
				hgvs_list = hgvs_correct_site[idd].split(';')
				hgvs_before, hgvs_after = hgvs_list[0],hgvs_list[1]
				if line_dic['Canonical.transcript'] == hgvs_before: 
					line_dic['Canonical.transcript'] = hgvs_after
					line_dic['AAChange.1']	= hgvs_after

			''' correct the gene name at 20210125'''
			if len(line_dic['Gene.refGene'].split(';')) > 1 and line_dic['Canonical.transcript'] != '.':
				gene_correct = line_dic['AAChange.1'].split(':')[0]
				print('Multigene name:',idd,'\t',line_dic['Gene.refGene'],'\t',gene_correct,'\t',line_dic['AAChange.1'],'\t','Please Check it!!!')
				line_dic['Gene.refGene'] = gene_correct
			line_dic['Gene_flag'] = True if line_dic['Gene.refGene'] in gene_dic else False

			''' correct the splicing site at 20210406 '''
			if line_dic['Func.refGene'] == 'splicing' and len(line_dic['AAChange.1'].split(":")) == 4:
				aa_change = line_dic['AAChange.1']
				line_dic['AAChange.1'] = aFC.correct_splicing(line_dic)

			''' correct for some sepcial site # BJ20CM002544''' 
			if line_dic['AAChange.1'] == 'TP53:NM_000546:exon10:c.1044_1045delinsTT:p.348_349delinsF*':
				line_dic['AAChange.1'] = 'TP53:NM_000546:exon10:c.1044_1045delinsTT:p.348_349delinsF' 
			if line_dic['AAChange.1'] == 'RECQL:NM_002907:exon2:c.2T>C:p.M1?':
				line_dic['AAChange.1'] = 'RECQL:NM_002907:exon2:c.2T>C:p.M1T'
			if line_dic['AAChange.1'] == 'APC:NM_001127511:exon1:c.1A>G:p.M1?':
				line_dic['AAChange.1'] = 'APC:NM_001127511:exon1:c.1A>G:p.M1V'
			if line_dic['AAChange.1'] == 'SH2B3:NM_005475:exon2:c.473_475delinsCCC:p.HT158PP':
				line_dic['AAChange.1'] = 'SH2B3:NM_005475:exon2:c.473_475delinsCCC:p.158_159delinsPP'
			if line_dic['AAChange.1'] == 'CHEK2:NM_007194:exon11:c.1116_1117delinsTG:p.K373E':
				line_dic['AAChange.1'] = 'CHEK2:NM_007194:exon11:c.1116_1117inv:p.K373E'
			line_dic['Canonical.transcript'] = line_dic['AAChange.1']

			line_dic['Description_OncoKB'] = '.'
			line_dic['Black_YCZ'] = '.'
			if ":p." in line_dic['Canonical.transcript']:
				p_site = line_dic['Canonical.transcript'].split('p.')[1]
				key_OncoKB = line_dic['Gene.refGene'] + ':' + p_site
				if key_OncoKB in oncokb_gene_dic:
					line_dic['Description_OncoKB'] = oncokb_gene_dic[key_OncoKB]
				if key_OncoKB in Black_YCZ_dic:
					line_dic['Black_YCZ'] = 'Y'

			line_dic['Ipen_YCZ'] = '-'
			if line_dic['Gene.refGene'] in Ipen_YCZ_dic:
				line_dic['Ipen_YCZ'] = Ipen_YCZ_dic[line_dic['Gene.refGene']]

			line_dic['Gene']			= line_dic['Gene.refGene']
			line_dic['VAF']				= line_dic['Var'].split("%")[0]
			line_dic['该突变cut-off值'] = '0.02'

		#	if line_dic['ExonicFunc.refGene'] in dic_annovar.keys():
		#		line_dic['Consequence'] = dic_annovar[line_dic['ExonicFunc.refGene']]
		#	else:
		#		line_dic['Consequence'] = line_dic['ExonicFunc.refGene']

		#####################fd################################################
			line_dic['Transcript_ann'],line_dic['cHGVS_ann'],line_dic['pHGVS_ann'],line_dic['Affected_Exon_ann']= '.','.','.','.'
			Canonical_list_ann = line_dic['Annovar_Trans'].split(":")
			if line_dic['Annovar_Trans'] != '.':
				if len(Canonical_list_ann) == 5:
					line_dic['Transcript_ann'] = Canonical_list_ann[1]
					line_dic['Affected_Exon_ann'] = Canonical_list_ann[-1]
					line_dic['cHGVS_ann'] = Canonical_list_ann[3]
				if len(Canonical_list_ann) == 6:
					line_dic['Transcript_ann'] = Canonical_list_ann[1]
					line_dic['Affected_Exon_ann'] = Canonical_list_ann[-1]
					line_dic['cHGVS_ann'] = Canonical_list_ann[3]
					line_dic['pHGVS_ann'] =one2three(line_dic['Annovar_Trans'].split(":")[-2])
	###########################################################################



			line_dic['Transcript'],line_dic['cHGVS'],line_dic['pHGVS'],line_dic['Affected_Exon']= '.','.','.','.'
			Canonical_list = line_dic['Canonical.transcript'].split(":")
			if line_dic['Canonical.transcript'] != '.':
				if len(Canonical_list) == 5: 
					line_dic['Transcript'] = Canonical_list[1]
					line_dic['Affected_Exon'] = Canonical_list[-1]
					line_dic['cHGVS'] = Canonical_list[3]
				if len(Canonical_list) == 6:
					line_dic['Transcript'] = Canonical_list[1]
					line_dic['Affected_Exon'] = Canonical_list[-1]
					line_dic['cHGVS'] = Canonical_list[3]

			if line_dic['Canonical.transcript'] != '.' and len(Canonical_list) != 1:
				if line_dic['Annot_Software'] == 'Annovar':
					if len(Canonical_list) == 6:
						#将一个字母的氨基酸转变为3个字母
						line_dic['pHGVS'] =one2three(line_dic['Canonical.transcript'].split(":")[-2])
					if Canonical_list[1] in line_dic['Vep_Trans_Three']:
						if len(Canonical_list) == 6:
							line_dic['pHGVS'] = line_dic['Vep_Trans_Three'].split(":")[-2]
					elif Canonical_list[1] in line_dic['SnpEff_Trans_Three']:
						if len(Canonical_list) == 6:
							line_dic['pHGVS'] = line_dic['SnpEff_Trans_Three'].split(":")[-2]
				if line_dic['Annot_Software'] == 'Vep':
					line_dic['Transcript'] = line_dic['Vep_Trans_Three'].split(":")[1]
					line_dic['pHGVS'] = one2three(line_dic['Canonical.transcript'].split(":")[-2])
					if len(Canonical_list) == 6:
						line_dic['pHGVS'] = line_dic['Vep_Trans_Three'].split(":")[-2]
				if line_dic['Annot_Software'] == 'SnpEff':
					line_dic['Transcript'] = line_dic['Transcript']
					if len(Canonical_list) == 6:
						line_dic['pHGVS'] = line_dic['SnpEff_Trans_Three'].split(':')[-2]
				if '%' in line_dic['pHGVS']:
					line_dic['pHGVS'] = re.sub("%3D","=",line_dic['pHGVS'])
				if len(Canonical_list) == 6 and line_dic['ExonicFunc.refGene']== 'synonymous SNV' and "=" not in line_dic['pHGVS']:
					tmp_a = re.match(r'(p.)([A-Za-z]+)([0-9]+)([A-Za-z]+)',line_dic['pHGVS'])
					if tmp_a:
						line_dic['pHGVS'] = tmp_a.groups()[0]+tmp_a.groups()[1]+tmp_a.groups()[2]+'='

			if line_dic['ExonicFunc.refGene'] in dic_annovar.keys():
				if line_dic['ExonicFunc.refGene']=='nonframeshift substitution'and'='not in line_dic['pHGVS'] and not line_dic['Synonymous_flag']:
					line_dic['Consequence'] = "Missense_substitution"
				elif line_dic['ExonicFunc.refGene']=='nonframeshift substitution'and('=' in line_dic['pHGVS'] or line_dic['Synonymous_flag']):
					line_dic['Consequence'] = "Synonymous_substitution"
				else:
					line_dic['Consequence'] = dic_annovar[line_dic['ExonicFunc.refGene']]
			elif 'splicing' in line_dic['ExonicFunc.refGene']:
				line_dic['Consequence'] = 'Splice_Site_mutation'
			else:
				line_dic['Consequence'] = line_dic['ExonicFunc.refGene']

			line_dic['Chr'] = line_dic['Chr']
			line_dic['CHR'] = line_dic['Chr'].split("chr")[1]
			line_dic['Gene'] = line_dic['Gene.refGene']
			line_dic['HGNC'] = 'Yes' if line_dic['Gene'] in dic_HGNC else 'No'
			line_dic = Trans(line_dic)
			############################## 06 write out ########################################################
			''' 06-1 write all and filtered variants to output : give for report '''
			all_filter_flag = aFC.write_all_filter_site(line_dic)
			if all_filter_flag :
				f_all_filter_out.write(aFC.write_out(line_dic,out_sub_title))

			''' 06-2 write the pass and filtered variants to output : give for report '''
			pass_filter_flag = aFC.write_pass_filter_site(line_dic)
			if pass_filter_flag == True:
				f_pass_filter_out.write(aFC.write_out(line_dic,out_sub_title))
				f_report_check_out.write(aFC.write_out(line_dic,out_check_title))
			
			line_dic['Adjust_PASS'] = pass_filter_flag
			''' 06-3 write the pass and filtered variants(var < 0.02) to output'''
			pass_lod_gene_flag = aFC.write_pass_lod_gene_site(line_dic)
			if pass_lod_gene_flag == True:
				f_pass_lod_gene_out.write(aFC.write_out(line_dic,out_check_title))
			''' 06-4 write all the variants to output : check for analysis '''
			f_all_out.write(aFC.write_out(line_dic,out_all_title))

			''' 06-5 write pass variants to output : for TMB analysis '''
			if line_dic['FILTER_flag'] and not line_dic['Black_flag'] and line_dic['Gene_flag']:
				f_pass_out.write(aFC.write_out(line_dic,out_all_title))

			#########核查是否遗漏#########################################################################
			''' add 20200227 cosmic hotsplot site '''
			hotspot_flag = aFC.filter_cosmic_hotspot(line_dic,cosmic_hotspot_dic)

