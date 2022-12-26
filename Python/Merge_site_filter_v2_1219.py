import argparse as ar
import re
import sys
import Dealpos
import numpy as np
import os

def write_del(file_pointer,line_dict):
	'''
	写出被删除的记录，方便查阅
	
	Args:
	file_pointer (file_write_buffer) : 文件写出指针
	line_dict (dict) :待写出的删除行字典
	'''
	
	write_list = []
	for i in line_dict.values() : write_list.append(str(i))
	write_line = "\t".join(write_list)+"\n"
	file_pointer.write(write_line)


def get_pass_flag(soft_flag):
	"""判断已经Call出的位点质量，若质量较高仍旧合并

	Args:
		soft_flag (str): Call 出软件的质量信息字符串

	Returns:
		bool: 是否为高质量位点
	"""
	flag_list = soft_flag.split("|")
	count_pass = 0
	for flag in flag_list:
		if flag == "PASS":
			count_pass+=1
	
	if count_pass >=2:
		return True
	return False


def table_diric(vcf_name):
	"""将Merge.vcf经Annovar注释后的表格文件txt字典化

	Args:
		vcf_name (str): Annovar.TXT文件名

	Returns:
		dict : annovar注释文件解析成的表格，依据是否有Merge标签分离成合并前后两个字典
	"""
	dict_raw = dict()
	dict_new = dict()
	with open(vcf_name,'r') as reader:
		index_split = reader.readline().strip().split("\t")
		for line in reader:
			line_split = line.split(sep="\t")
			line_dict = dict(zip(index_split,line_split))
			#转换
			idd_list = []
			if len(line_dict["Chr"].split("_")) >1:continue
			for i in ["Chr","Start","End","Ref","Alt"]:idd_list.append(line_dict[i])
			idd = "_".join(idd_list)
			if "Merge=Yes" in line:
				dict_new[idd] = line_dict
			else: dict_raw[idd] = line_dict
	return dict_raw,dict_new


def trans_right_idd(idd4):
	idd4_split = idd4.split("_")
	_,idd4_split = Dealpos.Dealpos(idd4_split[0],idd4_split[1],idd4_split[2],idd4_split[3])
	idd4_split = idd4_split.split("\t")
	if idd4_split[2] == "-":end_pos = idd4_split[1]
	else:end_pos = int(idd4_split[1])+len(idd4_split[2])-1
	list(idd4_split).insert(2,str(end_pos))
	idd5 = "_".join(idd4_split)
	return idd5


def get_split_site(info_string):
	split_buffer = info_string.split(";")
	idd4_list = []

	for i in split_buffer:
		if not i.startswith("MUTORI"):continue
		mutline = re.sub(pattern="MUTORI=",repl="",string=i)
	split_site = mutline.split("chr")
	del split_site[0]
	for ss in split_site:
		second_split = ss.split(",")
		if len(second_split) > 1:
			#多个等位突变
			ch_start_ref = second_split[0]
			idd4_list.append("chr"+re.sub(pattern=r"[.]",repl="_",string=ch_start_ref))
			ch_start_ref = ch_start_ref.split('.')[:-1]
			ch_start_ref = "_".join(ch_start_ref)
			ch_start_ref = "chr"+ch_start_ref
			for i in range(0,len(second_split)-2):
				anoth_idd = ch_start_ref+"_"+second_split[i+1]
				idd4_list.append(anoth_idd)
			pass
		else:
			#非多等位情况
			second_split = str(second_split[0])
			second_split = re.sub(pattern=r'[.]',repl="_",string=second_split)
			second_split = "chr"+second_split
			idd4_list.append(second_split)
	for i in range(0,len(idd4_list)): idd4_list[i] = re.sub(pattern=r'[*]',repl="-",string=idd4_list[i])

	return idd4_list


def extre_diff(mapping_list,total_dict,remain_dict,backup_file,report_line,diff_cutoff=5):
	"""对极端频率差距的合并位点进行鉴别
		
		设计思路：
			仅针对按p.判断应该保留合并位点的情况，所以只给出应该保留的分位点信息，修改原字典返回
		策略：
  			频率排序，如相邻频率差大于0.05的，设置断点。断点前为频率相近的位点，断点后位点合并到总匹配中。
			最终保留断点前的所有位点和合并后的位点。断点前合并位点深度应减去匹配到合并位点中的数量，重新计算深度
			获取各个位点深度（AD=ref,alt） ，建立临时字典，决定要去除哪些低频位点
			合并后位点频率对齐到分位点最低频，如果分位点中有报出线下且能的，会去除本应报出的分位点，此种情况不应合并
	
	Args:
		mapping_list(list) : 只含有匹配后的idd列表
		total_dict (dict) : 注释后的全局字典
		remain_dict (dict) : 当前未运行记录字典，用于识别是否被多个位点包含，防止提前删除
		diff_cutoff (float) : 频率断点最大差异值
		report_line (float) : 报出线，防止合并拉低临界报出低频位点频率到报出线下（默认0.01）
		backup_file (IOTextWarpper) : 写出文件指针
		
	return:
		total_dict (dict) , breakpoint (int) : 修正数值之后的总字典,切断状态提示（>=0断点位置，-1没有断点，其余特殊情况置空不处理）
	
	"""
	# line_list= [ref,alt]
	ar_list = []
	var_list = []
	total_dict_keys = total_dict.keys()

	for keys in mapping_list:
		if keys in total_dict_keys: 
			line_dict = total_dict[keys]
			ar_list.append(int(line_dict['Alt_Reads']))
			var = float(re.sub(pattern=r'%',repl="",string=line_dict['Var']))
			var_list.append(var)
		else:
			#if keys not in remain_dict.values():mapping_list.remove(keys)
			continue

	#如果var最大值高于报出线，最小值低于报出线，
	#为防止合并位点频率与最低分位点对齐无法通过过滤器，不合并，低频分位点会被其它过滤器自动去掉
	report_line = report_line*100
	for var in var_list:
		if max(var_list)> report_line and var<report_line:
			print("分位点位于报出线两侧，不合并")
			return total_dict,-2
		
	#依据Var列表对键值做局部排序，确定保留哪些
	var_list = np.array(var_list,dtype=float)
	sort_mappinglist = []
	ar_mappinglist = []
	#以数值化后的var列表顺序对idd列表和Alt_reads列表进行排序
	for i in np.argsort(var_list): 
		sort_mappinglist.append(mapping_list[i])
		ar_mappinglist.append(ar_list[i])
	# breakpoint 得到频率跳跃断点

	breakpoint = -1
	for i in range(0, len(var_list)-1):
		if abs(float(var_list[i]) - float(var_list[i+1])) > float(diff_cutoff):
			print("存在频率断点，对位点保留及频率进行重新整合",end="\t")
			breakpoint = i
			break
	if breakpoint == -1: 
		#没有断点，合并前位点频率相近，全部去除
		return total_dict,breakpoint
	else :
		#存在断点，删除断点之后的，合并前位点的alt_reads和Var重新计算
		for idd in sort_mappinglist[breakpoint-1:] :
			if idd not in total_dict.keys():continue
			line_raw = total_dict[idd]
			raw_alt = int(line_raw['Alt_Reads'])
			raw_cov = int(line_raw['Ref_Reads']) + int(line_raw['Alt_Reads'])
			raw_alt = raw_alt-ar_mappinglist[breakpoint]
			#断点之后的即使降到0以下也没关系，反正会移除
			if raw_alt <= 0 and idd not in sort_mappinglist[:breakpoint+1]:
				print("\tERROR:\t"+idd+" 深度异常，需要人工检查！\t")
				return total_dict,-2
			raw_var = raw_alt/int(raw_cov)
			raw_var = str(round(raw_var*100,2))+"%"
			line_raw['Var'] = raw_var
			line_raw['Alt_Reads'] = raw_alt
			total_dict[idd] = line_raw
		for idd in set(sort_mappinglist[:breakpoint]):
			print("移除断点后位点："+idd)
			if idd not in remain_dict.values() and idd in total_dict.keys():
				write_del(file_pointer=backup_file,line_dict=total_dict[idd])
				total_dict.pop(idd)
	return total_dict,breakpoint
	

def get_include(new):
	'''
	:func:返回两个位点包含关系字典
	:param raw: 分立位点字典
	:param new: 合并后位点字典
	:return: 合并前后位点包含关系字典
	'''

	include_dict = dict()
	for idd in new:
		line_dic = new[idd]
		oth_info = line_dic["Otherinfo11"]
		if re.findall(pattern=r'MUTORI=;',string=oth_info):continue
		if len(line_dic["Chr"].split("_")) >1:continue
		idd4_list = get_split_site(oth_info)
		for i in range(0,len(idd4_list)):
			idd_run = str(idd4_list[i])
			idd4_list[i] = trans_right_idd(idd_run)
		include_dict[idd] = idd4_list

	return include_dict


def deal_table(table_file_name):
	"""针对命令行中运行模式，将合并注释文件字典化

	Args:
		table_file_name (str): All.Annot.xls文件名称

	Returns:
		dict : All.Annot.xls文件转化成的字典
	"""

	outer = dict()
	tab = open(table_file_name,'r')
	dict_index = tab.readline().strip().split("\t")
	for line in tab:
		line_split = line.strip().split("\t")
		line_dict = dict(zip(dict_index,line_split))
		idd = "_".join([line_dict["Chr"],line_dict["Start"],line_dict["End"],line_dict["Ref"],line_dict["Alt"]])
		if idd not in outer.keys() : outer[idd] = line_dict
	tab.close()
	return outer


def get_aa_pos(trans_info):
	'''
	Parameters: 从经典转录本中返回氨基酸位置和可能影响片段的起始位置
	input:  经典转录本字符串
	output: AA positions (int format),probable start aa pos
	'''
	aa_pos = "0"
	forward_pos = 0
	if "p." not in trans_info or re.findall(pattern=r'p.$',string=trans_info):
		return 0,0
	canon_split = trans_info.split(":")
	if len(canon_split)==5:
		aa_pos = canon_split[4]
	else :return 0,0
	if not aa_pos.strip().split('_')[0] == aa_pos:
		forward_pos = aa_pos.strip().split("_")[0]
		forward_pos = re.sub(r'^\D+',repl="",string=forward_pos)
		aa_pos = aa_pos.strip().split('_')[1]
	aa_pos = re.sub(r'^\D+',repl="",string=aa_pos)
	aa_pos = re.search(pattern="^\d+",string=aa_pos).group(0)
	aa_pos = int(aa_pos)
	print("Trans info: "+ str(trans_info),end="\t")
	forward_pos = int(forward_pos)
	return aa_pos,forward_pos

def get_aa_multi_flag(f_aa_pos,f_pob_pos,s_aa_pos,s_pob_pos):
	"""判断两个相邻位点是否影响同一个氨基酸，或者影响的氨基酸区段是否存在交集

	Args:
		f_aa_pos (int): 前一个位点影响氨基酸的起始位置 
		f_pob_pos (int): （仅针对突变影响多个氨基酸，否则为0）  前一个位点影响氨基酸的终止位置
		s_aa_pos (int): 后一个位点影响氨基酸的起始位置
		s_pob_pos (int): （仅针对突变影响多个氨基酸，否则为0）  后一个位点影响氨基酸的终止位置

	Returns:
		bool: 相邻位点相互影响关系（相互影响为True）
	"""
	'''
	Input:两个p.的影响位置，及可能存在的影响区段的起始位置
	Output:是否具有生物学意义的Flag
	'''
	return_flag = False
	#都不是影响肽段区间，判断位置
	if not f_pob_pos and not s_pob_pos:
		if f_aa_pos==s_aa_pos:return_flag = True
	#前一个是区间，可能包含后一个
	if f_pob_pos and not s_pob_pos:
		f_range = range(f_pob_pos,f_aa_pos)
		if s_aa_pos in f_range:return_flag=True
	#前一个是单点，可能被后一个包含
	if s_pob_pos and not f_pob_pos:
		s_range = range(s_pob_pos,s_aa_pos)
		if f_aa_pos in s_range: return_flag=True
	#两个都是区间，可能会有交集
	if s_pob_pos and f_pob_pos:
		f_range = range(f_pob_pos, f_aa_pos)
		s_range = range(s_pob_pos, s_aa_pos)
		for i in f_range:
			if i in s_range:return_flag = True
	#如果任意一个位置没有注释到p. 则不应该合并
	if not f_aa_pos and not s_aa_pos:
		return_flag = False

	return return_flag


def merge_site_filter(total_dict, vcf_name,cutoff=5,backup_file="",report_line=0.01):
	"""嵌入流程中主要的入口函数

	Args:
		total_dict (dict): 行字典打包成的总字典，键值为idd
		vcf_name (str): merge.vcf 由Annovar注释的文件名
		cutoff (int): 判断频率断点的差异百分比

	Returns:
		dict: 过滤后的总字典，键值为idd
	"""

	# 给合并后位点最小优先级，意义未明一律不报

	# 获取配对信息
	raw, new = table_diric(vcf_name=vcf_name)
	mapping_dict = get_include(new=new)
	remain_dict = mapping_dict.copy()

	if backup_file=="":
		_,sample_name = os.path.split(vcf_name)
		sample_name = sample_name.split(".")[0]
		baclup_file = sample_name+".del_annot.xls"
		backup_file = os.path.join(os.path.dirname(vcf_name),baclup_file)
	bk_file = open(backup_file,'w')

	raw_idd_list = list(raw.keys())
	new_idd_list = list(new.keys())
	total_dict_keys = list(total_dict.keys())
	for mapping_key in mapping_dict:
		if mapping_key not in total_dict_keys: 
			#print(mapping_key+" is missing")
			continue
		raw_site_list = mapping_dict[mapping_key]
		raw_site_list.sort()
		# 如果只有一条匹配，弹出合并前记录，保留合并后结果
		if len(raw_site_list) == 1 and raw_site_list[0] in total_dict.keys():
			print(mapping_key+" 只有一个分位点，无法判断信息，不合并!")
			if mapping_key in new_idd_list: new_idd_list.remove(mapping_key)
			#if raw_site_list[0] in raw_idd_list: new_idd_list.remove(raw_site_list[0])
			continue

		# 多位点合并情况
		print(mapping_key, end="\t")
		print("：", end="")
		print(raw_site_list, end="\t")
		exist_flag = True
		for raw_site in raw_site_list:
			if raw_site not in total_dict.keys(): 
				print(" 部分合并前记录无法匹配，不合并！", end="\n")
				exist_flag = False
				break

		if not exist_flag: continue # 如果不是所有位点都在记录中存在，不好判断，不予合并
		
		#没有N_Var的Tumor only不适用于该策略，直接跳过
		#N_Var大于30%的视为胚系位点，从合并关系中移除（对应合位点是超高频体系位点合并真假存疑）
		for i in raw_site_list.copy():
			if "N_Var" not in total_dict[i].keys(): break
			nvar = total_dict[i]["N_Var"]
			if nvar == ".": break
			if float(re.sub(pattern=r'%',repl="",string=nvar)) > 30.0:
				print("胚系位点：{} 从合并关系中移除！".format(i),end=" ")
				raw_site_list.remove(i)
				
		if len(raw_site_list) <1:
			continue
		if total_dict[raw_site_list[0]]["Ref"] == "-" or total_dict[raw_site_list[0]]["Alt"] == "-" or len(total_dict[raw_site_list[0]]['Ref']) != len(total_dict[raw_site_list[0]]['Alt']):
			# 判断Indel
			print("第一个位点 " + raw_site_list[0] + " 是/包含Indel",end=" ")
			total_dict,break_flag = extre_diff(mapping_list=raw_site_list,total_dict=total_dict,remain_dict=remain_dict,diff_cutoff=cutoff,report_line = report_line,backup_file = bk_file)
			if break_flag == -1:
				print("分位点频率相近,合并")
				for raw_site in set(raw_site_list):
					if raw_site in raw_idd_list: raw_idd_list.remove(raw_site)
			continue

		# 20220906 方法改进，增加全局判断
		include_flag = True
		for k in range(0, len(raw_site_list) - 1):
			fk_dict = total_dict[raw_site_list[k]]
			sk_dict = total_dict[raw_site_list[k + 1]]
			fk_aa_pos, fk_pob_pos = get_aa_pos(fk_dict["Canonical.transcript"])
			sk_aa_pos, sk_pob_pos = get_aa_pos(sk_dict["Canonical.transcript"])
			# 判断是否包含，给flag
			running_flag = get_aa_multi_flag(f_aa_pos=fk_aa_pos, f_pob_pos=fk_pob_pos,
												 s_aa_pos=sk_aa_pos, s_pob_pos=sk_pob_pos)
			include_flag = running_flag and include_flag

		if include_flag:
			if mapping_key in raw_site_list:
				print("软件已Call出，",end=" ")
				#存在特殊情况，s1:[s2,s1] s2:[s1,s2] 第一次检索合并到s1，删除s2。导致第二次检索s2位点的时候报错
				if mapping_key not in total_dict.keys():
					print("但位点已被其它call出位点合并，仅保留第一次合并结果")
					continue
				called_bk = dict(total_dict[mapping_key])
				#total_dict.pop(mapping_key)
				total_dict,breakpoint = extre_diff(mapping_list=raw_site_list,total_dict=total_dict,remain_dict=remain_dict,diff_cutoff=cutoff,report_line = report_line,backup_file=bk_file)
				total_dict[mapping_key] = called_bk
				new_idd_list.remove(mapping_key)
				if breakpoint == -1:
					print("分位点频率相近，合并")
					for index in set(raw_site_list):
						if index != mapping_key and index in raw_idd_list: raw_idd_list.remove(index)
				continue
					# 存在被两个位点同时影响的氨基酸位点，去除raw,保留Merge
			print("相邻突变均影响同一氨基酸位置",end=" ")
			total_dict,break_flag = extre_diff(mapping_list=raw_site_list,total_dict=total_dict,remain_dict = remain_dict,diff_cutoff=cutoff,report_line = report_line,backup_file = bk_file)
			if break_flag==-1:
				print("，合并",end="\n")
				for j in raw_site_list:
					if j in raw_idd_list: raw_idd_list.remove(j)
				continue
		else:
			#PASS个数
			if mapping_key in total_dict.keys():
				pass_flag = get_pass_flag(total_dict[mapping_key]['Otherinfo10'])
			else: pass_flag = False

			if mapping_key in raw_site_list:

				#高质量位点已经Call出（Filter中大于两个PASS的），无条件合并高质量Call出位点
				if pass_flag:
					print("软件Call出高质量合位点",end="\t")
					total_dict,break_flag = extre_diff(mapping_list=raw_site_list,total_dict=total_dict,remain_dict=remain_dict,diff_cutoff=cutoff,report_line = report_line,backup_file = bk_file)
					if break_flag == -1:
						print("不存在频率断点，合并")
						for j in raw_site_list:
							if j in raw_site_list: raw_idd_list.remove(j)
						continue

				print("虽然软件已经Call出，但不影响同一氨基酸，去除Call出位点和合并位点")
				if mapping_key in total_dict.keys(): 
					write_del(bk_file,line_dict=total_dict[mapping_key])
					total_dict.pop(mapping_key)
			else:
				print("存在相邻突变不影响同一氨基酸，不合并")
				if mapping_key in new_idd_list: new_idd_list.remove(mapping_key)
		remain_dict.pop(mapping_key)

	out_dict_new, out_dict_raw = dict(), dict()
	for i in raw_idd_list:
		if i in total_dict.keys(): out_dict_raw[i] = total_dict[i]
	for i in new_idd_list:
		if i in total_dict.keys(): out_dict_new[i] = total_dict[i]
	
	#删掉的记录写回到备份文件中
	for i in raw.keys():
		if i not in raw_idd_list and i in total_dict.keys():write_del(file_pointer=bk_file,line_dict=total_dict[i])
	for i in new.keys():
		if i not in new_idd_list and i in total_dict.keys():write_del(file_pointer=bk_file,line_dict=total_dict[i])
	bk_file.close()

	out_dict = out_dict_new.copy()
	out_dict.update(out_dict_raw)
	return out_dict


if __name__ == "__main__":
	'''
	设计思路：
	1、获取Merge vcf文件中的Merge标签，区分新老位点
	2、解析上游合并脚本给出的MUTLINE合成关系，合并成字典
	3、判断合并后位点合理性，保留new记录则删除raw记录，反之亦然：
		3.1 第一个为非移码indel必须保留
		3.2 前后位点在注释中索引，影响同一个氨基酸的保留
		3.3 20220909策略更新，任何Indel都予以保留
	4、先写入raw对应的注释表格记录，后写入new对应的注释表格记录
	
	有关复杂包含关系：1+2=>3、3+4=>5等，树形判断。
	推断：
	1. 定义可知任意树的分支如不能合并，则上级的枝节点不能合并。合并位点出现先后顺序对过滤的结果无影响。
	2. 上级枝节点合并后删除下级枝节点，加判断如存在进下级枝节点的枝节点则删除下级的下级，不理会下级枝节点是否存在。（删除操作是最高优先级）位点先后顺序对保留的结果无影响。
	
	上级枝节点                  5
							 /   \ 
	下级枝节点               3     4
						   /  \
	下级枝节点的枝节点      1   2
	'''
	parser = ar.ArgumentParser(description="判断合并位点意义并更改写入记录")
	parser.add_argument("--anovar_table","-t",required=True,help="包含合并位点信息的annovar注释表格名称【annovar.txt】")
	parser.add_argument("--annot_table",'-a',required=True,help="合并位点注释后的文件【*.All_Annot.xls】")
	parser.add_argument("--breakpoint_cutoff","-b",required=False,help="频率断点百分比，默认为5",default=5)
	parser.add_argument("--backup_file","-bk",required=False,help="删除记录的备份文件名，默认为样本名+del_annot.xls",default="")
	parser.add_argument("--report_line","-r",required=False,help='报出线 百分之0.5 输入  0.005',default=0.01)
	parser.add_argument("--out_file",'-o',required=True,help="写出的文件全名 [Required]")
 
	try:
		arg = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	annovar_name = arg.anovar_table
	annot = arg.annot_table
	out_name = arg.out_file
	cutoff = arg.breakpoint_cutoff
	report_line = arg.report_line

	if arg.backup_file !="" :
		backup_file = arg.backup_file
	else:
		_,abs_filename = os.path.split(annot)
		sample_name = abs_filename.split('.')[0]
		backup_file = sample_name + ".del_annot.xls"

	#生成合并后字典、合并前字典和合并前后匹配字典
	#读取annovar表格的目的是获取Merge=Yes的标签和MUTLINE合成关系
	
	annot_dict = deal_table(annot)
	filted_dict = merge_site_filter(total_dict=annot_dict, vcf_name=annovar_name,cutoff=cutoff,backup_file=backup_file,report_line = report_line)

	writer = open(out_name,'w')
	header_col = open(annot,'r').readline()
	header_index = list(header_col.strip().split('\t'))
	header_str = header_col
	if header_str[-1]!='\n':
		writer.write(header_str+"\n")
	else: writer.write(header_str)
	
	for idds in filted_dict:
		line_dict = filted_dict[idds]
		write_str = ""
		for innerkeys in header_index:
			write_str = write_str+"\t"+str(line_dict[innerkeys])
		write_str = re.sub(pattern="^\t",repl="",string=write_str)
		print(write_str,file=writer)
	writer.close()
