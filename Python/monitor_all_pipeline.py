import pymysql
import re
import argparse
import datetime
import os
import send_email_myself

__author__ = "liujiaxing"
__mail__ = "jiaxingliu@chosenmedtech.com"
__date__ = "2022/11/24"

fo5 = "cat /data/Pipeline/PMC/TISSUSE/Task_Monitor/pipeline/others/stat_sample/所有*订单统计.csv /data/Pipeline/PMC/TISSUSE/Task_Monitor/pipeline/others/stat_sample/all_fat_localdata_order.xls /data/Pipeline/PMC/TISSUSE/Task_Monitor/pipeline/others/stat_sample/liaorui_order.xls|grep -a"


def get_and_write(sql_cmd,file_name):
	"""通过输入命令行获取数据库状态，然后把满足条件的写在缓冲文件中

	Args:
		sql_cmd (_string_): 需要执行的命令
		file_name (_string_) : 写出匹配记录的文件名

	Returns:
		_int_: case(-1,0,int>0)
			-1 : 获取信息失败等错误
			0 : 没有新的消息，可以直接跳过
			int>0 : 存在新消息
	"""

	try:
		db = pymysql.connect(host="192.168.1.122", port=3306, user='bmc', passwd='ChosenMed#2020', db='bmc', charset='utf8')
	except:
		print("Connect to SQL ERROR!")
		return -1

	#执行pySQL，返回为满足记录的Tuple，没有记录返回长度为1的元组
	cursor = db.cursor()
	cursor.execute(sql_cmd)
	result = cursor.fetchall()
	db.commit()
	db.close()
	if len(result) == 0:
		print("No records found !")
		return 0

	outer = open(file_name, 'a',encoding='utf-8')
	for i in range(len(result)):
		line = "\t".join(map(str, result[i]))
		outer.write(line+"\n")
	
	outer.close()
	return len(result)

def deal_table(file_name):
	"""解析从数据库中拉取的表格，返回字典

	Args:
		file_name (string): 缓冲文件名称

	Returns:
		dict: 表格文件解析成的字典，为节约内存，只写入时间、样本名称
	"""
	reader = open(file_name,'r')
	header = ["Date","Product","ID"]
	out_dict = {}
	for line in reader:
		line_split = line.split("\t")
		sample_id = line_split[3]
		sample_date = line_split[0]
		product_line = line_split[7]
		line_dict = dict(zip(header, [sample_date,product_line,sample_id]))
		out_dict[sample_id] = line_dict

	return out_dict

def deal_time(time,max_days,report_hours):
	"""处理输入的时间字符串，返回是否应报告时间窗口异常

	Args:
		time (string): SQL拉下表格中时间列字符串
		max_days (int): 最大日期差，有些未及时写入的时间过长的也不再检索报告
		report_hours (string): 最小报出时间 格式： hh:mm:ss

	Returns:
		bool : Flag 为true则有超期未报出 为false则无异常 
	"""
	# 有些很古早的样本 格式适用 2020-09-26-15-40-42 统一去除
 
	if len(time.split("-")) == 6:
		return False

	record_time = datetime.datetime.strptime(time,"%Y-%m-%d %H:%M:%S")
	now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	now = datetime.datetime.now().strptime(now,"%Y-%m-%d %H:%M:%S")

	diff_time = str(now - record_time)
	if "day" in diff_time:
		diff_days = int(diff_time.split(" ")[0])
	else:
		diff_time = "0 day, "+diff_time
		diff_days = 0

	if diff_days > max_days:
		#超过最大分析天数，不再追踪
		return False
	elif diff_days != 0:
		#超期一天以上
		return True

	print(diff_time)
	diff_hour = int(diff_time.split(" ")[2].split(":")[0])
	if diff_hour >= report_hours :
		return True

	return False


def get_finish_flag(sample_name):
    #为True表示已经完成，下同
	cmd = "{fo5} {sample}|cut -f3 ".format(sample = sample_name,fo5 = fo5)
	analysis_path = os.popen(cmd=cmd).readline()
	if not len(analysis_path) ==0 :
		print("发现分析文件夹未创建")
		return False
	record_file = "{path}/shell/log.txt".format(path = analysis_path)
	if not os.path.isfile(record_file):
		return False
	last_line_cmd = "tail -n1 " + record_file
	last_line = os.popen(last_line_cmd).readline()
	if "done" not in last_line:return False

	return True


def report_window(all_dict,product_line,max_days,report_hours):
	"""在一定时间窗口内报出所有未交付样本
		设计思路：
		所有未交付样本解析 => 过滤时间窗口 => 过滤完成flag
	
	Args:
	all_dict (dict): 解析SQL表格的字典，包含样本名、产品线名称、Pending时间


	Returns:
	list: 通过任何过滤器，在预警时间窗口内，真正应该注意到的样本

	"""
	
	#方便理解，拆分两个过滤流程，避免实时检索，一次性打包，创建深拷贝
	all_dict_keys = list(all_dict.keys())

	for key in all_dict_keys:
		line_dict = all_dict[key]
		time_string = line_dict['Date']
		if not deal_time(time=time_string,max_days=max_days,report_hours=report_hours) :
			all_dict.pop(key)
			continue

		line_dict = all_dict[key]
		product_desc = line_dict['Product']
		in_present_product = False
		for cel_product in product_line:
			if re.findall(pattern=cel_product,string=cel_product):
				in_present_product = True
		if not in_present_product: 
			print("{number} not present in product lines.".format(number=product_desc))
			continue
		
		#in_present_product 出现在现有产品中，产品类型固定不为空
		if get_finish_flag(sample_name=line_dict['ID']):
			all_dict.pop(key)
			continue
	
	return all_dict

def mail_to_myself(undone_samples):
    pass

if __name__ == "__main__":
	"""整体思路：
		解析多个SQL表格，全部写到一个文件中去 => 解析表格、筛选时间窗口 => 有时间窗口内的样本，每个产品设计方法检查完成flag => 有完成Flag的删除记录，没有的输出表格
		发送含有表格的邮件，邮件内容为未完成的订单名称+注意检查啥的
	"""
 
	parser = argparse.ArgumentParser("SQL 选择测试")
	parser.add_argument("--list", "-l", help="需要监控的产品线类型，应与库表中后缀名称一致",dest="list",nargs="+",required=True)
	parser.add_argument("--flag", "-f", help="待筛选的订单 Stat 标签",required=True)
	parser.add_argument("--verbose", "-v", help="是否反向选择 deflut=True ",required=False,type=bool,default=True)
	args = parser.parse_args()

	monit_line = args.list
	flag = args.flag
	out_file = "buff.xls"

	#现有产品线，后续扩充的话可以添加，汉字部分为对应关键词，只要可以唯一匹配即可
	product_line_sql_name = ['BLCAsubtype',"127_HT","blood",'bloodct','TME']
	product_line = ["膀胱癌分子分型","遗传性肿瘤","血液肿瘤综合","ctDNA","肿瘤免疫微环境"]
	product_dict = dict(zip(product_line_sql_name,product_line))

	for product in monit_line:
		if product not in product_dict.keys():
			print("Input product {product} cannot found in product lines".format(product = product))
			continue
		db_name = product
		#获取数据状态，写入缓冲文件
		if not args.verbose:
			sql = 'SELECT * FROM sampleInfo_{db} where State="{flag}"'.format(db=db_name, flag = flag)
		else:
			sql = 'SELECT * FROM sampleInfo_{db} where State!="{flag}"'.format(db=db_name, flag = flag)
		status = get_and_write(sql_cmd=sql, file_name=out_file)
		#print("Get records success from {}".format(product))


	checking_sample = deal_table(file_name=out_file)
	#给5个小时的分析窗口，7天的追溯极限
	report_sample = report_window(all_dict=checking_sample, product_line=product_dict,max_days=7,report_hours=5)
	
	if len(report_sample) > 0:
		report_sample_keys = list(report_sample.keys())
		for i in report_sample_keys:
			if re.findall(pattern="HE",string=i) and i in report_sample.keys():
				report_sample.pop(i)
				continue
			print("Sample {0}\tPending time: {1} \t Product:{2}".format(i,report_sample[i]['Date'],report_sample[i]['Product']))
	
	if len(report_sample) > 0:
		report_string = str("\t".join(report_sample.keys()))
		if os.path.isfile("sample_error") : send_email_myself.sendMail(report_string = report_string)
		os.popen("touch sample_error")
	else:
		now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		print("Checking all sample done at {time_string}  異常なし!".format(time_string=now))

	os.popen("rm {buff}".format(buff=out_file))

