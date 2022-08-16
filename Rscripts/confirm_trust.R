args_reader <- commandArgs()

#table_file_in <- args[6]
#write_file_out <- args[7]
args_reader <- as.array(args_reader)

print_help <- function(){
  cat("\n")
  cat("\t------------------------------------------------\n")
  cat("\tFunction: \tconfirm the low rate transcript id to the high rate one\n")
  cat("\tAuthor:liujiaxing\t date:20220815\n")
  cat("\tUsage: \tRscript confirm_trust.R input_file_name output_file_name\n")
  cat("\tinput_file_name: \txls file which needs to reorder\n")
  cat("\toutput_file_name: \toutput table file name, sep =<tab>\n")
  cat("\t-------------------------------------------------\n")
  cat("\n")
}

#添加了外显子数量不同时的从众策略（这种情况在实际数据分析中还发现的挺多）
get_right_exon <- function(running_table){
  pos_exon_info <-  running_table$Affected_Exon
  transcript_info <- running_table$Transcript
  uniq_transcript <- unique(transcript_info)
  
  total_exon <- c()
  affect_exon <- c()
  for(i in pos_exon_info){
    buff_e <- strsplit(i,split = "/")[[1]]
    if(length(buff_e)==2)
      {
      affect_exon <- append(affect_exon,buff_e[1])
      total_exon <- append(total_exon,buff_e[2])
    }else{
        affect_exon <- append(affect_exon,".")
        total_exon <- append(total_exon,".")
      }
  }
  
  for(ut in uniq_transcript){
    index_array <- which(transcript_info==ut)
    total_exon_this_transc <- total_exon[index_array]
    table_exon <- sort(table(total_exon_this_transc),decreasing = T)
    right_total_exon <- names(table_exon)[1]
    total_exon[index_array] <- right_total_exon
  }
  
  right_total_exon <- paste(affect_exon,total_exon,sep = "/")
  return(right_total_exon)
}

#相同基因，只要几种软件注释中报出更多频度转录本信息，即进行从众选择，提取报出数量更多的转录本
#频度依据从高到低的顺序进行匹配，一旦匹配，则进行数据修改，跳出当前位点
running_upper <- function(table_file_in, write_file_out) 
{
  reader <- read.table(table_file_in,header = T,sep = "\t",quote = "",stringsAsFactors = F)
  print(dim(reader))
  outer <- data.frame()
  gene_list <- reader$Gene
  unique_gene <- unique(gene_list)
  
  for (gene in unique_gene){
    df_subset <- reader[which(gene_list==gene),]
    transcript_pd <- sort(table(df_subset$Transcript),decreasing = T)
    if(length(transcript_pd>1)){
      #正在运行的基因频度表格长度
      running_tab = length(transcript_pd)
      #从大到小（低频-高频）循环
      while(running_tab){
        #得到当前运行频度优先级转录本在子表格中的匹配索引
        running_index = which(df_subset$Transcript==names(transcript_pd)[running_tab])
        for(ind in running_index){
          #输出三个软件给出的转录本信息
          annovar_trans_raw = df_subset$Annovar_Trans[ind]
          vep_trans_raw = df_subset$Vep_Trans[ind]
          snpeff_trans_raw = df_subset$SnpEff_annot[ind]
          
          annovar_trans = "."
          vep_trans = "."
          snpeff_trans = "."
          
          if(annovar_trans_raw!=".")
          {
            annovar_trans = strsplit(annovar_trans_raw,split = ":")[[1]][2]
          }
          if(vep_trans_raw!=".")
          {
            vep_trans = strsplit(vep_trans_raw,split = ":")[[1]][2]
          }
          if(snpeff_trans_raw!=".")
          {
            snpeff_trans = strsplit(df_subset[ind,]$SnpEff_annot,split = ":")[[1]][2]
          }
          #按频度从前到后的顺序检查三个转录本结果是否匹配到更高频度的结果中
          for(j in 1:running_tab){
            if(annovar_trans==names(transcript_pd)[j])
            {
              df_subset[ind,]$Transcript = names(transcript_pd)[j]
              break
            }
            else if(vep_trans==names(transcript_pd)[j])
            {
              df_subset[ind,]$Transcript = names(transcript_pd)[j]
              break
            }
            else if(snpeff_trans==names(transcript_pd)[j])
            {
              df_subset[ind,]$Transcript = names(transcript_pd)[j]
              break
            }
          }
        }
        #提示频度优先排序
        running_tab = running_tab-1
      }
    }
    #print(gene)
    outer <- rbind.data.frame(outer,df_subset)
  }
  
  colnames(outer) <- colnames(reader)
  get_right_exon(outer) -> right_exon_array
  #print(right_exon_array)
  outer$Affected_Exon <- right_exon_array
  write.table(outer,quote = F,file = write_file_out,row.names = F,sep="\t")
}


#主流程入口

if (file.exists(args_reader[6]))
{
  running_upper(table_file_in=args_reader[6],write_file_out= args_reader[7])
}else{
  print_help()
  cat("Input file not find!\n")
}
