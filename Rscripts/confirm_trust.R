args <- commandArgs()

table_file_in <- args[6]
write_file_out <- args[7]

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
write.table(outer,quote = F,file = write_file_out,row.names = F,sep="\t")
