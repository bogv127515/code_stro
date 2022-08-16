#export PERL5LIB=/software/Vep/src/ensembl/modules:/software/Vep/VEP_plugins-release-96:/software/Vep/src/ensembl-io/modules:/software/Vep/src/ensembl-variation/modules:/software/Vep/src/ensembl-funcgen/modules:/software/Vep/src/ensembl-compara/modules:${PERL5LIB} 

HELP:
	@echo Description:
	@echo Usage:
	@echo make -f makefile_zp CONFIG=*
	@echo "******************************"
	@echo Source:Liujiaxing
	@echo Editer:Liujiaxing
	@echo Date:2022/08/02
	@echo "******************************"

include $(CONFIG)

Sample_name=$(sample_name)
Project=$(project_name)
In_dir=$(input_dir)
Out_dir=$(output_dir)
Module=$(module)
	

ALL: Annovar SnpEff Vep

Vep:
	$(VEP) --cache --offline --everything --dir /software/Vep/cache --refseq --port 3337 --fasta /software/Vep/fasta/Homo_sapiens.GRCh37.dna.toplevel.fa --fork 4 --force_overwrite --no_stats --pick -i $(Out_dir)/$(Sample_name).all.vcf -o $(Out_dir)/$(Sample_name).Vep.annot.vcf

Annovar:
	$(PERL) $(table_annovar) $(Out_dir)/$(Sample_name).all.vcf $(ANNOVAR_DB) -outfile $(Out_dir)/$(Sample_name).Merage.all.Somatic.vcf.annovar --otherinfo -vcfinput -protocol refGeneWithVer,refGeneRegion,cytoBand,avsnp150,snp138,common_snp,clinvar_20220320,cosmic88_coding,cosmic_hotspots,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,exac03,gnomad211_genome,gnomad211_exome,esp6500siv2_all,cadd13gt20,gerp++gt2,ljb26_all,intervar_20180118,wgRna,targetScanS,gwasCatalog,rmsk -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -nastring . -polish --remove --thread 10 -buildver hg19 

SnpEff:
	$(JAVA) -Xmx4g -jar $(SnpEff) -v hg19 -hgvs $(Out_dir)/$(Sample_name).all.vcf -strict -noInteraction -noMotif -onlyProtein > $(Out_dir)/$(Sample_name).SnpEff.annot.vcf

Filter:
	@echo "======== $(Project) Annot and Filter Start at `date` ========"
	cp $(In_dir)/$(Sample_name).all.vcf $(Out_dir)/$(Sample_name).all.vcf
	sed -i '/^chr[0-9]*_/d' $(Out_dir)/$(Sample_name).all.vcf
	#删除特殊染色体
	$(PYTHON) /data/Project/liujiaxing/zp_pipeline/packages/annot/annotFilterSomatic.py -s $(Sample_name) -a $(Module) -vcf $(Out_dir)/$(Sample_name).Merage.all.Somatic.vcf.annovar.hg19_multianno.txt -vep $(Out_dir)/$(Sample_name).Vep.annot.vcf -eff $(Out_dir)/$(Sample_name).SnpEff.annot.vcf -fv 0.02 -fd 8 -g /fass3/BMC-to-GC/ZP/2022-somatic/bed/Illumina_gene_u.list -p 599 -IPY /data/Database/Drugable/ChosenMed/hg19_Ipen_YCZ.list -b /data/Pipeline/PMC/TISSUSE/Database/BlackList/Chosen599/IDT_599.sample_VS_var.mutation.xls
	@echo "======== $(Project) Annot and Filter Start at `date` ========"
	touch $(Out_dir)/$(Sample_name).end
	mv $(Sample_name).* $(Out_dir)/
	mv snpEff* $(Out_dir)
	sleep 1
