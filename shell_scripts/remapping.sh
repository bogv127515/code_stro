if [ -d "./tmp" ]; then
    rm -r ./tmp
    echo "Removing previous temporary directory !"
else
    mkdir -p ./tmp
    echo "Creating new temporary directory !"
fi

GATK4=/data/pipe_panel/software/GATK4/gatk-4.1.9.0/gatk
java_options='"-XX:ParallelGCThreads=10 -Xmx20G -Djava.io.tmpdir=./tmp"'
fasta=/data/pipe_panel/database/Genome/hg19/bwa_index/ucsc.hg19.fasta

chr=$(echo "$2" | cut -d':' -f1)
pos=$(echo "$2" | cut -d':' -f2)
start_pos=$(($pos-1000))
end_pos=$(($pos+1000))
echo -e "${chr}\t${start_pos}\t${end_pos}\ttest_gene">./dup.bed

${GATK4} HaplotypeCaller -L ./dup.bed -R ${fasta} -I $1 -ploidy 2 -ERC GVCF --bam-output remap.germline.bam --create-output-bam-md5 True -O remap.auto.g.vcf --tmp-dir tmp

rm -r ./tmp
rm dup.bed
