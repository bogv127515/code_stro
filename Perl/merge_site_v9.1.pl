#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long; 
use Text::NSP::Measures::2D::Fisher::left;

my ($distance,$vcf,$bam,$id);
GetOptions(
		"help|h"        => \&usage,
		"distance|d:s"  => \$distance,
		'vcf:s'         => \$vcf,
		'bam:s'         => \$bam,
		'id:s'          => \$id,
);
&usage if (not $bam or not $vcf or not $id);

$distance   ||= 5;
my $hg19     = "/data/Database/Genome/hg19/bwa_index/ucsc.hg19.fasta";
my $samtools = "/software/samtools/samtools-1.9/samtools";
my $bedtools = "/software/bedtools/bedtools2/bin/bedtools";
my $sortbed  = "/software/bedtools/bedtools2/bin/sortBed";

system("$sortbed -i $vcf >$id.sort.vcf");
system("$samtools stats --reference $hg19 --threads 10 $bam >$id.stats");
my $total_bases = my $mismatch_base = 1;

open FF,"$id.stats";
while(<FF>){
	if(/SN\s+bases\s+mapped\s+\(cigar\):\s+(\d+)/){$total_bases=$1}
	if(/SN\s+mismatches:\s+(\d+)/){$mismatch_base=$1}
}
#print "$total_bases-$mismatch_base\n";
close FF;

my $n=1;
my ($mut, $end1, $start1, $pos, $chr) = ("")x5;
my (%merge_info, %mutation, %vcf_info, %mutfa);
open FF,"$id.sort.vcf" || die $!;
while(<FF>){
	chomp;
	next if /^#/;
	my @terms = split /\t/,$_;
	my $len	  = length $terms[3];
	my $end	  = $terms[1]+$len-1;
	my @info  = split /:/,$terms[8];
	my %info;
	foreach my $i(0..$#info){$info{$info[$i]}=$i}
	my @mutinfo = split /:/,$terms[9];
	my $af;
	if($info{AF}){
         $af = $mutinfo[$info{AF}];
    }else{
        $mutinfo[$info{FREQ}] =~ s/%//;  $af = $mutinfo[$info{FREQ}];
    }
	$vcf_info{"$terms[0].$terms[1].$terms[3].$terms[4]"}++;

	if($terms[0] eq $chr and ($terms[1]-$end1<$distance)){
		next if $terms[1] == $start1;
		$mutation{$n}{mutline}{$_}++;
		$mutation{$n}{chr}=$terms[0];
		if(not $mutation{$n}{start}){
			my @mut=split /\t/,$mut;
			$mutation{$n}{start}=$mut[1];
		    my $len=length $mut[3];
			my $end1=$mut[1]+$len-1;
			$mutation{$n}{end}=$end1;
		}
		if($end> $mutation{$n}{end}){$mutation{$n}{end}=$end}
        
		if($mutation{$n}{AF}){
			if($mutation{$n}{AF} > $af){
				$mutation{$n}{AF} = $af;
				$mutation{$n}{info} = $_;
				$merge_info{"$terms[0].$mutation{$n}{start}"}=$_;
			}
		}else{
			$mutation{$n}{AF} = $af;
			$mutation{$n}{info} = $_;
			$merge_info{"$terms[0].$mutation{$n}{start}"}=$_;
		}
	}else{
		$n++;
		$mutation{$n}{AF}=$af;
		$mutation{$n}{info}=$_;
		$merge_info{"$terms[0].$terms[1]"}=$_;
	}
	$chr   = $terms[0];
	$pos   = $terms[1];
	$start1= $terms[1];
	$end1  = $terms[1]+$len-1;
	$mut   = $_;
	$mutfa{"$terms[0].$terms[1]"}{ref}  = $terms[3];
	$mutfa{"$terms[0].$terms[1]"}{alle} = $terms[4];
}
close FF;


open FOV,">$id.new.vcf"  || die $!;
open FOOV,">$id.ori.vcf" || die $!;
open FOI,">$id.mutinfo.xls" || die $!;
#open FOTS,">$id.tview.sh" || die $!;

foreach my $o(sort keys %mutation){
	print "Now $mutation{$o}{chr}:$mutation{$o}{start}-$mutation{$o}{end}\n";
	my $pinf="$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}";

	#########################################################################################
	open FA,"$samtools faidx $hg19 $mutation{$o}{chr}:$mutation{$o}{start}-$mutation{$o}{end} |" || die $!;
	#>chr8:31497795-31497805
    #ggaaggcggc
    my ($fasta, %fasta);
	while(<FA>){
		chomp;
		if(not /^>/){
			$fasta = $_;
			$fasta = uc $fasta;
			$_	   = uc $_;
			my @seq= split '',$_;
			my $p  = $mutation{$o}{start};
			foreach my $b(@seq){$fasta{$p}=$b; $p++;}
		}
	}
	close(FA);

	############################################################################################
	open FF, "$samtools view -h $bam $mutation{$o}{chr}:$mutation{$o}{start}-$mutation{$o}{end} | ";
	#open  FBA, ">$ARGV[2].$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}.bed.sam";
	my $read_num=0;
	my %mutreads;
	my $total_reads=0;
	while(<FF>){
		next if(/^@/);
		chomp;
		$total_reads++;
        next if /MD:Z:\d+\s+/;
		my @terms=split /\t/,$_;
		next if $terms[4] < 1;
		next if $terms[1] & 0x400;
		my $dis=0;
		if(/NM:i:(\d+)/){$dis=$1;}
	    next if /NM:i:1\s+/;
		next if /NM:i:0\s+/;
			
		$terms[9]=uc $terms[9];
		my $mutp=$terms[2];
		my $mutn=0;
		my %alle;
		my $alleseq;
		my $alignseq=$terms[9];
		my $seqcp=$terms[9];
		my $oriseq=$terms[9];
		#96M3I34M13H
		
		if($terms[5]=~m/^(\d+)S/){$seqcp=substr($terms[9],$1,)}
		my $pp=$terms[3];
		my $dp=1;
		my $deseq;
		while($terms[5]=~m/(\d+)([A-Z])/g){
			my $d=$1;
			my $m=$2;
			if($m eq "S" or $m eq  "H"){
				next;
			}elsif($m eq "M" ){
				$deseq.=substr($seqcp,$dp-1,$d);
				$pp+=$d-1;
				$dp+=$d;
			}elsif($m eq "I" ){
				my $seq=substr($seqcp,$dp-2,$d+1);
				$mutp.=".$pp";
				$alle{$pp}=$seq;
				$dp+=$d;
				$mutn++;
			}elsif($m eq "D" ){
				$deseq.=substr($seqcp,$dp-1,$d);
				$pp+=$d-1;
				$mutn++;
			}
		}
		
		my $align="";
		my $l1=length $terms[9];
		my $l2=length $deseq;
		$terms[9]=$deseq;
		if(/MD:Z:(\S+)/){
			$align=$1;
			$align=uc $align;
			my $s=$terms[3];
			my $loc=1;
			#11G2^CA136
			#124^AATTAAGAGA6^ATCTCCGA21
			#print "3\t$terms[0]\t$terms[5]\t$align\n";
			while($align=~m/(\d+)([A-Z]|\^[A-Z]+)/g){
				$loc+=$1;
				my $base=$2;					
				my $ml=$1;					
				$s+=$1;
				my $p=$s;
				$mutp.=".$p";
				$base=uc $base;
				my $len=length $base;
				if($base=~m/\^/){
					$base=~s/\^//;
					my $len=length $base;
					my $i=1;
					while($i<=$len){
						$alle{$p}="^";
						$alleseq.="^";
						$i++;
						$p++;
					}
				    $s+=$len;
					$loc+=$len;
				}else{
					my $i=0;
					my $seq=substr($terms[9],$loc-1,$len);
					my @seq=split '',$seq;
					while($i<$len){
						my $b=$seq[$i];
						$alle{$p}=$b;
						$alleseq.=$b;
						$mutn++ unless $b=~m/^N+$/;
						#print "4\t$terms[0] $terms[5] $oriseq $terms[9] $align  $p $alle{$p}\n";
						$i++;
						$p++;
					}
					$s+=$len;
					$loc+=$len;
					#print "4\t$terms[0].$terms[1].$terms[3].$terms[5] $terms[9] $align $base $p $alle{$p}\n";
				}
			}
		}
		if ($mutn>1){
			#print FBA "$_\n";
			$mutreads{$pinf}{$mutp}{$alleseq}{total}++;
			$mutreads{$pinf}{$mutp}{$alleseq}{info}.="$terms[0].$terms[3].$terms[5].$align ,";
			foreach my $a(keys %alle){$mutreads{$pinf}{$mutp}{$alleseq}{$a}=$alle{$a}}
		}
	}
	#close FBA;

	#system("$samtools sort -o $ARGV[2].$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}.bed.bam $ARGV[2].$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}.bed.sam");
    #system("$samtools index   $ARGV[2].$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}.bed.bam ");
    #print FOTS "$samtools tview -p   $mutation{$o}{chr}:$mutation{$o}{start} $ARGV[2].$mutation{$o}{chr}.$mutation{$o}{start}.$mutation{$o}{end}.bed.bam $hg19\n";
	
	my %result;
	my %ouputinfo;
	foreach my $p(sort keys %mutreads){
		foreach my $subp(sort keys %{$mutreads{$p}}){
			foreach my $seq(sort keys %{$mutreads{$p}{$subp}}){
			    if ($mutreads{$p}{$subp}{$seq}{total}<6){
					print FOOV "$p\t$subp\t$seq\t$mutreads{$p}{$subp}{$seq}{total}\n"; next;
				}
                my @inf=split /\./,$p;
				my $alle;
				foreach my $l($inf[1]..$inf[2]){
					if($mutreads{$p}{$subp}{$seq}{$l}){
						$alle.=$mutreads{$p}{$subp}{$seq}{$l};
					}else{
						$alle.=$fasta{$l};
					}
					my $p="$inf[0].$l";
					$ouputinfo{$p}=$merge_info{$p} if $merge_info{$p};
				}
				$alle=~s/\^//g;
				my $fre=$mutreads{$p}{$subp}{$seq}{total}/$total_reads;
				my $aligninfo=$mutreads{$p}{$subp}{$seq}{info};
				$aligninfo=~s/\s+//g;
				my @aligninfo=split /,/,$aligninfo;
				my %dupinf;
				foreach my $inf(@aligninfo){my @tmp=split /\./,$inf;$dupinf{$tmp[1].$tmp[2].$tmp[3]}++;}

				my $dupr=keys %dupinf;
				print FOI "$inf[0]\t$inf[1]\t$p,$subp,$seq\t$fasta\t$alle\t$mutreads{$p}{$subp}{$seq}{total}\t$mutreads{$p}{$subp}{$seq}{total}\t$fre\t$dupr\t$mutreads{$p}{$subp}{$seq}{info}\n";
                if($fasta ne $alle){
					$result{"$inf[0]\t$inf[1]\t.\t$fasta\t$alle"}{total}+=$dupr; 
				}
			}
		}
	}
	foreach my $m(sort keys %result){
		  my $fre=$result{$m}{total}/$total_reads;
		  my $normal_ref_dp=$total_bases-$mismatch_base;
          my $normal_total_dp=$total_bases;
          my $tumor_ref_dp=$total_reads-$result{$m}{total};
          my $tumor_total_dp=$total_reads;
         
         my $normal_twotailed_value = calculateStatistic(n11=>$tumor_ref_dp,
                                                     n1p=>$tumor_total_dp,
                                                     np1=>$tumor_ref_dp+ $normal_ref_dp,
                                                     npp=>$tumor_total_dp+$normal_total_dp);

		my @tmp=split /\t/,$m;
		if ($vcf_info{"$tmp[0].$tmp[1].$tmp[3].$tmp[4]"}){
			 print FOOV "$tmp[0].$tmp[1].$tmp[3].$tmp[4]\texit\n"; next; 
		}
		
		my $alignreads=`$samtools view $bam $mutation{$o}{chr}:$mutation{$o}{start}-$mutation{$o}{end} |awk '{print \$10}' | grep $tmp[4] | wc -l`;
		chomp $alignreads;
		if ($alignreads<5 or $result{$m}{total}<5){
			print FOOV "$m\talignreads-$alignreads\ttotal-$result{$m}{total}\n"; next;
		}
		
		if(length $tmp[3] == length $tmp[4]){
			my @seq1 = split '',$tmp[3];
			my @seq2 = split '',$tmp[4];
			my $diff = 0;
			foreach my $i(0..$#seq1){$diff++ if $seq1[$i] ne $seq2[$i]}
			if ($diff<2){print FOOV "$m\tdiff less 2\n";next}
		}
		
		my $outputinfo = $ouputinfo{"$tmp[0].$tmp[1]"};
		my @outputinfo = split /\t/,$outputinfo;
		
		print FOV "$m\t$outputinfo[5]\t$outputinfo[6]\t$outputinfo[7];Merge=Yes;MAF=$fre;MDP=$result{$m}{total};MPV=$normal_twotailed_value\t$outputinfo[8]\t$outputinfo[9]\t$outputinfo[10]\n";
	}
}
system("cat $vcf $id.new.vcf >$id.merged.vcf");


sub usage{
	print STDERR << "EOF";
	usage:[-vcf -bam -id -d]
		-h      :help message
		-vcf    :merged vcf file
		-bam    :recalbam file
		-id     :sample id
		-d	:distance of two mutation(default 5)

	Example:
		perl $0 -vcf Merged.vcf -bam recal.bam -id BJ22CM000993 -d 5

EOF
exit(1);
}
