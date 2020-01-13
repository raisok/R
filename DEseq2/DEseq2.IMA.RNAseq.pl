#!/usr/bin/perl


use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;
use threads;

my ($Group, $Diff, $Count, $Log2FC, $Padj, $Desc, $Rscript, $Outdir);
GetOptions(
	"group:s" => \$Group,
	"diff:s" => \$Diff,
	"count:s" => \$Count,
	"log2:s" => \$Log2FC,
	"padj:s" => \$Padj,
	"desc:s" => \$Desc,
	"Rpath:s" => \$Rscript,
	"outdir:s" => \$Outdir
);
if (!$Group || !$Diff || !$Count || !$Desc) {
	print STDERR <<USAGE;
Description: DESeq2 for Gene Differential Expression Analysis
Usage: perl $0 [options]
Options:
	* -group      file of samples' grouping, required while repeat comparison occured
                        format:
                           GroupA           SampleA,SampleC
                           GroupB           SampleB,SampleD
	* -diff       file of comparison plan
			format<Control Treat>:
			   SampleA,SampleC	    SampleB,SampleD	    # <for repeat comparison>
	* -count      file of multiple samples reads count
			format:
			Gene	SampleA		SampleB	SampleC	SampleD
			geneID		5		7	4	6
			* note the reads of Sample has three repeat,read num can't be zero
	  -log2       cutoff for log2 fold changes, default: 1
	  -padj       cutoff for adjust-pvalue, default: 0.05
	* -desc       gene description file, header required, format:GeneID<tab>Desc1
	  -Rpath      path of Rscript, default: $Bin/Rscript
	  -outdir     directory of output files, default: current directory
E.g.:
	perl $0 -group group.list -diff diff.list -count count.list -desc description.txt

USAGE
	exit;
}

$Log2FC  ||= 1;
$Padj    ||= 0.05;
$Rscript ||= "$Bin/Rscript";
$Outdir  ||= "./";

my %group = ();
my %desc  = ();
my %list  = ();

if ($Group) {
	open my $fh_group, "<", $Group or die $!;
	while (<$fh_group>) {
		chomp;
		my ($group, $samples) = split; 
		$group{$samples} = $group;	# $group{"HBRR1,HBRR2,HBRR3"} = HBRR
	}
	close $fh_group;
}

if ($Desc) {
	open my $fh_desc, "<", $Desc or die $!;
	chomp(my$desc_header = <$fh_desc>);
	while (<$fh_desc>) {
		chomp;
		my @a = split /\s+/,$_,2;
		$desc{$a[0]} = $a[1];
	}
	close $fh_desc;
}


my @threads = ();
open my $fh_diff, "<", $Diff or die $!;		#CompareList.txt: HBRR1,HBRR2,HBRR3	UHRR1,UHRR2,UHRR3
while (my $line = <$fh_diff>) {
	my $th = async {
		chomp $line;
		my ($control_name, $treat_name) = split /\s+/, $line; #HBRR1,HBRR2,HBRR3  UHRR1,UHRR2,UHRR3
		my @control = split /,/, $control_name; #HBRR1,HBRR2,HBRR3
		my @treat   = split /,/, $treat_name;   #UHRR1,UHRR2,UHRR3
		if (@control > 1){($group{$control_name}) ? $control_name = $group{$control_name} : die "Can't find $control_name in group.list\n";}   # $control_name  = HBRR
		if (@treat > 1){($group{$treat_name}) ? $treat_name = $group{$treat_name} : die "Can't find $treat_name in group.list\n";}             # $treat_name = UHRR
		my @factors = (map("Control",(1..@control)),map("Treat",(1..@treat))); # Control Control Control Treat Treat Treat
		my @names   = (@control,@treat); # HBRR1,HBRR2,HBRR3,UHRR1,UHRR2,UHRR3
	
		my $output_a = "$Outdir/$control_name-VS-$treat_name.deseq2.total_genes";
		my $output_b = "$Outdir/$control_name-VS-$treat_name.differentially_expressed_genes";
		#my $output_a = "$Outdir/$control_name_vs_$treat_name.total_genes.xls";
		#my $output_b = "$Outdir/$control_name_vs_$treat_name.differentially_expressed_genes.xls";
		my $rcode = "$Outdir/$control_name-VS-$treat_name.deseq2.R";
		
		
		my $type = join(",",map("\"$_\"",@factors)); # Control,Control,Control,Treat,Treat,Treat
		my $samples = join(",",map("\"$_\"",@names)); # HBRR1,HBRR2,HBRR3,UHRR1,UHRR2,UHRR3
		open my $fh_rcode, ">", $rcode or die $!;
		print $fh_rcode <<RCODE;
library("DESeq2")

data <- read.csv("$Count",header=T,sep=",",check.names=F)
countData <- subset(data, select = c($samples))

keep <- apply(countData,1,sum) >= 10
countData <- countData[keep,]
keepData <- data[keep,]

condition <- c($type)
coldata <- data.frame(condition)
rownames(coldata) <- c($samples)

dds <- DESeqDataSetFromMatrix(countData=countData,colData=coldata,design=~condition)
dds2 <- DESeq(dds)
res <- results(dds2,contrast=c("condition","Treat","Control"))

res2 <- cbind(keepData[,1],apply(countData,1,sum),countData,counts(dds2,normalized=TRUE),as.matrix(res))
colnames(res2)[1:2] <- c("gene_id","total_count")
rename <- paste(colnames(res2)[(2+length(condition)+1):(2+length(condition)*2)],"(norm)",sep="")
colnames(res2)[(2+length(condition)+1):(2+length(condition)*2)] <- rename

resOrdered <- res2[order(res2\$padj),]
write.table(resOrdered,file="$output_a",sep="\\t",quote=F,row.names=F)

DEGs_FC <- abs(resOrdered\$log2FoldChange) >= "$Log2FC" & resOrdered\$padj < "$Padj" & !is.na(resOrdered\$padj)
write.table(resOrdered[DEGs_FC,],file="$output_b",sep="\\t",quote=F,row.names=F)

RCODE
		system("$Rscript $rcode");

		die "Can't find DESeq2 output file, please check...\n" unless (-f $output_a && -f $output_b);
		open my $fh_diffexp, ">", "$Outdir/$control_name\_vs\_$treat_name.total_genes.xls" or die $!;
		open my $fh_diffexpfilter, ">", "$Outdir/$control_name\_vs\_$treat_name.differentially_expressed_genes.xls" or die $!;
		
		open my $fh_deg, "<", $output_b or die $!;
		my $degheader=<$fh_deg>;
		chomp($degheader);
		my $degtitle=(split/\s+/,$degheader,2)[1];
		print $fh_diffexpfilter "gene_id\tdescription\t$degtitle\n";
		while (<$fh_deg>) {
			chomp;
			my @a = split /\s+/,$_,2;
			if (exists $desc{$a[0]}){
				print $fh_diffexpfilter "$a[0]\t$desc{$a[0]}\t$a[1]\n";
			}else{
				print $fh_diffexpfilter "$a[0]\tNA\t$a[1]\n";
			}
		}
		
		open my $fh_totalgene, "<", $output_a or die $!;
		my $totalgeneheader=<$fh_totalgene>;
		chomp($totalgeneheader);
		my $totaltitle=(split/\s+/,$totalgeneheader,2)[1];
		print $fh_diffexp "gene_id\tdescription\t$totaltitle\n";
		while (<$fh_totalgene>) {
			chomp;my @b = split /\s+/,$_,2;
			if (exists $desc{$b[0]}){
				print $fh_diffexp "$b[0]\t$desc{$b[0]}\t$b[1]\n";
			}else{
				print $fh_diffexp "$b[0]\tNA\t$b[1]\n";
			}
		}
		close $fh_deg;
		close $fh_totalgene;
		close $fh_diffexp;
		close $fh_diffexpfilter;
		return 0;
	};
	push @threads, $th;
}

sub flt_to_pct 
{
        sprintf( "\%d", shift );
}

foreach my $th (@threads) {
	die "Thread Error\n" unless ($th -> join() == 0);
}
exit;
