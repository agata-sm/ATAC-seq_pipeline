#!c:/perl/bin/perl.exe


# script to generate gtf with coordinates of annotated peaks
# proj 4833
# 27 - 28 i 2020

# Agata Smialowska


use warnings;
use strict;
use diagnostics;
use Getopt::Long;


my $script_name="peaks2gff.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/annotated_peaks.bed\n";
	print "--outfile_gtf: /path/to/gtf\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile_gtf=s'		=>	\(my $outfile_gtf)
	) or die "Error in command line arguments";


	#outfile info contains more info on genes annotated to peaks
	my $outfile_info="$outfile_gtf\_genes_info.tab"; 



	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	open (OUTFILE_GTF, ">","$outfile_gtf") or die "Cannot open input file $outfile_gtf: $!"; 
	open (OUTFILE_INFO, ">","$outfile_info") or die "Cannot open input file $outfile_info: $!"; 

	my $header_info="peakid\tchr\tstart\tend\tensembl_gene_5\tname_gene_5\tbiotype_gene_5\tdistance_gene_5\tensembl_gene_3\tname_gene_3\tdistance_gene_3\tbiotype_gene_3";
	print OUTFILE_INFO "$header_info\n";


	my $peak_number=1;

	while(<INFILE>){
		chomp $_;
		my @line=split/\|/,$_;
		my $dist5=$line[2];
		my $gene5=$line[1];
		my $dist3=$line[4];
		my $gene3=$line[3];

		$gene5=~s/\"//g;
		$gene3=~s/\"//g;


		my @peak_bed=split/\s/,$line[0];
		my $start=$peak_bed[1]+1;

		my $peak_name="peak_$peak_number";
		$peak_number+=1;

		my $gtf_line1="$peak_bed[0]\tmacs2\tatacseq_peak\t$start\t$peak_bed[2]\t.\t.\t.";
		

		my $dist;
		my $gene_id;
		my $gene_name;
		my $gene_biotype;

		my $gene5_id;
		my $gene5_name;
		my $gene5_biotype;

		my $gene3_id;
		my $gene3_name;
		my $gene3_biotype;

		if ( ($dist5 ne qw /NA/) && ($dist3 ne qw /NA/) ) {
			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr)=split/\t/,$gene5;
			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr)=split/\t/,$gene3;

			my @gene5_info=split/;\s/,$gene5_attr;
			my @gene3_info=split/;\s/,$gene3_attr;


			$gene5_id=$gene5_info[0];
			$gene5_id=~s/gene_id\s+//;
			$gene3_id=$gene3_info[0];
			$gene3_id=~s/gene_id\s+//;


			$gene5_biotype=$gene5_info[4];
			$gene5_biotype=~s/gene_biotype\s+//;
			$gene3_biotype=$gene3_info[4];
			$gene3_biotype=~s/gene_biotype\s+//;
			$gene3_biotype=~s/;//g;
			$gene5_biotype=~s/;//g;


			$gene5_name=$gene5_info[2];
			$gene5_name=~s/gene_name\s+//;
			$gene3_name=$gene3_info[2];
			$gene3_name=~s/gene_name\s+//;


			if(abs($dist5)<abs($dist3)){
				$dist=$dist5;
				$gene_id=$gene5_id;
				$gene_name=$gene5_name;
				$gene_biotype=$gene5_biotype;


			}elsif (abs($dist3)<abs($dist5)) {
				$dist=$dist3;
				$gene_id=$gene3_id;
				$gene_name=$gene3_name;
				$gene_biotype=$gene3_biotype;

			}elsif (abs($dist5)==abs($dist3)){
				$dist=$dist5;


				$gene_id="$gene5_id\_$gene3_id";

				$gene_biotype="$gene5_biotype\_$gene3_biotype";

				$gene_name="$gene5_name\_$gene3_name";

			} 
		}elsif( ($dist5 eq qw /NA/) && ($dist3 ne qw /NA/) ){

			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr)=split/\t/,$gene3;
			my @gene3_info=split/;\s/,$gene3_attr;
			$gene3_id=$gene3_info[0];
			$gene3_id=~s/gene_id\s+//;
			$gene3_biotype=$gene3_info[4];
			$gene3_biotype=~s/gene_biotype\s+//;
			$gene3_biotype=~s/;//g;
			$gene3_name=$gene3_info[2];
			$gene3_name=~s/gene_name\s+//;
	

			$dist=$dist3;
			$gene_id=$gene3_id;
			$gene_name=$gene3_name;
			$gene_biotype=$gene3_biotype;

			$gene5_id="NA";
			$gene5_name="NA";
			$gene5_biotype="NA";

		}elsif( ($dist5 ne qw /NA/) && ($dist3 eq qw /NA/) ){

			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr)=split/\t/,$gene5;
			my @gene5_info=split/;\s/,$gene5_attr;
			$gene5_id=$gene5_info[0];
			$gene5_id=~s/gene_id\s+//;
			$gene5_biotype=$gene5_info[4];
			$gene5_biotype=~s/gene_biotype\s+//;
			$gene5_biotype=~s/;//g;
			$gene5_name=$gene5_info[2];
			$gene5_name=~s/gene_name\s+//;

			$dist=$dist5;
			$gene_id=$gene5_id;
			$gene_name=$gene5_name;
			$gene_biotype=$gene5_biotype;

			$gene3_id="NA";
			$gene3_name="NA";
			$gene3_biotype="NA";


		
		}

		my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";


		print OUTFILE_GTF "$gtf_line1\t$gtf_line2\n";

		print OUTFILE_INFO "$peak_name\t$peak_bed[0]\t$start\t$peak_bed[2]\t$gene5_id\t$gene5_name\t$dist5\t$gene5_biotype\t$gene3_id\t$gene3_name\t$dist3\t$gene3_biotype\n";

	}

}



