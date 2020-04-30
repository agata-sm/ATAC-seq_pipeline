#!c:/perl/bin/perl.exe


# script to generate gtf with coordinates of annotated peaks (output of BEDOPS closest-feature)
# change from previous versions: v3 selects genes where peaks are in the TSS region 
# (i.e. genes on (-) strand on the 5' peak side and  genes on (+) strand on the 3' side)
# additional output
# bed and gtf files of genes with peaks within select distance from their TSS (2kb is recommended)
# bed files of genomic coordinates for genes selected by peak distance, padded by desired value (only on TSS end)

# TSS_dist is calculate relative to TSS with gene strand taken into account! - is upstream (i.e. outside of gene)  and + is downstream (i.e. within gene)

# proj 4833
# 29 iv 2020

# author: Agata Smialowska


use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use POSIX;

my $script_name="peaks2gff_v3.0.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/annotated_peaks.bed\n";
	print "--outfile_gtf: /path/to/gtf\n";
	print "--peak_distance: cutoff distance gene TSS to the annotated peak\n";
	print "--gene_coord_distance: distance added to genomic coordinates on gene 5'end (upstream of TSS) for genes selected by --peak_distance\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile_gtf=s'		=>	\(my $outfile_gtf_all_peaks),
		'peak_distance=s'		=>	\(my $peak_dist),
		'gene_coord_distance=s'		=>	\(my $padded_dist)
	) or die "Error in command line arguments";

	my $outfile_gtf_base;
	if ($outfile_gtf_all_peaks=~m/(.+).gtf/){
		$outfile_gtf_base=$1;
	}else{
		print "outfile_gtf file name does not have a gtf extension. Please provide a valid name.";
		print "arguments:\n";
		print "--infile: /path/to/annotated_peaks.bed\n";
		print "--outfile_gtf: /path/to/gtf\n";
		print "--peak_distance: cutoff distance gene TSS to the annotated peak\n";
		print "--gene_coord_distance: distance added to genomic coordinates on gene 5'end (upstream of TSS) for genes selected by --peak_distance\n";
		exit;
	}

	#outfile info contains more info on genes annotated to peaks
	my $outfile_info="$outfile_gtf_base\_genes_info.tab"; 
	
	#outfile_gtf_CO_peaks coordinates of PEAKS within 2kh of genes
	my $outfile_gtf_CO_peaks="$outfile_gtf_base\.peaks\.TSS_dist_$peak_dist\.gtf";

	#outfiles with gene coordinates
	my $outfile_bed_genes="$outfile_gtf_base\.genes\.TSS_dist_$peak_dist\.bed";
	my $outfile_gtf_genes="$outfile_gtf_base\.genes\.TSS_dist_$peak_dist\.gtf";

	#outfile with gene coordinates padded (gene_coord_distance) on the TSS side (5´of the gene in + and 3' of the gene on -)
	my $outfile_bed_padded_genes="$outfile_gtf_base\.genes\.TSS_dist_$peak_dist\.gene_coordinates_$padded_dist\.bed";

	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	open (OUTFILE_GTF_PEAKS, ">","$outfile_gtf_all_peaks") or die "Cannot open output file $outfile_gtf_all_peaks: $!"; 
	open (OUTFILE_GTF_PEAKS_CO, ">","$outfile_gtf_CO_peaks") or die "Cannot open output file $outfile_gtf_CO_peaks: $!"; 
	open (OUTFILE_INFO, ">","$outfile_info") or die "Cannot open output file $outfile_info: $!"; 
	my $header_info="peakid\tchr\tstart\tend\tensembl_gene_5\tname_gene_5\tbiotype_gene_5\tdistance_gene_5\tensembl_gene_3\tname_gene_3\tbiotype_gene_3\tdistance_gene_3";
	print OUTFILE_INFO "$header_info\n";

	open (OUTFILE_BED_GENES, ">", "$outfile_bed_genes") or die "Cannot open output file $outfile_bed_genes: $!";
	open (OUTFILE_GTF_GENES, ">","$outfile_gtf_genes") or die "Cannot open output file $outfile_gtf_genes: $!";
	open (OUTFILE_BED_PADDED_GENES, ">","$outfile_bed_padded_genes") or die "Cannot open output file $outfile_bed_padded_genes: $!";

	my $peak_number=1; #give a unique ID to peaks: consecutive number

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

		my $gtf_peaks_line1="$peak_bed[0]\tmacs2\tatacseq_peak\t$start\t$peak_bed[2]\t.\t.\t.";

		my $peak_midpnt=$start+ceil(($peak_bed[2]-$start)/2);



		if ( ($dist5 ne qw /NA/) && ($dist3 ne qw /NA/) ) {
			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr,$gene5_score,$gene5_strand)=split/\t/,$gene5;
			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr,$gene3_score,$gene3_strand)=split/\t/,$gene3;
			
			my @gene5_info=split/;\s/,$gene5_attr;
			my @gene3_info=split/;\s/,$gene3_attr;

			my $gene5_id=$gene5_info[0];
			$gene5_id=~s/gene_id\s+//;
			my $gene3_id=$gene3_info[0];
			$gene3_id=~s/gene_id\s+//;


			my $gene5_biotype=$gene5_info[4];
			$gene5_biotype=~s/gene_biotype\s+//;
			my $gene3_biotype=$gene3_info[4];
			$gene3_biotype=~s/gene_biotype\s+//;
			$gene3_biotype=~s/;//g;
			$gene5_biotype=~s/;//g;


			my $gene5_name=$gene5_info[2];
			$gene5_name=~s/gene_name\s+//;
			my $gene3_name=$gene3_info[2];
			$gene3_name=~s/gene_name\s+//;


			print OUTFILE_INFO "$peak_name\t$peak_bed[0]\t$start\t$peak_bed[2]\t$gene5_id\t$gene5_name\t$gene5_biotype\t$dist5\t$gene3_id\t$gene3_name\t$gene3_biotype\t$dist3\n";


		# if genes on (-) strand 5' of the peak and genes on (+) strand 3' of the peak > only these will be evaluated
			my $TSS_dist_5=$gene5_end-$peak_midpnt;
			my $TSS_dist_3=$peak_midpnt-$gene3_start;

			# select genes on (-) strand on the 5' peak side and  genes on (+) strand on the 3' side to evaluate peaks close to TSS rather than gene end
			if ( ($gene5_strand eq qw /-/) && ($gene3_strand eq qw /+/) ) { #both genes in correct orientation to be a candidate #7317

				if(abs($TSS_dist_5)<abs($TSS_dist_3)){
					&get_gene($gene5,$dist5,$TSS_dist_5 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
				}
				elsif(abs($TSS_dist_5)>abs($TSS_dist_3)){
					&get_gene($gene3,$dist3,$TSS_dist_3 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
				}
				elsif(abs($TSS_dist_5)==abs($TSS_dist_3)){
					&get_2_genes($gene5,$dist5,$TSS_dist_5,$gene3,$dist3,$TSS_dist_3 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
				}



			} #if ( ($gene5_strand eq qw /-/) && ($gene3_strand eq qw /+/) ) { 
			else{
				if ($gene5_strand eq qw /-/){
					&get_gene($gene5,$dist5,$TSS_dist_5 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
				}
				elsif($gene3_strand eq qw /+/){
					&get_gene($gene3,$dist3,$TSS_dist_3 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
				}
			}
		} #if ( ($dist5 ne qw /NA/) && ($dist3 ne qw /NA/) ) {
		elsif( ($dist5 eq qw /NA/) && ($dist3 ne qw /NA/) ){

			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr,$gene3_score,$gene3_strand)=split/\t/,$gene3;
			
			my @gene3_info=split/;\s/,$gene3_attr;

			my $gene3_id=$gene3_info[0];
			$gene3_id=~s/gene_id\s+//;

			my $gene3_biotype=$gene3_info[4];
			$gene3_biotype=~s/gene_biotype\s+//;
			$gene3_biotype=~s/;//g;

			my $gene3_name=$gene3_info[2];
			$gene3_name=~s/gene_name\s+//;

			my $gene5_id="NA";
			my $gene5_name="NA";
			my $gene5_biotype="NA";

			print OUTFILE_INFO "$peak_name\t$peak_bed[0]\t$start\t$peak_bed[2]\t$gene5_id\t$gene5_name\t$gene5_biotype\t$dist5\t$gene3_id\t$gene3_name\t$gene3_biotype\t$dist3\n";
	
			my $TSS_dist_5="NA";
			my $TSS_dist_3=$peak_midpnt-$gene3_start;

			if($gene3_strand eq qw /+/){
				&get_gene($gene3,$dist3,$TSS_dist_3 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
			}

		}
		elsif( ($dist5 ne qw /NA/) && ($dist3 eq qw /NA/) ){

			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr,$gene5_score,$gene5_strand)=split/\t/,$gene5;
			
			my @gene5_info=split/;\s/,$gene5_attr;

			my $gene5_id=$gene5_info[0];
			$gene5_id=~s/gene_id\s+//;


			my $gene5_biotype=$gene5_info[4];
			$gene5_biotype=~s/gene_biotype\s+//;
			$gene5_biotype=~s/;//g;

			my $gene3_id="NA";
			my $gene3_name="NA";
			my $gene3_biotype="NA";

			my $gene5_name=$gene5_info[2];
			$gene5_name=~s/gene_name\s+//;


			print OUTFILE_INFO "$peak_name\t$peak_bed[0]\t$start\t$peak_bed[2]\t$gene5_id\t$gene5_name\t$gene5_biotype\t$dist5\t$gene3_id\t$gene3_name\t$gene3_biotype\t$dist3\n";

			my $TSS_dist_5=$gene5_end-$peak_midpnt;
			my $TSS_dist_3="NA";

			if ($gene5_strand eq qw /-/){
				&get_gene($gene5,$dist5,$TSS_dist_5 ,$gtf_peaks_line1,$peak_name,$peak_midpnt,$peak_dist,$padded_dist);
			}

		}

	}#while

}





close(INFILE);
close(OUTFILE_INFO);
close(OUTFILE_GTF_PEAKS);
close(OUTFILE_GTF_PEAKS_CO);
close(OUTFILE_GTF_GENES);
close(OUTFILE_BED_GENES);
close(OUTFILE_BED_PADDED_GENES);


sub get_gene {

	my $gene=$_[0];
	my ($gene_chr,$gene_start,$gene_end,$gene_attr,$gene_score,$gene_strand)=split/\t/,$gene;
	my @gene_info=split/;\s/,$gene_attr;

	my $gene_id=$gene_info[0];
	$gene_id=~s/gene_id\s+//;

	my $gene_biotype=$gene_info[4];
	$gene_biotype=~s/gene_biotype\s+//;
	$gene_biotype=~s/;//g;

	my $gene_name=$gene_info[2];
	$gene_name=~s/gene_name\s+//;

	my $gene_dist=$_[1];
	my $TSS_dist=$_[2];
	
	my $gtf_peaks_line1=$_[3];	
	my $peak_name=$_[4];
	my $peak_midpnt=$_[5];
	my $peak_dist=$_[6];
	my $padded_dist=$_[7];

	#my $gene_coords="$gene_start::$gene_end";
	my $gtf_peaks_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_strand $gene_strand; gene_biotype $gene_biotype; gene_dist $gene_dist; peak_midpnt: $peak_midpnt; TSS_dist $TSS_dist";

	print OUTFILE_GTF_PEAKS "$gtf_peaks_line1\t$gtf_peaks_line2\n";

  	if (abs($TSS_dist)<=$peak_dist){

  		print OUTFILE_BED_GENES "$gene\n";
  		my $gene_start_gtf=$gene_start+1;
	 	my $gtf_line_gene="$gene_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene_end\t.\t$gene_strand\t.\t$gene_attr";
	 	print OUTFILE_GTF_GENES "$gtf_line_gene\n";

	 	#padded gene coordinates
	 	my @gene_bed=split/\t/,$gene;
	 	if($gene_bed[5] eq qw /+/){
 			$gene_bed[1]=$gene_bed[1]-$padded_dist;
	 	}
	 	elsif($gene_bed[5] eq qw /-/){
	 		$gene_bed[2]=$gene_bed[2]+$padded_dist;
	 	}

		my $line_padded_genes_bed=join("\t",@gene_bed);
	 	print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";
	 	print OUTFILE_GTF_PEAKS_CO "$gtf_peaks_line1\t$gtf_peaks_line2\n";
	 }
}


sub get_2_genes {

	my $gene_5=$_[0];
	my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr,$gene5_score,$gene5_strand)=split/\t/,$gene_5;
	my @gene5_info=split/;\s/,$gene5_attr;

	my $gene5_id=$gene5_info[0];
	$gene5_id=~s/gene_id\s+//;

	my $gene5_biotype=$gene5_info[4];
	$gene5_biotype=~s/gene_biotype\s+//;
	$gene5_biotype=~s/;//g;

	my $gene5_name=$gene5_info[2];
	$gene5_name=~s/gene_name\s+//;

	my $gene5_dist=$_[1];
	my $TSS5_dist=$_[2];

	my $gene_3=$_[3];
	my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr,$gene3_score,$gene3_strand)=split/\t/,$gene_3;
	my @gene3_info=split/;\s/,$gene3_attr;

	my $gene3_id=$gene3_info[0];
	$gene3_id=~s/gene_id\s+//;

	my $gene3_biotype=$gene3_info[4];
	$gene3_biotype=~s/gene_biotype\s+//;
	$gene3_biotype=~s/;//g;

	my $gene3_name=$gene3_info[2];
	$gene3_name=~s/gene_name\s+//;

	my $gene3_dist=$_[4];
	my $TSS3_dist=$_[5];

	my $gtf_peaks_line1=$_[6];	
	my $peak_name=$_[7];
	my $peak_midpnt=$_[8];
	my $peak_dist=$_[9];
	my $padded_dist=$_[10];

	my $gene_id="$gene5_id\_$gene3_id";
	my $gene_biotype="$gene5_biotype\_$gene3_biotype";
	my $gene_name="$gene5_name\_$gene3_name";
	my $gene_strand="$gene5_strand\_$gene3_strand";
	my $gene_coords="$gene5_start::$gene5_end\_$gene3_start::$gene3_end";
	my $gene_dist="$gene5_dist\_$gene3_dist";
	my $TSS_dist=$TSS5_dist;

	my $gtf_peaks_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_coordinates $gene_coords; gene_strand $gene_strand; gene_biotype $gene_biotype; gene_dist $gene_dist; peak_midpnt: $peak_midpnt; TSS_dist $TSS_dist";

	print OUTFILE_GTF_PEAKS "$gtf_peaks_line1\t$gtf_peaks_line2\n";

  	if (abs($TSS_dist)<=$peak_dist){

 
		print OUTFILE_BED_GENES "$gene_3\n";
		print OUTFILE_BED_GENES "$gene_5\n";
		my $gene_start_gtf3=$gene3_start+1;
		my $gtf_line_gene3="$gene3_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf3\t$gene3_end\t.\t$gene3_strand\t.\t$gene3_attr";
		print OUTFILE_GTF_GENES "$gtf_line_gene3\n";

		my $gene_start_gtf5=$gene5_start+1;
		my $gtf_line_gene5="$gene5_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf5\t$gene5_end\t.\t$gene5_strand\t.\t$gene5_attr";
		print OUTFILE_GTF_GENES "$gtf_line_gene5\n";

#		my $gtf_peaks_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_coordinates $gene_coords; gene_strand $gene_strand; gene_biotype $gene_biotype; gene_dist $dist; peak_midpnt: $peak_midpnt; TSS_dist $TSS_dist";
		print OUTFILE_GTF_PEAKS_CO "$gtf_peaks_line1\t$gtf_peaks_line2\n";


		#padded gene coordinates
		my @gene_bed_3=split/\t/,$gene_3;
		if($gene_bed_3[5] eq qw /+/){
			$gene_bed_3[1]=$gene_bed_3[1]-$padded_dist;
		}
		elsif($gene_bed_3[5] eq qw /-/){
			$gene_bed_3[2]=$gene_bed_3[2]+$padded_dist;
		}

		my $line_padded_genes_bed_3=join("\t",@gene_bed_3);
		print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed_3\n";

		my @gene_bed_5=split/\t/,$gene_5;
		if($gene_bed_5[5] eq qw /+/){
			$gene_bed_5[1]=$gene_bed_5[1]-$padded_dist;
		}
		elsif($gene_bed_5[5] eq qw /-/){
			$gene_bed_5[2]=$gene_bed_5[2]+$padded_dist;
		}

		my $line_padded_genes_bed_5=join("\t",@gene_bed_5);
		print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed_5\n";
	 }
}




exit;
