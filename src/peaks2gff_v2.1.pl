#!c:/perl/bin/perl.exe


# script to generate gtf with coordinates of annotated peaks
# additional output
# bed and gtf files of genes with peaks within select distance (2kb is recommended)
# bed files of genomic coordinates for genes selected by peak distance, padded by desired value (only on TSS end)

# proj 4833
# 27 - 28 i 2020

# author: Agata Smialowska


use warnings;
use strict;
use diagnostics;
use Getopt::Long;


my $script_name="peaks2gff_v2.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/annotated_peaks.bed\n";
	print "--outfile_gtf: /path/to/gtf\n";
	print "--peak_distance: cutoff distance gene to the annotated peak\n";
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
		print "--peak_distance: cutoff distance gene to the annotated peak\n";
		print "--gene_coord_distance: distance added to genomic coordinates on gene 5'end (upstream of TSS) for genes selected by --peak_distance\n";
		exit;
	}

	#outfile info contains more info on genes annotated to peaks
	my $outfile_info="$outfile_gtf_base\_genes_info.tab"; 
	
	#outfile_gtf_CO_peaks coordinates of PEAKS within 2kh of genes
	my $outfile_gtf_CO_peaks="$outfile_gtf_base\.peaks_$peak_dist\.gtf";

	#outfiles with gene coordinates
	my $outfile_bed_genes="$outfile_gtf_base\.genes_$peak_dist\.bed";
	my $outfile_gtf_genes="$outfile_gtf_base\.genes_$peak_dist\.gtf";

	#outfile with gene coordinates padded (gene_coord_distance) on the TSS side (5´of the gene in + and 3' of the gene on -)
	my $outfile_bed_padded_genes="$outfile_gtf_base\.genes_$peak_dist\.gene_coordinates_$padded_dist\.bed";

	open (INFILE, "<","$infile") or die "Cannot open input file $infile: $!"; 

	open (OUTFILE_GTF_PEAKS, ">","$outfile_gtf_all_peaks") or die "Cannot open output file $outfile_gtf_all_peaks: $!"; 
	open (OUTFILE_GTF_PEAKS_CO, ">","$outfile_gtf_CO_peaks") or die "Cannot open output file $outfile_gtf_CO_peaks: $!"; 
	open (OUTFILE_INFO, ">","$outfile_info") or die "Cannot open output file $outfile_info: $!"; 
	my $header_info="peakid\tchr\tstart\tend\tensembl_gene_5\tname_gene_5\tbiotype_gene_5\tdistance_gene_5\tensembl_gene_3\tname_gene_3\tdistance_gene_3\tbiotype_gene_3";
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
			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr,$gene5_score,$gene5_strand)=split/\t/,$gene5;
			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr,$gene3_score,$gene3_strand)=split/\t/,$gene3;

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
				my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";

				if ($dist5<$peak_dist){
					print OUTFILE_BED_GENES "$gene5\n";
					my $gene_start_gtf=$gene5_start+1;
					my $gtf_line_gene="$gene5_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene5_end\t.\t$gene5_strand\t.\t$gene5_attr";
					print OUTFILE_GTF_GENES "$gtf_line_gene\n";

					#padded gene coordinates
					my @gene_bed=split/\t/,$gene5;
					if($gene_bed[5] eq qw /+/){
						$gene_bed[1]=$gene_bed[1]-$padded_dist;
					}
					elsif($gene_bed[5] eq qw /-/){
						$gene_bed[2]=$gene_bed[2]+$padded_dist;
					}

					my $line_padded_genes_bed=join("\t",@gene_bed);
					print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";

					print OUTFILE_GTF_PEAKS_CO "$gtf_line1\t$gtf_line2\n";

				}

			}elsif (abs($dist3)<abs($dist5)) {
				$dist=$dist3;
				$gene_id=$gene3_id;
				$gene_name=$gene3_name;
				$gene_biotype=$gene3_biotype;
				my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";

				if ($dist3<$peak_dist){
					print OUTFILE_BED_GENES "$gene3\n";
					my $gene_start_gtf=$gene3_start+1;
					my $gtf_line_gene="$gene3_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene3_end\t.\t$gene3_strand\t.\t$gene3_attr";
					print OUTFILE_GTF_GENES "$gtf_line_gene\n";

					#padded gene coordinates
					my @gene_bed=split/\t/,$gene3;
					if($gene_bed[5] eq qw /+/){
						$gene_bed[1]=$gene_bed[1]-$padded_dist;
					}
					elsif($gene_bed[5] eq qw /-/){
						$gene_bed[2]=$gene_bed[2]+$padded_dist;
					}

					my $line_padded_genes_bed=join("\t",@gene_bed);
					print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";

					print OUTFILE_GTF_PEAKS_CO "$gtf_line1\t$gtf_line2\n";

				}



			}elsif (abs($dist5)==abs($dist3)){
				$dist=$dist5;


				$gene_id="$gene5_id\_$gene3_id";

				$gene_biotype="$gene5_biotype\_$gene3_biotype";

				$gene_name="$gene5_name\_$gene3_name";

				if ($dist3<$peak_dist){
					print OUTFILE_BED_GENES "$gene3\n";
					print OUTFILE_BED_GENES "$gene5\n";
					my $gene_start_gtf=$gene3_start+1;
					my $gtf_line_gene1="$gene3_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene3_end\t.\t$gene3_strand\t.\t$gene3_attr";
					print OUTFILE_GTF_GENES "$gtf_line_gene1\n";

					my $gene_start_gtf2=$gene5_start+1;
					my $gtf_line_gene2="$gene5_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene5_end\t.\t$gene5_strand\t.\t$gene5_attr";
					print OUTFILE_GTF_GENES "$gtf_line_gene2\n";


					my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";
					print OUTFILE_GTF_PEAKS_CO "$gtf_line1\t$gtf_line2\n";


					#padded gene coordinates
					my @gene_bed=split/\t/,$gene3;
					if($gene_bed[5] eq qw /+/){
						$gene_bed[1]=$gene_bed[1]-$padded_dist;
					}
					elsif($gene_bed[5] eq qw /-/){
						$gene_bed[2]=$gene_bed[2]+$padded_dist;
					}

					my $line_padded_genes_bed=join("\t",@gene_bed);
					print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";

					#padded gene coordinates
					my @gene_bed_5=split/\t/,$gene5;
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

		}elsif( ($dist5 eq qw /NA/) && ($dist3 ne qw /NA/) ){

			my ($gene3_chr,$gene3_start,$gene3_end,$gene3_attr,$gene3_score,$gene3_strand)=split/\t/,$gene3;
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


			if ($dist<$peak_dist){
				print OUTFILE_BED_GENES "$gene3\n";
				my $gene_start_gtf=$gene3_start+1;
				my $gtf_line_gene="$gene3_chr\tensembl\tgene_atacseq_peak_2kb\t$gene_start_gtf\t$gene3_end\t.\t$gene3_strand\t.\t$gene3_attr";
				print OUTFILE_GTF_GENES "$gtf_line_gene\n";

				#padded gene coordinates
				my @gene_bed=split/\t/,$gene3;
				if($gene_bed[5] eq qw /+/){
					$gene_bed[1]=$gene_bed[1]-$padded_dist;
				}
				elsif($gene_bed[5] eq qw /-/){
					$gene_bed[2]=$gene_bed[2]+$padded_dist;
				}

				my $line_padded_genes_bed=join("\t",@gene_bed);
				print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";

				my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";
				print OUTFILE_GTF_PEAKS_CO "$gtf_line1\t$gtf_line2\n";

			}


		}elsif( ($dist5 ne qw /NA/) && ($dist3 eq qw /NA/) ){

			my ($gene5_chr,$gene5_start,$gene5_end,$gene5_attr,$gene5_score,$gene5_strand)=split/\t/,$gene5;
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


			if ($dist<$peak_dist){
				print OUTFILE_BED_GENES "$gene5\n";
				my $gene_start_gtf=$gene5_start+1;
				my $gtf_line_gene="$gene5_chr\tensembl\tgene_atacseq_peak_2kb_atacseq_peak\t$gene_start_gtf\t$gene5_end\t.\t$gene5_strand\t.\t$gene5_attr";
				print OUTFILE_GTF_GENES "$gtf_line_gene\n";

				#padded gene coordinates
				my @gene_bed=split/\t/,$gene5;
				if($gene_bed[5] eq qw /+/){
					$gene_bed[1]=$gene_bed[1]-$padded_dist;
				}
				elsif($gene_bed[5] eq qw /-/){
					$gene_bed[2]=$gene_bed[2]+$padded_dist;
				}

				my $line_padded_genes_bed=join("\t",@gene_bed);
				print OUTFILE_BED_PADDED_GENES "$line_padded_genes_bed\n";

				my $gtf_line2="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";
				print OUTFILE_GTF_PEAKS_CO "$gtf_line1\t$gtf_line2\n";



			}

		
		}

		my $gtf_line2_b="peak_id $peak_name; gene_id $gene_id; gene_name $gene_name; gene_biotype $gene_biotype";

		print OUTFILE_GTF_PEAKS "$gtf_line1\t$gtf_line2_b\n";

		print OUTFILE_INFO "$peak_name\t$peak_bed[0]\t$start\t$peak_bed[2]\t$gene5_id\t$gene5_name\t$dist5\t$gene5_biotype\t$gene3_id\t$gene3_name\t$dist3\t$gene3_biotype\n";




	}

}


close(INFILE);
close(OUTFILE_INFO);
close(OUTFILE_GTF_PEAKS);
close(OUTFILE_GTF_PEAKS_CO);
close(OUTFILE_GTF_GENES);
close(OUTFILE_BED_GENES);
close(OUTFILE_BED_PADDED_GENES);

exit;
