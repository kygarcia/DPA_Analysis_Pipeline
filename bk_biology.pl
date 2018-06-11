#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

use Parallel::ForkManager;									#USE:
use File::Slurp;														#USE:
use GD::Graph::lines;												#USE:
use GD::Graph::bars;												#USE:
use Statistics::TTest;											#USE:
use Statistics::R;													#USE:
use Statistics::ChisqIndep;									#USE:
use List::Util qw(max sum min shuffle);			#USE:
use Math::NumberCruncher;										#USE:

																						#B<SCRIPT INFO: "bk_biology", "0.0.22", "20180611">
my ($varNam, $varRev, $varUpd, $varDev)		 = ("bk_biology", "0.0.22", "20180611", 0);

unless($varDev eq 1){												#B<GROUP: Miscellaneous>	
	sub bk_biology_version {									#METHOD: B<bk_biology_version> Tested:20131218. Time:0 sec.
		return 																	"$varNam\tv$varRev\tUpdated:\t$varUpd";
		
																						#ARGS: None
																						#RETURN: Scalar with Script Information.
																						#DESCRIPTION: Provides versioning information for the script.
	}
	
	sub bk_biology_blastToHash{								#METHOD: B<bk_biology_blastToHash> Tested:20131219. Time:54.5K hits/sec. 
		my ($file)															= (@_);
		my (%outBlast, @hit, $key);
		open																		(my $oIN, "<", "$file") 
			or die																"Could not open Input File - $!";
				while																(<$oIN>) { 
					chomp 														$_;
					if																($_ =~ /\# Query:/){ my @tmp1 = split " ",$_; $key = $tmp1[2]; }
					elsif															($_ !~ /\#/) { push @hit, "$_"; }
					elsif															($_ =~ /\# BLAST/) {
						if															(@hit>0) {
							my @tmp1 											= split "\t",$hit[0];
							push													@tmp1, @hit."";
							$outBlast{$key}							 = [@tmp1];
						}
						$key 														= ""; 
						@hit														= ();
					}
				}	
				close 															$oIN;
		return 																	%outBlast;
		
																						#ARGS: BLAST .out file path. *Out Tabular Formats only.
																						#RETURN: Hash of BLAST results by query id.
																						#DESCRIPTION: This module loads a BLAST tabular out file and
																						#loads it into a hash for access to the hit information by 
																						#query id. In most cases here the query id is the read ID but
																						#can be the name of any sequence in a fasta.
	}
}

unless($varDev eq 1){												#B<GROUP: Sequence Operations>	
	sub bk_biology_reverseComplement{					#METHOD: B<bk_biology_reverseComplement> Tested:20140206. Time:0 sec. 
		my ($dna)																= (@_);
		my $revcomp 														= reverse($dna);
		$revcomp 																=~ tr/ACGTacgt/TGCAtgca/;
		$revcomp 																= join "",$revcomp;
		return 																	$revcomp;

																						#ARGS: None.
																						#RETURN: Returns the reverse complement of the given string.
																						#DESCRIPTION: This module accepts a string of DNA characters and
																						#returns its reverse complement. *No error checking is done for
																						#Ns or non-ATCG characters.
	}
	
	sub bk_biology_removeDuplicates{					#METHOD: B<bk_biology_removeDuplicates> 
		my (%args)															= (@_);

		my $fork_manager												= new Parallel::ForkManager($args{Core});	
		my $out_metrics													= $args{Metrics_File};
				print $out_metrics									"DUPLICATE REMOVAL\nSAMPLE\tPAIRS REMOVED\tMAX THRESHOLD\n";
		foreach my $sample											(@{$args{"Samples"}}){
			$fork_manager													-> start and next;
				my (%read1)													= bk_biology_fastqToHash($args{Path}."/".$sample."/".$sample."_R1.inp");
				my (%read2)													= bk_biology_fastqToHash($args{Path}."/".$sample."/".$sample."_R2.inp");	
				my (%read_out1)											= (%read1);
				my (%read_out2)											= (%read2);
				if 																	($args{Method} eq "ends"){
					foreach my $key										(keys(%read1)){
						$read1{$key}[1]									= substr($read1{$key}[1],1,25);
						$read2{$key}[1]									= substr($read2{$key}[1],length($read2{$key}[1])-25,length($read2{$key}[1]));
					}
				}
				if 																	($args{Method} eq "nomer"){
					foreach my $key										(keys(%read1)){
						$read1{$key}[1]									= substr($read1{$key}[1],50,length($read1{$key}[1]));
						$read2{$key}[1]									= substr($read2{$key}[1],1,length($read2{$key}[1])-50);
					}
				}
				my (%unique);
						my (%exact, @temp_exact);
						foreach my $key								 (keys %read1){
							if														($read1{$key}[1] and $read2{$key}[1]){ 
								my $concatemer							= $read1{$key}[1].$read2{$key}[1];
								$exact{$key}								= $concatemer;								
								push												@temp_exact, $concatemer;
							}
						}													 
						my (@tmp_unique)								= bk_basic_uniqueArray(@temp_exact);										
						foreach my $reads							 (@tmp_unique){ $unique{$reads} = [0, ""]; }			
						foreach my $key								 (keys %exact){
							$unique{$exact{$key}}[0]++; 
							$unique{$exact{$key}}[1]			= "$key $unique{$exact{$key}}[1]";					
						}
				my (@totals, $threshold)						= (0) x 1000000, 0;
						my (@counts);
						foreach my $key								 (keys %unique){
							$totals[$unique{$key}[0]]++;
							push													@counts, $unique{$key}[0];
						}
						while														($totals[@totals-1] eq 0) { pop(@totals); }
						my ($deviation)									= Math::NumberCruncher::StandardDeviation(\@counts);
						$threshold											= bk_basic_divide
																						( numerator						=> sum(@counts),
																							denominator					=> @counts."",
																							significant_digits 	=> 2 						)+($deviation*3); 
						$threshold												= bk_basic_round
																						( significant_digits	=> 0,
																							number							=> $threshold	 );		
						my @data;
							for														(my $i = 0; $i < @totals; $i++) {	
								$data[0][$i]								= $i;
								$data[1][$i]								= $totals[$i];
							}	 
							my $graph 										= GD::Graph::bars->new(800, 400);
								$graph											-> set( 
									x_label									 	=> 'Number of Replicates',
									y_label									 	=> 'Reads',
									title											=> 'Reads Vs. Replicates',
									y_max_value								=> $data[1][2],
									x_max_value								=> 100,
									y_tick_number							=> 10,
									y_label_skip							=> 2,
									x_tick_number							=> 100,
									x_label_skip							=> 5,
									dclrs			 								=> ['lblue'],
									bar_spacing								=> 2,
									transparent								=> 0,
									) or die $graph						-> error;
								my $gd 											= $graph -> plot(\@data) or die $graph -> error;
								open												(IMG, '>'.$args{Path}."/_images/".$sample.'_dups.png') or die $!;
								binmode 										IMG;
								print 											IMG $gd -> png;	
				my (%rem_read1, %rem_read2);
						foreach my $key								 (keys %unique){
							if														($unique{$key}[0]>$threshold){ 
								my @keys										= split " ", $unique{$key}[1];
								for(my $i=$threshold-1; $i<@keys; $i++){
									$rem_read1{$keys[$i]}			= $read_out1{$keys[$i]};
									$rem_read2{$keys[$i]}			= $read_out2{$keys[$i]};
									delete 										$read_out1{$keys[$i]};
									delete 										$read_out2{$keys[$i]};								 
								}
							}
						} 
				bk_biology_hashToSequence					 (	Hash_Ref	=> \%read_out1,
																							Out_File	=> $args{Path}."/".$sample."/".$sample."_R1.inp",
																							Seq_Type	=> "fastq"																				);		
				bk_biology_hashToSequence					 (	Hash_Ref	=> \%read_out2,
																							Out_File	=> $args{Path}."/".$sample."/".$sample."_R2.inp",
																							Seq_Type	=> "fastq"																				); 
				bk_biology_hashToSequence					 (	Hash_Ref	=> \%rem_read1,
																							Out_File	=> $args{Path}."/_sequence/".$sample."_dup_R1.bs",
																							Seq_Type	=> "fastq"																				);
				bk_biology_hashToSequence					 (	Hash_Ref	=> \%rem_read2,
																							Out_File	=> $args{Path}."/_sequence/".$sample."_dup_R2.bs",
																							Seq_Type	=> "fastq"																				);										 
				my $rem_keys												= keys(%rem_read1)."";
				my $remaining												= keys(%read1)."";
				print $out_metrics 									"$sample\t$rem_keys\t$threshold\n";
			$fork_manager													-> finish; 
		}
		$fork_manager														-> wait_all_children;

																						#ARGS: 
																						#	Metrics_File			=> $oMET,
																						#	 Path 						=> $varPat,
																						#	 Core 						=> $varCor,		
																						#		Core 						=> $varCor,	
																						# 	Samples 				=> \@arrSAM,											
																						#RETURN: None.
																						#DESCRIPTION: This module removes exact duplicates over a specified
																						#maximum depth (default 3).
	}
	
	sub bk_biology_AAtoNuc {									#METHOD: B<bk_biology_AAtoNuc>#needs ARGS
		my $tmpAA 															= shift; 
		my $tmpNuc 															= "";
		if 	  																	($tmpAA =~ /A/) { $tmpNuc="GCG"; } 
		elsif 																	($tmpAA =~ /F/) { $tmpNuc="TTT"; }
		elsif 																	($tmpAA =~ /L/) { $tmpNuc="CTG"; }
		elsif 																	($tmpAA =~ /I/) { $tmpNuc="ATT"; }
		elsif 																	($tmpAA =~ /V/) { $tmpNuc="GTG"; }
		elsif 																	($tmpAA =~ /P/) { $tmpNuc="CCG"; }
		elsif 																	($tmpAA =~ /T/) { $tmpNuc="ACC"; }
		elsif 																	($tmpAA =~ /Y/) { $tmpNuc="TAT"; }
		elsif 																	($tmpAA =~ /H/) { $tmpNuc="CAT"; }
		elsif 																	($tmpAA =~ /Q/) { $tmpNuc="CAG"; }
		elsif 																	($tmpAA =~ /N/) { $tmpNuc="AAC"; }
		elsif 																	($tmpAA =~ /K/) { $tmpNuc="AAA"; }
		elsif 																	($tmpAA =~ /D/) { $tmpNuc="GAT"; }
		elsif 																	($tmpAA =~ /E/) { $tmpNuc="GAA"; }
		elsif 																	($tmpAA =~ /C/) { $tmpNuc="TGC"; }
		elsif 																	($tmpAA =~ /W/) { $tmpNuc="TGG"; }
		elsif 																	($tmpAA =~ /R/) { $tmpNuc="CGC"; }
		elsif 																	($tmpAA =~ /S/) { $tmpNuc="AGC"; }
		elsif 																	($tmpAA =~ /G/) { $tmpNuc="GGC"; }
		elsif 																	($tmpAA =~ /M/) { $tmpNuc="ATG"; }
		return 																	$tmpNuc;
}

	sub bk_biology_readMateCorrection {				#METHOD: B<bk_biology_readMateCorrection> #needs ARGS
		my (%args)															= (@_);	
		print $args{Sample}." ".time."\n";
		my ($out_cor)														= $args{Metrics};
		my (%read1)															= bk_biology_fastqToHash(@{$args{Reads}}[0]);
		my (%read2)															= bk_biology_fastqToHash(@{$args{Reads}}[1]);
		my (%corr_fas);
		my (@avg_error, @avg_length);
		my (@cnt_length)												= (0) x 250;
		my (%indel_r1, %indel_r2);
		my (%nomatch_r1, %nomatch_r2);
			my $fork_manager											= new Parallel::ForkManager($args{Core});
			foreach my $key												(keys %read1){
					my ($r1)													= $read1{$key}[1];
					my ($r2)													= bk_biology_reverseComplement($read2{$key}[1]);
					my ($cat)													= "";
					my ($flg_hit)											= 0;
					for																(my $i=0; $i<(length($r1)-7); $i++){
						my ($search)										= substr($r1,$i,7);
						if															($r2 =~ /$search/){ 
							my (@tmp_r1)									= split /$search/,$r1;
							my (@tmp_r2)									= split /$search/,$r2;
							if														($tmp_r1[1]){ $r1 = $search.$tmp_r1[1]; }
							if														($tmp_r2[1]){ $r2 = $search.$tmp_r2[1]; }
							my (@arr_r1)									= split //,$r1;
							my (@arr_r2)									= split //,$r2;
							my (@reads)										= ((@arr_r1+0),(@arr_r2+0));
							my ($length)									= max(@reads);
							my (@mm)											= (0) x 250;
							for														(my $y=0; $y<$length; $y++){
								if													($arr_r1[$y] and $arr_r2[$y]){
									if												($arr_r1[$y] =~ /[ATCGatcg]/ and $arr_r2[$y] =~ /[ATCGatcg]/ and $arr_r1[$y] eq $arr_r2[$y]){
										$mm[$y]									= 1;
									}
								}else{
									$mm[$y]										= 0;
								}
							}
							while													($mm[@mm-1] eq 0){ pop @mm; }
							my @err												= grep { $_ eq 0 } @mm;
							if														(@err <= 10){
								push 												@avg_error, (@err+0);
								my ($out_seq)								= "";
								for													(my $y=0; $y<@mm; $y++){		
									if												($y < (@mm*0.6)){
										$out_seq								= $out_seq.$arr_r1[$y];
									}else{
										$out_seq								= $out_seq.$arr_r2[$y];
									}
								}		
								push 												@avg_length, length($out_seq);
								$cnt_length									[length($out_seq)]++;
								$corr_fas{$key}[1]					= $out_seq;
							}else{
								$indel_r1{$key}							= $read1{$key};
								$indel_r2{$key}							= $read2{$key};
								delete											$read1{$key};
								delete											$read2{$key};
							}
							$flg_hit											= 1;
							last;
						}
					}
					if																($flg_hit ne 1){
							$nomatch_r1{$key}							= $read1{$key};
							$nomatch_r2{$key}							= $read2{$key};
							delete												$read1{$key};
							delete												$read2{$key};
					}else															{ $flg_hit = 0; }
			}
		bk_biology_hashToSequence								(	Hash_Ref	=> \%corr_fas,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_corr.fasta",
																							Seq_Type	=> "fasta"																						);		
		bk_biology_hashToSequence								(	Hash_Ref	=> \%indel_r1,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_R1_indel.bs",
																							Seq_Type	=> "fastq"																						);	
		bk_biology_hashToSequence								(	Hash_Ref	=> \%indel_r2,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_R2_indel.bs",
																							Seq_Type	=> "fastq"																						);	
		bk_biology_hashToSequence								(	Hash_Ref	=> \%nomatch_r1,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_R1_nomatch.bs",
																							Seq_Type	=> "fastq"																						);	
		bk_biology_hashToSequence								(	Hash_Ref	=> \%nomatch_r2,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_R2_nomatch.bs",
																							Seq_Type	=> "fastq"																						);	
		my ($avg_err)														= bk_basic_divide
																								( numerator						=> sum(@avg_error),
																									denominator					=> @avg_error+0,
																									significant_digits 	=> 2 				 				);	
		my ($avg_len)														= bk_basic_divide
																								( numerator						=> sum(@avg_length),
																									denominator					=> @avg_length+0,
																									significant_digits 	=> 2 				 				);
		my ($dev_err)														= bk_basic_round(	significant_digits => 2,
																															number => Math::NumberCruncher::StandardDeviation(\@avg_error));
		my ($dev_len)														= bk_basic_round(	significant_digits => 2,
																															number => Math::NumberCruncher::StandardDeviation(\@avg_length));
		print $out_cor													"$args{Sample}\t".keys(%corr_fas)."\t".keys(%nomatch_r1)."\t".keys(%indel_r1)."\t$avg_err\t$dev_err\t$avg_len\t$dev_len\n";
		foreach my $value												(@cnt_length){ $value = $value."\n"; }
		write_file															($args{Path}."/".$args{Sample}."/align_length.txt", @cnt_length);
		my @data;
			for																		(my $i = 0; $i < @cnt_length; $i++) {	
				$data[0][$i]												= $i;
				$data[1][$i]												= $cnt_length[$i];
			}	 
			my $graph 														= GD::Graph::bars->new(800, 400);
				$graph															-> set( 
					x_label									 					=> 'Length bp',
					y_label									 					=> 'Reads',
					title															=> 'Reads Vs. Length',
					y_max_value												=> max(@{$data[1]}),
					x_max_value												=> 251,
					y_tick_number											=> 10,
					y_label_skip											=> 2,
					x_tick_number											=> 100,
					x_label_skip											=> 5,
					dclrs			 												=> ['lblue'],
					bar_spacing												=> 2,
					transparent												=> 0,
					) or die $graph										-> error;
				my $gd 															= $graph -> plot(\@data) or die $graph -> error;
				open																(IMG, '>'.$args{Path}."/_images/".$args{Sample}.'_len.png') or die $!;
				binmode 														IMG;
				print 															IMG $gd -> png;	
	}

}

unless($varDev eq 1){												#B<GROUP: FASTQ Operations>
	sub bk_biology_randInp {	 								#METHOD: B<bk_biology_randInp>
		my (%args) 															= (@_);
		my $forkManager													= new Parallel::ForkManager($args{Core});
		my $oMET																= $args{Metrics_File};
				print $oMET													"Sequence\nSample\tTotal\tInput\n";
		foreach my $sample 											(@{$args{Files}}){
			$forkManager													-> start and next;
				my $folder													= $sample;
				if																	($args{Folder}){ $folder = $args{Folder}; }
				my @indeces													= <$args{Path}/$folder/*_I*uz>;				
				my @reads														= <$args{Path}/$folder/*_R*uz>;		
				if																	($args{Folder}){ @reads = grep { $_ =~ /$sample/ } @reads; }
				if																	($args{Folder}){ @indeces = grep { $_ =~ /$sample/ } @indeces; }
				my %read1														= bk_biology_fastqToHash($reads[0]);
				my (%read2, %index1, %index2);
				if 																	($reads[1]) 	{ %read2 	= bk_biology_fastqToHash($reads[1]); }
				if 																	($indeces[0]) { %index1 = bk_biology_fastqToHash($indeces[0]); }
				if																	($indeces[1]) { %index2 = bk_biology_fastqToHash($indeces[1]); }				
				my (%out_r1, %out_r2);
				my (%out_i1, %out_i2);		
				my $cnt_crap												= 0;
				my $cnt_num													= 0;
				if																	($reads[1]){
					if 																($args{Input_Number} eq 0){ $args{Input_Number} = keys(%read1); }
					foreach my $key 									(shuffle (keys %read1)){ 
						if															($read1{$key}[1]=~/NNNNNNNNNN/ or $read2{$key}[1]=~/NNNNNNNNNN/){ $cnt_crap++; next; }
						if															($cnt_num > $args{Input_Number}){ last; }
						if															($read1{$key} and $read2{$key}){
						if															(length($read1{$key}[1]) > 10 and length($read2{$key}[1]) > 10){
							$out_r1{$key}									= $read1{$key};
							$out_r2{$key}									= $read2{$key};
							if 														($indeces[0]) { $out_i1{$key} = $index1{$key}; }
							if 														($indeces[1]) { $out_i2{$key} = $index2{$key}; }
						}
						}
						$cnt_num++;
					}
				}else{
					if 																($args{Input_Number} eq 0){ $args{Input_Number} = keys(%read1); }
					foreach my $key 									(shuffle (keys %read1)){ 
						if															($cnt_num > $args{Input_Number}){ last; }
						$out_r1{$key}										= $read1{$key};
						if 															($indeces[0]) { $out_i1{$key} = $index1{$key}; }
						$cnt_num++;
					}
				}
				print "$cnt_crap sequences removed from $sample for being total or largely garbage.\n";
				bk_biology_hashToSequence						(	Hash_Ref	=> \%out_r1,
																							Out_File	=> $args{Path}."/".$folder."/".$sample."_R1.inp",
																							Seq_Type	=> "fastq"																				);	
				if 																	($reads[1]) {
					bk_biology_hashToSequence					(	Hash_Ref	=> \%out_r2,
																							Out_File	=> $args{Path}."/".$folder."/".$sample."_R2.inp",
																							Seq_Type	=> "fastq"																				);
				}
				if 																	($indeces[0]) {
					bk_biology_hashToSequence					(	Hash_Ref	=> \%out_i1,
																							Out_File	=> $args{Path}."/".$folder."/".$sample."_R1.ind",
																							Seq_Type	=> "fastq"																				);	
				}
				if 																	($indeces[1]) {
					bk_biology_hashToSequence					(	Hash_Ref	=> \%out_i2,
																							Out_File	=> $args{Path}."/".$folder."/".$sample."_R2.ind",
																							Seq_Type	=> "fastq"																				);
				}
					print $oMET											 	$sample."\t".(keys (%read1))."\t".(keys (%out_r1))."\n";
			$forkManager													-> finish; 
		}
		$forkManager														-> wait_all_children;
		my @inMET 															= bk_basic_uniqueArray(read_file($args{Path}."/_metrics/seq.met"));
		write_file 															($args{Path}."/_metrics/seq.met", @inMET);
		
																						#ARGS: 
																						#	 Metrics_File => \@outBGE,
																						#	 Path => $varPat,
																						#	 Core => $varCor,
																						#	 Files => @Files,
																						#	 Input_Number => $varInp
																						#RETURN: None.
																						#DESCRIPTION: Accepts a list of files and reduces them to requested number	
																						#of fastq reads, access beginning and end counts and reports to a metrics 
																						#file. Process runs in parallel by the number of core passed.
	}
	
	sub bk_biology_fastqToHash {							#METHOD: B<bk_biology_fastqToHash> Tested:20131218. Time:67K reads(fastq)/sec.
		my ($fastq) 														= (@_); 
		
		my ($name,$seq,$qual,$nam2) 						= ("") x 4;
		my ($cntSeq) 														= (0);
		my 																			(%outSeq);
		open my																 $iUZ, '<', "$fastq" or die "Can't open $fastq: $!";
				while																(<$iUZ>){ 
					chomp 														$_;
					if																($cntSeq eq 0) { 				
						my @tmpNam 											= split " ", $_; 
						$nam2														= $tmpNam[1]; 
						if															($fastq =~ /SRR/){ 
							if														($fastq =~ /\_R1/){
								$nam2												= "1:N:0:8";
							}else{
								$nam2												= "2:N:0:8";
							}
							shift @tmpNam; 
							$tmpNam[0]										=~ s/\./\:/g;
						}
						$name														= $tmpNam[0];
						$name														=~ s/@//; 
						$cntSeq++; }
					elsif 														($cntSeq eq 1) { $seq = $_; $cntSeq++; }
					elsif															($cntSeq eq 3) { 
						$qual 													= $_; 
						if 															(substr($qual,0,1) eq "@"){
							$qual													= "A".substr($qual,1,length($qual));
						}
						$outSeq{$name} 									= [$nam2, $seq, $qual];
						$cntSeq													= 0; 
						$name 													= $seq = $qual = $nam2 = ""; }
					else															{ $cntSeq++; }		
				} 
				close 															$iUZ;
		foreach my $key 												(keys %outSeq){ if ($outSeq{$key}[1]){}else{ delete $outSeq{$key}; }}
		return 																	%outSeq;
		
																						#ARGS: A fastq file path.
																						#RETURN: A hash containing the fastq from the file.
																						#DESCRIPTION: This module loads a fastq sequence file
																						#into a hash for easy read name retieval in other methods.
	}
	
	sub bk_biology_fastaToHash {							#METHOD: B<bk_biology_fastaToHash> Tested:20131218. Time:67K reads(fastq)/sec.
		my ($fasta) 														= (@_); 
		my ($name,$seq, $nam2) 									= ("") x 3;
		my 																			(%outSeq);
		open my																 	$iUZ, '<', "$fasta" or die "Can't open $fasta: $!";
				while																(<$iUZ>){ 
					chomp $_;
					if																(substr($_,0,1) =~ /\>/) { 	
						my @tmpNam 											= split " ", $_; 
						$name														= $tmpNam[0];
						$nam2														= $tmpNam[1];
						$seq														= "";
					}else{ 	
						$seq														.= $_;
						$outSeq{$name} 									= [$nam2, $seq];
					}
				} 
				close 															$iUZ;
		foreach my $key 												(keys %outSeq){ if ($outSeq{$key}[1]){}else{ delete $outSeq{$key}; }}
		return 																	%outSeq;
		
																						#ARGS: A fastq file path.
																						#RETURN: A hash containing the fastq from the file.
																						#DESCRIPTION: This module loads a fastq sequence file
																						#into a hash for easy read name retieval in other methods.
	}
	
	sub bk_biology_hashToSequence{						#METHOD: B<bk_biology_hashToSequence>
		my (%args)															= (@_);

		my ($out_file);
				open $out_file, 										">", $args{Out_File} 
					or die														"Could not open file - $!";
				my (%tmp_hash)											= %{$args{Hash_Ref}};
				foreach my $key										 (sort { $a cmp $b } keys %tmp_hash){
					my ($out_line)										= "";			
					if																($args{Seq_Type} eq "fasta"){
						my $name_2											= "";
						if 															($tmp_hash{$key}[0]){ $name_2 = $tmp_hash{$key}[0]; }
						if 															($tmp_hash{$key}[1]){
							if														($key =~ />/){
								$out_line 									= "$key $name_2\n$tmp_hash{$key}[1]\n";	
							}else{
								$out_line 									= "\>$key $name_2\n$tmp_hash{$key}[1]\n";	
							}
						}
					}elsif														($args{Seq_Type} eq "fastq"){ 
						my $name_2											= "";
						if 															($tmp_hash{$key}[0]){ $name_2 = $tmp_hash{$key}[0]; }
						if 															($tmp_hash{$key}[1]){
							$out_line 										= "\@$key $name_2\n$tmp_hash{$key}[1]\n+\n$tmp_hash{$key}[2]\n"; 
						}
					}		
					print $out_file $out_line;
				}
				
																						#ARGS: 
																						#	Hash_Ref	=> \%hash,
																						# Out_File	=> "path",
																						# Seq_Type	=> "fasta"
																						#RETURN: None.
																						#DESCRIPTION: This module writes sequence information from a hash
																						#to a file of a specified type. (fasta or fastq).		
	}	
}

unless($varDev eq 1){												#B<GROUP: PHAGE Methods> 
	sub bk_biology_PhageArray {								#METHOD: B<bk_biology_PhageArray> #needs ARGS
			my (%args) 															= (@_);
			my (%arr_lib);
			my (@oligos);
			my (@oligos1);
			
			my $out_met															= $args{Metrics};
			foreach my $file 												(@{$args{Alignments}}){	
				my ($file_name)												= bk_basic_removeFileType(bk_basic_returnFile($file));
				my (@tmp_file)												= read_file($file);
				my ($tmp_ssn)													= ("");
					foreach my $line										(@tmp_file){
						chomp															$line;
						if																($line =~ /\>/){
							$line														=~ s/>//;
							$line														=~ s/ |-|\.|\|/_/g;
							$tmp_ssn												= $file_name.":".$line;
									chomp												$tmp_ssn;
						}else{
							my (@tmp_amino) 								= split //,$line;
							my (@tmp_seq);
									for 												(my $i=0; $i<@tmp_amino; $i++){
										if												($tmp_amino[$i] !~ /-/){
											push										@tmp_seq, [ $tmp_amino[$i], ($i+1) ];
										}
									}
									pop @tmp_seq;
									for													(my $i=0; $i<@tmp_seq; $i++){
										if												(@tmp_seq-$i>=$args{Window}){
											my $tmp_end 						= $i+($args{Window}-1);
											my @tmp_arr 						= @tmp_seq[$i..$tmp_end];
											my $tmp_name 						= ">".$tmp_ssn.":".$tmp_arr[0][1].":".$tmp_arr[$args{Window}-1][1].":".$i.":".$tmp_end;
											my $tmp_ase 						= join "", map $_->[0], @tmp_arr;
											my $tmp_nse 						= "";
											for											(my $i=0; $i<@tmp_arr; $i++){
												$tmp_nse 							= $tmp_nse.bk_biology_AAtoNuc($tmp_arr[$i][0])
											}
											$arr_lib{$tmp_name}[0]	= $tmp_ase;
											push 										@oligos, [$tmp_name, $tmp_ase, $tmp_nse]; 
											$arr_lib{$tmp_name}[0]	= $tmp_nse;
										}
										$i=$i+$args{Slide}-1;
									}
							$tmp_ssn = "";								
						}
					}
			}
			print $out_met													"ARRAY\nTASK\tCOUNT\t% TOTAL\n";
			print $out_met													"Alignments\t".@{$args{Alignments}}."\n";
			print $out_met													"Total Oligos\t".keys(%arr_lib)."\n";
			my @uni_oligos													= bk_basic_unique2DArray(1, @oligos);
			my @fin_oligos													= reverse(grep { $_->[1] !~ /\*|\?|X/ } @uni_oligos);
			my $dup_rem															= (@oligos-@uni_oligos);
			my $stp_rem															= (@uni_oligos-@fin_oligos);
			my $dup_per															= bk_basic_divide
																								( numerator						=> $dup_rem,
																									denominator					=> @oligos+0,
																									significant_digits 	=> 4 				 );		
			my $stp_per															= bk_basic_divide
																								( numerator						=> $stp_rem,
																									denominator					=> @oligos+0,
																									significant_digits 	=> 4 				 );
			my $fin_per															= bk_basic_divide
																								( numerator						=> @fin_oligos+0,
																									denominator					=> @oligos+0,
																									significant_digits 	=> 4 				 );
			print $out_met													"Duplicates Removed\t".$dup_rem."\t".$dup_per."\n";
			print $out_met													"Partial/Stop Removed\t".$stp_rem."\t".$stp_per."\n";
			print $out_met													"Array Size\t".@fin_oligos."\t".$fin_per."\n";
			for																			(my $i=0; $i<@fin_oligos; $i++){	
				my @tmp_cnt 													= grep { $_->[1] eq $fin_oligos[$i][1] } @oligos;
				$fin_oligos[$i][0] 										= $fin_oligos[$i][0].":".@tmp_cnt;
			}
			print "Dup counts fixed.\n";
			my (@out_oligos, @out_aa, @out_nuc);
					for																	(my $i=0; $i<@fin_oligos; $i++){										
						push 															@out_oligos, 	$fin_oligos[$i][0]."\t".$fin_oligos[$i][1]."\t".$fin_oligos[$i][2]."\n";									
						push 															@out_nuc, 		$fin_oligos[$i][0]."\n".$fin_oligos[$i][2]."\n";		
						push 															@out_aa, 			$fin_oligos[$i][0]."\n".$fin_oligos[$i][1]."\n";	
					}
					my @exclude_vals										= split ",", $args{Exclude};
					foreach my $val											(@exclude_vals){
						@out_oligos												= grep { $_ !~ /$val/ } @out_oligos;
						@out_nuc													= grep { $_ !~ /$val/ } @out_nuc;
						@out_aa														= grep { $_ !~ /$val/ } @out_aa;
					}
					write_file													($args{Path}."/_array/oligos.txt",@out_oligos);
					write_file													($args{Path}."/_library/ref.fasta",@out_nuc);
					write_file													($args{Path}."/_library/refp.fasta",@out_aa);
					bk_wrappers_formatBlast							(	Path_Fasta => $args{Path}."/_library/ref.fasta",
																								Database_Type => "nucl" );
					bk_wrappers_formatBlast							(	Path_Fasta => $args{Path}."/_library/refp.fasta",
																								Database_Type => "prot" );
	}
	
	sub bk_biology_PhageCleaning {						#METHOD: B<bk_biology_PhageCleaning> #needs ARGS
		my (%args) 															= (@_);		
		
		my (%reads)															= bk_biology_fastaToHash("$args{Path}/_sequence/$args{Sample}_corr.fasta");
		my (%chimera, %short, %good);
		my (%rep1, %rep2, %rep3);
		my (%bad_rep);
		my ($out_met)														= $args{Metrics};
		my (@cnt_insert)												= (0) x 250;
		foreach my $key													(keys(%reads)){
			my ($name)														= ($key);
			$name																	=~ s/\>//;
			if																		(index($reads{$key}[1],"GAATTC") > index($reads{$key}[1],"AAGCTT")){
				$reads{$key}[1]											= bk_biology_reverseComplement($reads{$key}[1]);
			}
			$reads{$key}[1]												=~ s/GAATTCG|AAGCTT/_/g;			
			my @tmp_insert												= split "_", $reads{$key}[1];
			if 																		(@tmp_insert eq 3){ 
				$good{$name}[0]											= $reads{$key}[0];
				$good{$name}[1]											= $tmp_insert[1];
				if																	($reads{$key}[1] =~ /GGTTCC/){
					if																($reads{$key}[1] =~ /GTACTC/){
						$rep1{$name}[0]									= $reads{$key}[0];
						$rep1{$name}[1]									= $tmp_insert[1];
					}elsif														($reads{$key}[1] =~ /CCTATC/ or $reads{$key}[1] =~ /CTTACG/){
						$bad_rep{$name}[0]							= $reads{$key}[0];
						$bad_rep{$name}[1]							= $tmp_insert[1];						
					}
				}elsif 															($reads{$key}[1] =~ /CCAATT/){
					if																($reads{$key}[1] =~ /CCTATC/){
						$rep2{$name}[0]									= $reads{$key}[0];
						$rep2{$name}[1]									= $tmp_insert[1];
					}elsif														($reads{$key}[1] =~ /GTACTC/ or $reads{$key}[1] =~ /CTTACG/){
						$bad_rep{$name}[0]							= $reads{$key}[0];
						$bad_rep{$name}[1]							= $tmp_insert[1];
					}
				}elsif 															($reads{$key}[1] =~ /GTCCAG/){ 
					if 																($reads{$key}[1] =~ /CTTACG/){ 
						$rep3{$name}[0]									= $reads{$key}[0];
						$rep3{$name}[1]									= $tmp_insert[1];
					}elsif														($reads{$key}[1] =~ /GTACTC/ or $reads{$key}[1] =~ /CCTATC/){ 
						$bad_rep{$name}[0]							= $reads{$key}[0];
						$bad_rep{$name}[1]							= $tmp_insert[1];						
					}
				}
																						$cnt_insert[length($tmp_insert[1])]++;
			}elsif 																(@tmp_insert > 3){ 
				$chimera{$name}[0]									= $reads{$key}[0];
				$chimera{$name}[1]									= $reads{$key}[1];				
			}else{
				$short{$name}[0]										= $reads{$key}[0];
				$short{$name}[1]										= $reads{$key}[1];	
			}
		}
		my ($assigned)													= keys(%rep1) + keys(%rep2) + keys(%rep3) + 0;
		my ($bad_pairing)												= keys(%bad_rep) + 0;
		my ($per_rep)														= bk_basic_divide
																								( numerator						=> $assigned,
																									denominator					=> keys(%good)+0,
																									significant_digits 	=> 4 				 				);
		my ($per_bad)														= bk_basic_divide
																								( numerator						=> $bad_pairing,
																									denominator					=> keys(%good)+0,
																									significant_digits 	=> 4 				 				);
		my ($per_good)													= bk_basic_divide
																								( numerator						=> keys(%good)+0,
																									denominator					=> keys(%reads)+0,
																									significant_digits 	=> 4 				 				);
		my ($per_short)													= bk_basic_divide
																								( numerator						=> keys(%short)+0,
																									denominator					=> keys(%reads)+0,
																									significant_digits 	=> 4 				 				);
		my ($per_chimera)												= bk_basic_divide
																								( numerator						=> keys(%chimera)+0,
																									denominator					=> keys(%reads)+0,
																									significant_digits 	=> 4 				 				);
		print $out_met 													"$args{Sample}\t".keys(%reads)."\t".keys(%good)."\t$per_good\t".keys(%short)."\t$per_short\t".keys(%chimera)."\t$per_chimera\t$assigned\t$per_rep\t$bad_pairing\t$per_bad\n";
		if																			($assigned>0){
			bk_biology_hashToSequence							(	Hash_Ref	=> \%rep1,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_insert_rep1.fasta",
																							Seq_Type	=> "fasta"																												);	
			bk_biology_hashToSequence							(	Hash_Ref	=> \%rep2,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_insert_rep2.fasta",
																							Seq_Type	=> "fasta"																													);	
			bk_biology_hashToSequence							(	Hash_Ref	=> \%rep3,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_insert_rep3.fasta",
																							Seq_Type	=> "fasta"																												);				
			bk_biology_hashToSequence							(	Hash_Ref	=> \%bad_rep,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_bad_pairing.fasta",
																							Seq_Type	=> "fasta"																												);				
		}else{
			bk_biology_hashToSequence							(	Hash_Ref	=> \%good,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_insert.fasta",
																							Seq_Type	=> "fasta"																												);	
		}
		bk_biology_hashToSequence								(	Hash_Ref	=> \%chimera,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_chimeric.fasta",
																							Seq_Type	=> "fasta"																												);	
		bk_biology_hashToSequence								(	Hash_Ref	=> \%short,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_noinsert.fasta",
																							Seq_Type	=> "fasta"																												);	
		my @cnt_length													= read_file($args{Path}."/".$args{Sample}."/align_length.txt");
		foreach my $value												(@cnt_insert){ $value = $value."\n"; }
		write_file															($args{Path}."/".$args{Sample}."/insert_length.txt", @cnt_insert);
		foreach my $value												(@cnt_length){ chomp $value; }
		my @data;
			for																		(my $i = 0; $i < @cnt_length; $i++) {	
				$data[0][$i]												= $i;
				$data[1][$i]												= $cnt_length[$i];
				$data[2][$i]												= $cnt_insert[$i];
			}	 
			my $graph 										= GD::Graph::bars->new(800, 400);
				$graph											-> set( 
					x_label									 	=> 'Length bp',
					y_label									 	=> 'Reads',
					title											=> 'Reads Vs. Length',
					y_max_value								=> max(@{$data[1]}),
					x_max_value								=> 251,
					y_tick_number							=> 10,
					y_label_skip							=> 2,
					x_tick_number							=> 100,
					x_label_skip							=> 5,
					dclrs			 								=> ['lblue','red'],
					bar_spacing								=> 3,
					transparent								=> 0,
					) or die $graph						-> error;
				$graph 											-> set_legend("Corrected","Insert");
				my $gd 											= $graph -> plot(\@data) or die $graph -> error;
				open												(IMG, '>'.$args{Path}."/_images/".$args{Sample}.'_insert.png') or die $!;
				binmode 										IMG;
				print 											IMG $gd -> png;	

	}	
	
	sub bk_biology_PhageLibrary {							#METHOD: B<bk_biology_PhageLibrary> #needs ARGS
		my (%args) 															= (@_);
		bk_wrappers_paraTabBlast								( Core 					=> $args{Core},
																							Blast_Program => "blastx",
																							Database 			=> $args{Path}."/_library/refp.fasta",
																							GI_List 			=> "",
																							Query_Fasta 	=> $args{Path}."/_sequence/".$args{Sample}.".fasta",
																							Output_Format => "qseqid qlen qstart qend sseqid slen sstart send qcovs ppos",
																							Number_Hits 	=> 1,
																							Task 					=> "" 																													);
		my (%blast_in)													= bk_biology_blastToHash($args{Path}."/_sequence/".$args{Sample}.".fasta.out");
		my (%reads)															= bk_biology_fastaToHash($args{Path}."/_sequence/".$args{Sample}.".fasta");
		my (%nucl)															= bk_biology_fastaToHash("$args{Path}/_library/ref.fasta");		
		my (%prot)															= bk_biology_fastaToHash("$args{Path}/_library/refp.fasta");		
		my (%results);
		my (%seq_bstart, %seq_good, %seq_short);
		my (@cnt_hit)														= (0) x 250;
		my (@cnt_bstart)												= (0) x 250;
		my (@cnt_short)													= (0) x 250;
		my ($out_met)														= $args{Metrics};
		foreach my $key													(keys(%prot)){
			my (@initialize)											= (0) x 14;
			$results{$key}												= [@initialize];
			$results{$key}[12]										= $prot{$key}[1];
			$results{$key}[13]										= $nucl{$key}[1];
		}
		foreach my $key													(keys(%blast_in)){		
			my ($name)														= ">".$key;
			my ($sequence)												= $blast_in{$key}[0];
			my ($hit)															= ">".$blast_in{$key}[4];
			$results{$hit}[0]++;	
			# if																		($blast_in{$key}[6]-1 % 3 ne 0 or $blast_in{$key}[2]-1 % 3 ne 0){
				# print bk_biology_nucToAA($reads{$name}[1])."\n";
				# $results{$hit}[12]++;
				# $seq_bstart{$key}[0]								= "";
				# $seq_bstart{$key}[1]								= $reads{$name}[1];
																						# $cnt_bstart[length($reads{$name}[1])]++;
			if																		(abs($blast_in{$key}[7]-($blast_in{$key}[6]-1)) eq $blast_in{$key}[5]){
				$results{$hit}[1]++;
				$results{$hit}[2]++;
																						$cnt_hit[length($reads{$name}[1])]++;			
				if																	($blast_in{$key}[9] > 99.9){ 
					$results{$hit}[3]++;
				}else{
					$results{$hit}[4]++;
				}
				$seq_good{$key}[0]									= "";
				$seq_good{$key}[1]									= $reads{$name}[1];				
			}elsif																(abs($blast_in{$key}[7]-($blast_in{$key}[6]-1)) > ($blast_in{$key}[5]*0.67)){
				$results{$hit}[1]++;
				$results{$hit}[5]++;
				if																	($blast_in{$key}[9] > 99.9){ 
					$results{$hit}[6]++;
				}else{
					$results{$hit}[7]++;
				}
				$seq_good{$key}[0]									= "";
				$seq_good{$key}[1]									= $reads{$name}[1];	
			}elsif																(abs($blast_in{$key}[7]-($blast_in{$key}[6]-1)) > ($blast_in{$key}[5]*0.33)){
				$results{$hit}[1]++;
				$results{$hit}[8]++;
				if																	($blast_in{$key}[9] > 99.9){ 
					$results{$hit}[9]++;
				}else{
					$results{$hit}[10]++;
				}
				$seq_good{$key}[0]									= "";
				$seq_good{$key}[1]									= $reads{$name}[1];	
			}else{
				$results{$hit}[11]++;
				$seq_short{$key}[0]									= "";
				$seq_short{$key}[1]									= $reads{$name}[1];
																						$cnt_short[length($reads{$name}[1])]++;				
			}
		}		
		# my (@out_res)														= "OLIGO\tTOTAL\tHIT\tPTM\tPPM\tPMM\tGTM\tGPM\tGMM\tOTM\tOPM\tOMM\tSHORT\tBSTART\tPROT\tNUCL\n";
		my (@out_res)														= "OLIGO\tTOTAL\tHIT\tPTM\tPPM\tPMM\tGTM\tGPM\tGMM\tOTM\tOPM\tOMM\tSHORT\tPROT\tNUCL\n";
		for my $key 														( sort { $results{$b}[1] <=> $results{$a}[1] } keys %results ) {
			my ($out_values)											= join "\t", @{$results{$key}};
			push 																	@out_res, $key."\t".$out_values."\n";
		}
		my ($per_good)													= bk_basic_divide
																							( numerator						=> keys(%seq_good)+0,
																								denominator					=> keys(%reads)+0,
																								significant_digits 	=> 4 				 				);
		my ($per_short)													= bk_basic_divide
																							( numerator						=> keys(%seq_short)+0,
																								denominator					=> keys(%reads)+0,
																								significant_digits 	=> 4 				 				);
		# my ($per_bstart)												= bk_basic_divide
																							# ( numerator						=> keys(%seq_bstart)+0,
																								# denominator					=> keys(%reads)+0,
																								# significant_digits 	=> 4 				 				);
		my ($per_hits)													= bk_basic_divide
																							( numerator						=> keys(%blast_in)+0,
																								denominator					=> keys(%reads)+0,
																								significant_digits 	=> 4 				 				);
		# print $out_met $args{Sample}."\t".keys(%seq_good)."\t$per_good\t".keys(%seq_short)."\t$per_short\t".keys(%seq_bstart)."\t$per_bstart\t".keys(%blast_in)."\t$per_hits\t".keys(%reads)."\n";
		print $out_met $args{Sample}."\t".keys(%seq_good)."\t$per_good\t".keys(%seq_short)."\t$per_short\t".keys(%blast_in)."\t$per_hits\t".keys(%reads)."\n";
		write_file															( $args{Path}."/_reports/".$args{Sample}.".hit", @out_res );
		bk_biology_hashToSequence								(	Hash_Ref	=> \%seq_good,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_hit.fasta",
																							Seq_Type	=> "fasta"																						);
		bk_biology_hashToSequence								(	Hash_Ref	=> \%seq_short,
																							Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_short.fasta",
																							Seq_Type	=> "fasta"																						);	 
		# bk_biology_hashToSequence								(	Hash_Ref	=> \%seq_bstart,
																							# Out_File	=> $args{Path}."/_sequence/".$args{Sample}."_bstart.fasta",
																							# Seq_Type	=> "fasta"																						);	
		my (@tmp_name)													= split /\_insert/,$args{Sample};
		my @cnt_insert													= read_file($args{Path}."/".$tmp_name[0]."/insert_length.txt");
		foreach my $value												(@cnt_insert){ chomp $value; }
		my @data;
			for																		(my $i = 0; $i < @cnt_insert; $i++) {	
				$data[0][$i]												= $i;
				$data[1][$i]												= $cnt_insert[$i];
				$data[2][$i]												= $cnt_hit[$i];
				$data[3][$i]												= $cnt_short[$i];
				# $data[4][$i]												= $cnt_bstart[$i];
			}	 
			my $graph 										= GD::Graph::lines->new(800, 400);
				$graph											-> set( 
					x_label									 	=> 'Length bp',
					y_label									 	=> 'Reads',
					title											=> 'Reads Vs. Length',
					# line_types								=> [1, 1, 1, 1, 1],
					line_types								=> [1, 1, 1, 1],
					line_width								=> 1,
					y_max_value								=> max(@{$data[1]}),
					x_max_value								=> 251,
					y_tick_number							=> 10,
					x_label_position 					=> 1/2, 
					x_labels_vertical 				=> 1, 
					y_label_skip							=> 2,
					x_tick_number							=> 100,
					x_label_skip							=> 5,
					dclrs			 								=> ['lblue','red','lgreen','orange'],
					bar_spacing								=> 2,
					transparent								=> 0,
					box_axis 									=> 0, 
					r_margin 									=> 15, 
					) or die $graph						-> error;
				# $graph 											-> set_legend("Insert","Hit","Short","Bad Start");
				$graph 											-> set_legend("Insert","Hit","Short");
				$graph											-> set_legend_font(GD::gdMediumBoldFont);
				my $gd 											= $graph -> plot(\@data) or die $graph -> error;
				open												(IMG, '>'.$args{Path}."/_images/".$args{Sample}.'_final.png') or die $!;
				binmode 										IMG;
				print 											IMG $gd -> png;	
	}
	
	sub bk_biology_PhageAnalysis {						#METHOD: B<bk_biology_PhageAnalysis> #needs ARGS
		my (%args) 															= (@_);
		my (%summaries);
		my (%align_histo);
		my (%filter);
		my (%alignments)												= bk_biology_FastaAlignToHash( Path => $args{Path} );
		my (%out_metric);
		my ($out_met)														= $args{Metrics};
		my ($exp_thresh)												= $args{Exp};
		foreach my $sample											(@{$args{Sample}}){
			my (@hit_files)												= <$args{Path}/_reports/$sample*hit>;
			my $sam_name													= $sample;
			if																		($sample =~ /_lib|_LIB|_Lib/){ $sam_name = "LIBRARY"; }
			my (@tmp_files);
			foreach my $read_in										(@hit_files){
				my (@in_file)												= read_file($read_in);
						chomp														$in_file[0];
				my (%file);
				my (@headings)											= split /\t/,$in_file[0];
						shift														@in_file;
				foreach my $line										(@in_file){
					chomp $line;
					my (@values)											= split /\t/,$line;
					for																(my $i=0; $i<@headings; $i++){
					 $file{$values[0]}{$headings[$i]} = $values[$i];
					}
				}
				push 																@tmp_files, \%file;
			}
			my (%sam_sum);
			my (@r1_hit, @r2_hit, @r3_hit);
			my ($r1_tot, $r2_tot, $r3_tot)				= (0) x 3;
			foreach my $key												(keys %{$tmp_files[0]}){ push @r1_hit, $tmp_files[0]{$key}{HIT}; } 
					$r1_tot														= sum(@r1_hit);
			if 																		($tmp_files[1]){										
				foreach my $key											(keys %{$tmp_files[1]}){ 
					push 															@r2_hit, $tmp_files[1]{$key}{HIT}; 
				} 
				$r2_tot															= sum(@r2_hit); 
			}
			if 																		($tmp_files[2]){	
				foreach my $key											(keys %{$tmp_files[2]}){ 
					push 															@r3_hit, $tmp_files[2]{$key}{HIT}; 
				} 
				$r3_tot															= sum(@r3_hit);
			}
			foreach my $key												(keys %{$tmp_files[0]}){
				$sam_sum{$key}{OLIGO}								= $tmp_files[0]{$key}{OLIGO};
				my (@line_loc)											= split /\:/,$key;
				$sam_sum{$key}{ALIGN}								= $line_loc[0];
				$sam_sum{$key}{ALIGN}								=~ s/>//;
				$sam_sum{$key}{START}								= $line_loc[2];
				$sam_sum{$key}{END}									= $line_loc[3];
				$sam_sum{$key}{PROT}								= $tmp_files[0]{$key}{PROT};
				$sam_sum{$key}{NUCL}								= $tmp_files[0]{$key}{NUCL};
				$sam_sum{$key}{REP1_CNT}						= $tmp_files[0]{$key}{HIT};
				$sam_sum{$key}{REP1_PER}						= sprintf("%.10g", bk_basic_divide
																							( numerator						=> $tmp_files[0]{$key}{HIT},
																								denominator					=> $r1_tot,
																								significant_digits 	=> 8 				 								)+0.00000001);
				$sam_sum{$key}{REP2_CNT}						= 0;
				$sam_sum{$key}{REP2_PER}						= 0;
				$sam_sum{$key}{REP3_CNT}						= 0;
				$sam_sum{$key}{REP3_PER}						= 0;
				if																	($tmp_files[1]){ 
					$sam_sum{$key}{REP2_CNT}					= $tmp_files[1]{$key}{HIT};	
					$sam_sum{$key}{REP2_PER}					= bk_basic_divide
																							( numerator						=> $tmp_files[1]{$key}{HIT},
																								denominator					=> $r2_tot,
																								significant_digits 	=> 8 				 								)+0.00000001;					
				}
				if																	($tmp_files[2]){ 
					$sam_sum{$key}{REP3_CNT}					= $tmp_files[2]{$key}{HIT};
					$sam_sum{$key}{REP3_PER}					= bk_basic_divide
																							( numerator						=> $tmp_files[2]{$key}{HIT},
																								denominator					=> $r3_tot,
																								significant_digits 	=> 8 				 								)+0.00000001;						
				}
			my (@cnts)														= ($sam_sum{$key}{REP1_CNT},$sam_sum{$key}{REP2_CNT},$sam_sum{$key}{REP3_CNT});
			$sam_sum{$key}{CNT_AVG}								= bk_basic_round
																							(	significant_digits 	=> 0,
																								number 							=> Math::NumberCruncher::Mean(\@cnts));
			}
			$summaries{$sam_name}									= \%sam_sum;
		}
		my (%library)														= %{$summaries{LIBRARY}};
				delete 															$summaries{LIBRARY};	
		foreach my $key													(keys %library){
			my (@reps)														= ($library{$key}{REP1_PER},$library{$key}{REP2_PER},$library{$key}{REP3_PER});
			$library{$key}{LIB_AVG}								= bk_basic_round
																							(	significant_digits 	=> 8,
																								number 							=> Math::NumberCruncher::Mean(\@reps));
			$library{$key}{LIB_STD}								= bk_basic_round
																							(	significant_digits 	=> 8,
																								number 							=> Math::NumberCruncher::StandardDeviation(\@reps));
			my (@line_loc)												= split /\:/, $key;
			my ($start_loc)												= $line_loc[2];
			my ($end_loc)													= $line_loc[3];
			$library{$key}{CLUS}									= 0;	
			foreach my $line											(keys %library){ 
				my @test_loc												= split /\:/, $line;
				my ($start_lin)											= $test_loc[2];
				my ($end_lin)												= $test_loc[3];
					if																	($library{$key}{ALIGN} eq $library{$line}{ALIGN}){				
					if 																	($start_loc ~~ [$start_lin..$end_lin] or $end_loc ~~ [$start_lin..$end_lin]){
						$library{$key}{CLUS}++;	
					}
					}
			}
		}
		my (@align_files)													= <$args{Path}/_array/*.fasta>;
		foreach my $file													(@align_files){
			my (@lines)															= read_file($file);
					foreach 														(@lines){ chomp; }
			my (@names)															= grep { $_ =~ /\>/ } @lines;
			my (@seqs)															= grep { $_ !~ /\>/ } @lines;
			my (@tmp_histo)													= (0) x length($seqs[0]);
			my ($align_name)												= uc(bk_basic_removeFileType(bk_basic_returnFile($file)));
			$align_histo{$align_name}{LIB}{UNI}			= [@tmp_histo];
			$align_histo{$align_name}{LIB}{ACT}			= [@tmp_histo];
			$align_histo{$align_name}{LIB}{EXP}			= [@tmp_histo];
			foreach my $sample											(keys %summaries){
				$align_histo{$align_name}{$sample}		{ACT}	= [@tmp_histo];
				$align_histo{$align_name}{$sample}		{EXP}	= [@tmp_histo];
			}
		}
		foreach my $sample												(keys %summaries){
			foreach my $key													(keys %{$summaries{$sample}}){
				my (@tmp_values)											= split /\:/, $key;
						$tmp_values[0]										=~ s/\>//g;
				for																		(my $i=$tmp_values[2]; $i<$tmp_values[3]; $i++){
					if																	($summaries{$sample}{$key}{CNT_AVG} >= 1){
						$align_histo{uc($tmp_values[0])}			{$sample}{ACT}[$i]++;
					}
				}				
			}
		}
		my (@lib_hit);
		my ($lib_tot)														= 0;
		foreach my $key													(keys %library){ push @lib_hit, $library{$key}{CNT_AVG}; } 
				$lib_tot														= sum(@lib_hit);
		foreach my $sample											(keys %summaries){
			my (@sam_hit);
			my ($sam_tot)													= 0;
			foreach my $key												(keys %{$summaries{$sample}}){ push @sam_hit, $summaries{$sample}{$key}{CNT_AVG}; } 
					$sam_tot													= sum(@sam_hit);			
			foreach my $key												(keys %{$summaries{$sample}}){
				$summaries{$sample}{$key}{LIB_CNT}	= bk_basic_round
																							(	significant_digits 	=> 4,
																								number 							=> sprintf("%.10g", $library{$key}{CNT_AVG}));
				$summaries{$sample}{$key}{LIB_PER}	= bk_basic_round
																							(	significant_digits 	=> 4,
																								number 							=> sprintf("%.10g", $library{$key}{LIB_AVG}));
				$summaries{$sample}{$key}{REP1_EXP}	= bk_basic_divide
																							( numerator						=> $summaries{$sample}{$key}{REP1_PER},
																								denominator					=> $library{$key}{LIB_AVG},
																								significant_digits 	=> 8 				 															);
				
				$summaries{$sample}{$key}{REP2_EXP}	= bk_basic_divide
																							( numerator						=> $summaries{$sample}{$key}{REP2_PER},
																								denominator					=> $library{$key}{LIB_AVG},
																								significant_digits 	=> 8 				 															);
				$summaries{$sample}{$key}{REP3_EXP}	= bk_basic_divide
																							( numerator						=> $summaries{$sample}{$key}{REP3_PER},
																								denominator					=> $library{$key}{LIB_AVG},
																								significant_digits 	=> 8 				 															);
				my (@reps)													= (	$summaries{$sample}{$key}{REP1_EXP},
																								$summaries{$sample}{$key}{REP2_EXP},
																								$summaries{$sample}{$key}{REP3_EXP});
				$summaries{$sample}{$key}{EXP_AVG}	= bk_basic_round
																							(	significant_digits 	=> 8,
																								number 							=> Math::NumberCruncher::Mean(\@reps)); #3.2.0 zero error fixed with 8 sig digits
				if 																		($summaries{$sample}{$key}{EXP_AVG} < $exp_thresh){ #3.2.0 expression threshold
					$summaries{$sample}{$key}{TTEST}	= "-";					
					$summaries{$sample}{$key}{CHI2}		= "-";					
					$summaries{$sample}{$key}{CLUS}		= "-";					
					$summaries{$sample}{$key}{LOG_CUMM}= "-";					
					$summaries{$sample}{$key}{LOG_CLUS}= "-";					
					$filter{$sample."_DR"}{$key}			= $summaries{$sample}{$key};
					delete 															$summaries{$sample}{$key}; 
					next;
				}
				$summaries{$sample}{$key}{EXP_STDEV}= bk_basic_round
																							(	significant_digits 	=> 4,
																								number 							=> Math::NumberCruncher::StandardDeviation(\@reps));
				my $ttest 													= new Statistics::TTest; 
				$ttest															-> set_significance(95);
				my (@x_val)													= ($library{$key}{REP1_PER},$library{$key}{REP2_PER},$library{$key}{REP3_PER});
				my (@y_val)													= ($summaries{$sample}{$key}{REP1_PER},$summaries{$sample}{$key}{REP2_PER},$summaries{$sample}{$key}{REP3_PER});
				$ttest															-> load_data(\@x_val,\@y_val); 
				$ttest															-> set_significance(95);
				my $out_p														= sprintf("%.10g", $ttest->{t_prob});
				$out_p															= bk_basic_round
																							(	significant_digits 	=> 4,
																								number 							=> $out_p);
				$summaries{$sample}{$key}{TTEST}		= $out_p;
				my $n1															= $library{$key}{CNT_AVG}; 
				my $n2															=	$summaries{$sample}{$key}{CNT_AVG}; 
				my $np1															= $lib_tot-$n1;	
				my $np2															= $sam_tot-$n2;		  
				my @obs 														= ([$n1, $n2], [$np1, $np2]);
				my $chi 														= new Statistics::ChisqIndep;
				$chi																-> load_data(\@obs);
				$summaries{$sample}{$key}{CHI2}			= bk_basic_round
																							(	significant_digits 	=> 4,
																								number 							=> sprintf("%.10g", $chi->{p_value}));
				if 																	($summaries{$sample}{$key}{CHI2} > 0.05 or $summaries{$sample}{$key}{TTEST} > 0.05){ 	
					$summaries{$sample}{$key}{CLUS}		= "-";					
					$summaries{$sample}{$key}{LOG_CUMM}= "-";					
					$summaries{$sample}{$key}{LOG_CLUS}= "-";					
					$filter{$sample."_STAT"}{$key}		= $summaries{$sample}{$key};
					delete 															$summaries{$sample}{$key}; 
					next;
				}			
				$summaries{$sample}{$key}{REP1_PER}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP1_PER}));
				$summaries{$sample}{$key}{REP2_PER}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP2_PER}));
				$summaries{$sample}{$key}{REP3_PER}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP3_PER}));
				$summaries{$sample}{$key}{REP1_EXP}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP1_EXP}));
				$summaries{$sample}{$key}{REP2_EXP}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP2_EXP}));
				$summaries{$sample}{$key}{REP3_EXP}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> sprintf("%.10g", $summaries{$sample}{$key}{REP3_EXP}));
			}
			$out_metric{$sample}{SIG}								= keys(%{$summaries{$sample}});
			$out_metric{$sample}{REM_RED}						= keys(%{$filter{$sample."_DR"}});
			$out_metric{$sample}{REM_STAT}					= keys(%{$filter{$sample."_STAT"}});
			$out_metric{$sample}{TOTAL}							= $out_metric{$sample}{SIG}+$out_metric{$sample}{REM_RED}+$out_metric{$sample}{REM_STAT};
		}
		foreach my $sample												(keys %summaries){		
			my (@out_summary)												= "OLIGO\tLIB_CNT\tLIB_PER\tREP1_CNT\tREP2_CNT\tREP3_CNT\tAVG_CNT\t".
																								"REP1_PER\tREP2_PER\tREP3_PER\tREP1_EXP\tREP2_EXP\tREP3_EXP\tEXP_AVG\tEXP_STDEV\t".
																								"TTEST\tCHI2\tPROT\tNUCL\n";
			my %tmp_hash														= %{$summaries{$sample}};
			foreach my $key 												( sort { $tmp_hash{$b}->{EXP_AVG} <=> $tmp_hash{$a}->{EXP_AVG} } keys %tmp_hash) {
				push																	@out_summary, "$key\t".
																														sprintf("%.0f", $summaries{$sample}{$key}{LIB_CNT})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{LIB_PER})."\t".
																														sprintf("%.0f", $summaries{$sample}{$key}{REP1_CNT})."\t".
																														sprintf("%.0f", $summaries{$sample}{$key}{REP2_CNT})."\t".
																														sprintf("%.0f", $summaries{$sample}{$key}{REP3_CNT})."\t".
																														sprintf("%.0f", $summaries{$sample}{$key}{CNT_AVG})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{REP1_PER})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{REP2_PER})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{REP3_PER})."\t".
																														sprintf("%.2f", $summaries{$sample}{$key}{REP1_EXP})."\t".
																														sprintf("%.2f", $summaries{$sample}{$key}{REP2_EXP})."\t".
																														sprintf("%.2f", $summaries{$sample}{$key}{REP3_EXP})."\t".
																														sprintf("%.2f", $summaries{$sample}{$key}{EXP_AVG})."\t".
																														sprintf("%.2f", $summaries{$sample}{$key}{EXP_STDEV})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{TTEST})."\t".
																														sprintf("%.4f", $summaries{$sample}{$key}{CHI2})."\t".
																																						$summaries{$sample}{$key}{PROT}."\t".
																																						$summaries{$sample}{$key}{NUCL}."\n";
			}
			write_file															($args{Path}."/_reports/".$sample.".sum", @out_summary);
		}
		##################CLUSTERING###################
		my (%hsh_align);
		my (@files)																= <$args{Path}/_array/*.fasta>;
		foreach my $link 													(@files){
			my %tmp_hsh 														= bk_biology_fastaToHash($link);
			my $tmp_name														= bk_basic_removeFileType(bk_basic_returnFile($link));
			$hsh_align{$tmp_name}										= \%tmp_hsh;
		}		
		my (%clusters);
		foreach my $sample												(keys %summaries){	
			my (%position);
			foreach my $link 													(@files){
				my (@file)															= read_file($link);
				my ($file_name)													= bk_basic_removeFileType(bk_basic_returnFile($link));
				chomp																		$file[1];
				my (@tmp_pos)														= (0) x length($file[1]);
				$position{uc($file_name)}								= [@tmp_pos];
			}		
			my %tmp_hash														= %{$summaries{$sample}};		
			foreach my $key 												( keys %tmp_hash) {
				my (@line_loc)												= split /\:/, $key;
				$line_loc[0]													=~ s/\>//g;
				for 																	(my $x = ($line_loc[2]-1); $x < $line_loc[3]; $x++) {
																							$position{uc($line_loc[0])}[$x]++;
				}
				$tmp_hash{$key}{REPEAT}								=0;
				$tmp_hash{$key}{ALL_REP}							="";
			}
			##### Noise Removal #####
			foreach my $pos													( keys %position){
				for																		(my $i = (@{$position{$pos}})-1; $i >= 0; $i--){
					if 																	($position{$pos}[$i] >= 3){}else{ $position{$pos}[$i] = 0; }
				}
			}
			#####
			my $total_clus													= 1;
			foreach my $key 												( keys %position) {
				my ($flg_print, $cnt_clus, $start_pos)= (0) x 3;
				for																		(my $i = 0; $i < @{$position{$key}}; $i++){
					if																	($flg_print eq 1 and $position{$key}[$i] eq 0){
						my $clus_name											= sprintf("%08d", $total_clus).":".$sample.":".$key.":".$start_pos.":".$i.":".$cnt_clus;
						$clusters{$sample}{$clus_name}		{START} 		= $start_pos;
						$clusters{$sample}{$clus_name}		{END} 			= $i;
						$clusters{$sample}{$clus_name}		{ALIGN} 		= $key;
						$clusters{$sample}{$clus_name}		{LIB_CNT} 	= 0;
						$clusters{$sample}{$clus_name}		{TEST_CNT}	= 0;
						$clusters{$sample}{$clus_name}		{EXP_CUMM} 	= 0;
						$clusters{$sample}{$clus_name}		{EXP_CLUS} 	= 0;
						$clusters{$sample}{$clus_name}		{CLUS} 			= 0;
						$clusters{$sample}{$clus_name}		{LOG_CUMM} 	= 0;
						$clusters{$sample}{$clus_name}		{LOG_CLUS} 	= 0;
						$clusters{$sample}{$clus_name}		{SCORE} 		= 0;
						$clusters{$sample}{$clus_name}		{LOGO} 			= (["Oligo","EXP","ALIGNMENT"]);
						$clusters{$sample}{$clus_name}		{IDS} 			= "";
						$clusters{$sample}{$clus_name}		{MAX_DEPTH} = $cnt_clus;
						$clusters{$sample}{$clus_name}		{CLUS_NUM}	= $total_clus;
						$clusters{$sample}{$clus_name}		{REP_CLUS}	= "";
						$total_clus++;
						$flg_print												= 0;                                        
						$cnt_clus													= 0;
						$start_pos												= 0;
					}else{
						if 																($position{$key}[$i] > 0 and $flg_print eq 0){ 
							$flg_print 											= 1;
							$start_pos											= $i+1; 
							$cnt_clus												= $position{$key}[$i];
						}else{
							if 															($position{$key}[$i] > $cnt_clus){ $cnt_clus = $position{$key}[$i]; }
						}
					}
				}
			}
			foreach my $key 												( keys %tmp_hash) {
				my (@line_loc)												= split /\:/, $key;
				$line_loc[0]													=~ s/\>//g;
				my (@reps)														= (	$summaries{$sample}{$key}{REP1_PER},
																									$summaries{$sample}{$key}{REP2_PER},
																									$summaries{$sample}{$key}{REP3_PER});
				$summaries{$sample}{$key}{TEST_PER}		= bk_basic_round
																								(	significant_digits 	=> 4,
																									number 							=> Math::NumberCruncher::Mean(\@reps));
				foreach my $clus											( keys %{$clusters{$sample}}){
					if																	(uc($line_loc[0]) eq uc($clusters{$sample}{$clus}{ALIGN})){
					##Include all hits within range {START}-{END}
						if																(($line_loc[2] ~~ [$clusters{$sample}{$clus}{START}..$clusters{$sample}{$clus}{END}] or $line_loc[3] ~~ [$clusters{$sample}{$clus}{START}..$clusters{$sample}{$clus}{END}]) or 
																								($clusters{$sample}{$clus}{START} ~~ [$line_loc[2]..$line_loc[3]])){
					######
							$clusters{$sample}{$clus}				{CLUS}++;
							$clusters{$sample}{$clus}				{LIB_CNT} 	= $clusters{$sample}{$clus}{LIB_CNT} + $summaries{$sample}{$key}{LIB_CNT};
							$clusters{$sample}{$clus}				{TEST_CNT} 	= $clusters{$sample}{$clus}{TEST_CNT} + $summaries{$sample}{$key}{CNT_AVG};
							$clusters{$sample}{$clus}				{EXP_CUMM} 	= $clusters{$sample}{$clus}{EXP_CUMM} + sprintf("%.1f", $summaries{$sample}{$key}{EXP_AVG});
							$clusters{$sample}{$clus}				{IDS}				= $clusters{$sample}{$clus}{IDS}.$key.":".sprintf("%.1f", $summaries{$sample}{$key}{EXP_AVG}).",";
							$tmp_hash{$key}	{REPEAT}++;
							$tmp_hash{$key}{ALL_REP}										= $tmp_hash{$key}{ALL_REP}.$clusters{$sample}{$clus}{CLUS_NUM}."&";
						}
					}
				}		
			}
		my @repeatedSegments;
		my @allRepeats;
		foreach my $key														( keys %tmp_hash){
				if 																		($tmp_hash{$key}{REPEAT} >1){
						push															@repeatedSegments, $key;
						push 															@allRepeats ,$tmp_hash{$key}{ALL_REP};
				}		
		}

		for																				(my $i = 0; $i<@repeatedSegments;$i++){
				foreach my $clus											( keys %{$clusters{$sample}}){
				if 																		($clusters{$sample}{$clus}{IDS} =~ /$repeatedSegments[$i]/){
						my (@tmp_keys)										= split /\,/, $clusters{$sample}{$clus}{IDS};
						my (@tmp_rep)											= split /\,/,$clusters{$sample}{$clus}{REP_CLUS};
						$clusters{$sample}{$clus}{REP_CLUS}="";
						for 															(my $tmpID = 0; $tmpID < @tmp_keys; $tmpID++){					
									if ($tmp_keys[$tmpID] 				=~ /$repeatedSegments[$i]/){
																								$clusters{$sample}{$clus}{REP_CLUS} = $clusters{$sample}{$clus}{REP_CLUS}.$allRepeats[$i].",";}
									else{
												if($tmp_rep[$tmpID]){
												if($tmp_rep[$tmpID] ne ""){
												$clusters{$sample}{$clus}{REP_CLUS} = $clusters{$sample}{$clus}{REP_CLUS}.$tmp_rep[$tmpID].",";}}
												else{
												$clusters{$sample}{$clus}{REP_CLUS} = $clusters{$sample}{$clus}{REP_CLUS}.",";}}
						}
						$clusters{$sample}{$clus}{REP_CLUS} =~ s/$clusters{$sample}{$clus}{CLUS_NUM}&//g;					
				}
		}
}
			
##########


			foreach my $clus												( keys %{$clusters{$sample}} ) {
				$clusters{$sample}{$clus}{EXP_CLUS}		=	bk_basic_divide
																								( numerator						=> ($clusters{$sample}{$clus}{TEST_CNT}+1),
																									denominator					=> ($clusters{$sample}{$clus}{LIB_CNT}+1),
																									significant_digits 	=> 4 				 											  				);
				$clusters{$sample}{$clus}{LOG_CUMM}		= bk_basic_round
																							(	significant_digits 	=> 2,
																								number 							=> bk_basic_log10(($clusters{$sample}{$clus}{EXP_CUMM}/$clusters{$sample}{$clus}{CLUS})));
				$clusters{$sample}{$clus}{LOG_CLUS}		= bk_basic_round
																							(	significant_digits 	=> 2,
																								number 							=> bk_basic_log10($clusters{$sample}{$clus}{EXP_CLUS}));
				my ($find_consensus)									= bk_biology_phageAlignment							
																								( IDs					=>	$clusters{$sample}{$clus}{IDS},
																								  REP					=>	$clusters{$sample}{$clus}{REP_CLUS},#flagging
																									Sample 			=>	$sample,
																									Cluster 		=>	$clus,													
																									Path				=>	$args{Path}, 											
																									Alignments	=>	\%alignments										);
				$clusters{$sample}{$clus}{SCORE}			= $clusters{$sample}{$clus}{LOG_CUMM} + $clusters{$sample}{$clus}{LOG_CLUS} + $clusters{$sample}{$clus}{MAX_DEPTH};
				if																		($clusters{$sample}{$clus}{LIB_CNT} < 10){ $clusters{$sample}{$clus}{SCORE} -= 1; }
				if																		($clusters{$sample}{$clus}{TEST_CNT} < 10){ $clusters{$sample}{$clus}{SCORE} -= 1; }
				if																		($find_consensus ne 0){ $clusters{$sample}{$clus}{SCORE} -= 2; }
				# $clusters{$sample}{$clus}{LOG_CUMM}			=	bk_basic_divide
																								# ( numerator						=> $clusters{$sample}{$clus}{CLUS},
																									# denominator					=> $library{$clus}{CLUS},
																									# significant_digits 	=> 4 				 											  );
				# $clusters{$sample}{$clus}{LOG_CLUS}			= bk_basic_round
																								# (	significant_digits 	=> 4,
																									# number 							=> sprintf("%.10g", $clusters{$sample}{$clus}{EXP_CUMM}*$clusters{$sample}{$clus}{LOG_CUMM}));
				# $clusters{$sample}{$clus}{SCORE}			= bk_basic_round
																								# (	significant_digits 	=> 4,
																									# number 							=> sprintf("%.10g", $clusters{$sample}{$clus}{EXP_CLUS}*$clusters{$sample}{$clus}{LOG_CUMM}));																							
			}

			$out_metric{$sample}{CLUS}								= keys(%{$clusters{$sample}});			
		}	
		foreach my $sample												(keys %out_metric){
			print $out_met													"$sample\t$out_metric{$sample}{TOTAL}\t$out_metric{$sample}{SIG}\t$out_metric{$sample}{REM_RED}\t$out_metric{$sample}{REM_STAT}\t$out_metric{$sample}{CLUS}\n";
		}
		foreach my $sample												(keys %clusters){		
			my (@out_summary)												= "CLUSTER\tALIGN\tLIB_CNT\tTEST_CNT\tEXP_CUMM\tEXP_CLUS\tCLUS\tLOG_CUMM\tLOG_CLUS\tSCORE\n";
			my %tmp_hash														= %{$clusters{$sample}};
			foreach my $key 												( sort { $tmp_hash{$b}->{SCORE} <=> $tmp_hash{$a}->{SCORE} } keys %tmp_hash) {
				push																	@out_summary, "$key\t".
																																						$clusters{$sample}{$key}{ALIGN}."\t".
																														sprintf("%.0f", $clusters{$sample}{$key}{LIB_CNT})."\t".
																														sprintf("%.0f", $clusters{$sample}{$key}{TEST_CNT})."\t".
																														sprintf("%.2f", $clusters{$sample}{$key}{EXP_CUMM})."\t".
																														sprintf("%.2f", $clusters{$sample}{$key}{EXP_CLUS})."\t".
																														sprintf("%.0f", $clusters{$sample}{$key}{CLUS})."\t".
																														sprintf("%.2f", $clusters{$sample}{$key}{LOG_CUMM})."\t".
																														sprintf("%.2f", $clusters{$sample}{$key}{LOG_CLUS})."\t".
																														sprintf("%.2f", $clusters{$sample}{$key}{SCORE})."\n";
			}
			write_file															($args{Path}."/_reports/".$sample.".fin", @out_summary);
		}
		foreach my $key														(keys %library){
			my (@tmp_values)												= split /\:/, $key;
					$tmp_values[0]											=~ s/\>//g;
			for																			(my $i=$tmp_values[2]; $i<$tmp_values[3]; $i++){
				$align_histo{uc($tmp_values[0])}			{LIB}{UNI}[$i]++;
				if																		($library{$key}{CNT_AVG} >= 1){
					$align_histo{uc($tmp_values[0])}		{LIB}{ACT}[$i]++;
					$align_histo{uc($tmp_values[0])}		{LIB}{EXP}[$i] = $align_histo{uc($tmp_values[0])}{LIB}{EXP}[$i]+$library{$key}{CNT_AVG};
				}
			}
		}
		foreach my $sample												(keys %clusters){
			foreach my $key													(keys %{$clusters{$sample}}){
				my (@tmp_values)											= split /\:/, $key;
						$tmp_values[0]										=~ s/\>//g;
				for																		(my $i=$tmp_values[3]; $i<$tmp_values[4]; $i++){
					if																	($clusters{$sample}{$key}{TEST_CNT} >= 1){
					$align_histo{uc($tmp_values[2])}	{$sample}{EXP}[$i] = $align_histo{uc($tmp_values[2])}{$sample}{EXP}[$i]+$clusters{$sample}{$key}{TEST_CNT};
					}
				}				
			}
		}
		foreach my $align													(keys %align_histo){
			foreach my $sample											(keys %clusters){
				my (@data);
						my $length												= @{$align_histo{$align}{LIB}{UNI}};
						for																(my $i=0; $i<$length; $i++) {	
							$data[0][$i]										= $i;
							$data[1][$i]										= $align_histo{$align}{LIB}{UNI}[$i];
							$data[2][$i]										= $align_histo{$align}{LIB}{ACT}[$i];
							$data[3][$i]										= $align_histo{$align}{$sample}{ACT}[$i];
						}
				my $mygraph 													= GD::Graph::lines->new(800, 300);
						$mygraph													-> set(
																									x_label		 				=> 'Position',
																									y_label						=> 'Count',
																									title			 				=> $align.":".$sample." Diversity Coverage Map",
																									line_types				=> [1, 1, 1],
																									line_width				=> 1,
																									x_tick_number 		=> 10, 
																									x_label_skip 			=> 2, 
																									x_label_position 	=> 1/2, 
																									y_tick_number 		=> 10, 
																									y_label_skip 			=> 2, 
																									y_min_value 			=> 0, 
																									y_max_value 			=> max(@{$data[1]}),  
																									x_max_value 			=> $length,  
																									dclrs			 				=> ['blue', 'red', 'orange'],
																									box_axis 					=> 0, 
																									r_margin 					=> 15, 
																									transparent 			=> 0, 														) or 
																									warn $mygraph->error;
						$mygraph													-> set_legend_font(GD::gdMediumBoldFont);
						$mygraph													-> set_legend('LIB', 'LIB-ACTUAL',$sample.'-ACTUAL');
						my $myimage 											= $mygraph -> plot(\@data) or die $mygraph->error;
						my ($out_name)										= '>'.$args{Path}."/_images/".$sample."_".$align.'_lib.png';
						open															(IMG, $out_name) or die $!;
						binmode 													IMG;
						print 														IMG $myimage -> png;	
				my (@data1);
						my $length1												= @{$align_histo{$align}{LIB}{UNI}};
						for																(my $i=0; $i<$length1; $i++) {	
							$data1[0][$i]										= $i;
							$data1[1][$i]										= $align_histo{$align}{LIB}{EXP}[$i];
							$data1[2][$i]										= $align_histo{$align}{$sample}{EXP}[$i];
						}
				my $mygraph1													= GD::Graph::lines->new(800, 300);
						$mygraph1													-> set(
																									x_label		 				=> 'Position',
																									y_label						=> 'Count',
																									title			 				=> $align.":".$sample." Expression Coverage Map",
																									line_types				=> [1, 1],
																									line_width				=> 1,
																									x_tick_number 		=> 10, 
																									x_label_skip 			=> 2, 
																									x_label_position 	=> 1/2, 
																									x_labels_vertical => 1, 
																									y_tick_number 		=> 10, 
																									y_label_skip 			=> 2, 
																									y_min_value 			=> 0, 
																									y_max_value 			=> max(@{$data1[1]},@{$data1[2]}), 
																									x_max_value 			=> $length, 
																									dclrs			 				=> ['blue', 'red'],
																									box_axis 					=> 0, 
																									r_margin 					=> 15, 
																									transparent 			=> 0, 														) or 
																									warn $mygraph1->error;
						$mygraph1													-> set_legend_font(GD::gdMediumBoldFont);
						$mygraph1													-> set_legend('LIB-EXPRESSION',$sample.'-EXPRESSION');
						my $myimage1 											= $mygraph1 -> plot(\@data1) or die $mygraph1->error;
						my ($out_name1)										= '>'.$args{Path}."/_images/".$sample."_".$align.'_exp.png';
						open															(IMG, $out_name1) or die $!;
						binmode 													IMG;
						print 														IMG $myimage1 -> png;								
			}
		}
	}	
	
	sub bk_biology_FastaAlignToHash {					#METHOD: B<bk_biology_FastaAlignToHash> #needs ARGS
			my (%args) 															= (@_);

			my (%all);
			my (@align_files)												= <$args{Path}/_array/*.fasta>;
			foreach my $file												(@align_files){
				my (@in_align)												= read_file($file);
						until															($in_align[@in_align-1] ne ""){ pop @in_align; }
						foreach														(@in_align){ chomp; }
						my (@names)												= grep { $_ =~ /\>/ } @in_align;
						my (@sequence)										= grep { $_ !~ /\>/ } @in_align;
						my ($align_name)									= bk_basic_removeFileType(bk_basic_returnFile($file));
				for																		(my $i=0; $i<@names; $i++){
					$names[$i]													=~ s/\>//;
					$names[$i]													=~ s/\_| /\-/g;
					$all{$align_name}{uc($names[$i])}		= $sequence[$i];	
				}
			}
			return %all;
	}
	
	sub bk_biology_phageAlignment {						#METHOD: B<bk_biology_phageAlignment> #needs ARGS
		my (%args) 															= (@_);
		
		my (%tmp_align)													= %{$args{Alignments}};
		my ($cluster)														= $args{Cluster};
				$cluster														=~ s/\:| |\_/\-/g;
				$cluster														=~ s/\>//;
		my (@tmp_keys)													= split /\,/, $args{IDs};
		my (@tmp_reps)													= split /\,/, $args{REP};
		my (@key_parse);
		foreach my $key													(@tmp_keys){
			my (@tmp_key)													= split /\:/, $key;
					$tmp_key[0]												=~ s/\>//;
			push																	@key_parse, [@tmp_key];
		}
		my ($start)															= min(map { $_->[2] } @key_parse);
				$start--;
		my ($stop)															= max(map { $_->[3] } @key_parse); 
		my (@align_out);
		for																			(my $i=0; $i<@key_parse; $i++){
			$key_parse[$i][1]											=~ s/\_| /\-/g;
			my (@tmp_sequence)										= split //, substr($tmp_align{$key_parse[$i][0]}{uc($key_parse[$i][1])}, $start, ($stop-$start));
			my ($gap_start)												= $key_parse[$i][2]-($start+1); 
			my ($gap_end)													= (@tmp_sequence-($stop-$key_parse[$i][3]))-1;
			for																		(my $y=0; $y<@tmp_sequence; $y++){
				if																	($y < $gap_start or $y > $gap_end){ 
					$tmp_sequence[$y]									= " ";
				}
			}
			my (@tmp_name)												= @{$key_parse[$i]};
			my ($out_exp)													= $tmp_name[@tmp_name-1];
			pop 																	@tmp_name;
			my ($out_name)												= join "_", @tmp_name;
			my $out_rep														= $tmp_reps[$i];
			if($out_rep){
			$out_rep															=~ s/^(.*)&/$1/;
			$out_rep														  =~ s/&/,/g;}
			else{$out_rep = " ";}
			push																	@align_out, [$out_name,$out_exp,$out_rep,@tmp_sequence,$gap_start];
		}
		my (@con_seq)														= (" ") x (@{$align_out[0]}-1);
		push 																		@align_out, ["Consensus"," ",@con_seq];
		my (@align_rev);
		my ($ret_con)														= 0;
		for																			(my $i=0; $i<@align_out; $i++){
			push																	@{$align_rev[$i]}, $align_out[$i][0];
			push																	@{$align_rev[$i]}, $align_out[$i][1];
			push 																	@{$align_rev[$i]}, $align_out[$i][2];
		}
		for																			(my $i=3; $i<@{$align_out[0]}; $i++){
			my (@bases)														= map { $_ -> [$i] } @align_out;
			pop																		@bases;
			my (@gap)															= grep { $_ =~ /\-/ } @bases;
			my (@tot_aa)													= grep { $_ !~ /\-/	} @bases;
			my $score															= 0;
			if 																		(@gap > 0 and @bases > 0){ $score = @gap/@bases; }
			if																		($score ne 1){
				for																	(my $y=0; $y<@align_out-1; $y++) { 
					push															@{$align_rev[$y]}, $align_out[$y][$i];
				}
				my (@uni_base)											= bk_basic_uniqueArray(@bases);
				my ($con_amino)											= " ";
				foreach my $aa											(@uni_base){
					if 																($aa ne "-" and $aa ne "*"){
						my (@aminos)										= grep { $_ =~ /$aa/ } @bases;
						my $aa_score										= 0;
						if 															(@aminos > 0 and @tot_aa > 2){ $aa_score = @aminos/@tot_aa; }
						if															($aa_score > 0.59){
							$con_amino										= $aa;
							$ret_con											= -1;
							if														(@tmp_keys eq 1){ $con_amino = " "; }
						}
					}
				}
				push																@{$align_rev[@align_rev-1]}, $con_amino;
			}
		}
		my (@out_file);
		foreach																	(@align_rev){
			my $start															= pop @{$_};
			my $out_seq														= join "", @{$_}[3..(@{$_}-1)]; 
			$out_seq															=~ s/\s+$//; 
			my $virus = "";
			my @splitUp														= split /_/, @{$_}[0];
			if																		($splitUp[1]){	
																						$splitUp[1]	=~ s/NC-/NC/g; 
																						$splitUp[1]	=~ s/YP-/YP/g; 
				my @NameArray												= split /-/, $splitUp[1];
				if 																	($NameArray[1])	{$virus = 	$NameArray[1];}} 
				##########

			push																	@out_file, ["@{$_}[0]\t@{$_}[1]\t@{$_}[2]\t$out_seq", $start, $virus];
		}
		my $consensus														= pop @out_file;
		@out_file																= sort { $a->[2] cmp $b->[2] || $a->[1] <=> $b->[1]} @out_file;
		my @out_file1														= ($cluster."\nOLIGO\tEXP\tREPEAT\tALIGNMENT\n");
		foreach																	(@out_file){
			push																	@out_file1, @{$_}[0]."\n";
		}
		push																		@out_file1, @{$consensus}[0];
		$cluster																=~ s/\:| |\_/\-/g;
		write_file															($args{Path}."/_sequence/".$args{Sample}."_".$cluster."_cluster.aln", @out_file1);
		
		return $ret_con;
	}
}

1;