#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

use Getopt::Long;														#USE:
use Pod::Usage;															#USE:
use Cwd;																		#USE:
use Parallel::ForkManager;									#USE:
use File::Slurp;														#USE:
use File::Copy;															#USE:

																						#B<SCRIPT INFO: "DPA_Analysis", "3.2.1", "20180611"> 
my ($var_nam, $var_rev, $var_upd)						= ("DPA_Analysis", "3.2.1", "20180611");

my ($var_hlp, $var_man, $var_ver, $var_dev) = (0) x 4;
my ($var_are, $var_lib, $var_cus)						= (0) x 3;
my ($var_sil, $var_inp,	$var_dup) 					= (0) x 3;
my ($var_exp)																= (1);
my ($var_cor, $var_pat, $var_tim) 					= (20, getcwd(), time);
my ($var_sli, $var_osi, $var_arg) 					= (7, 36, "");#default 36 arrray
my ($var_exc) 															= ("");
my ($var_vir)																= "EBOV";

my (@args_out)															= (@ARGV);

my (@arr_sam, @arr_ftq, @arr_ind);

my ($out_log, $out_met, $log_nam, $out_seq);
my ($out_ind, $out_dup, $out_cor, $out_cle);
my ($out_exp, $out_lib);

unless ($var_dev eq 1){											#B<SECTION: Logging and Script Operation>
	GetOptions																('help|?' => \$var_hlp, man			=> \$var_man,
																						"ver"			=> \$var_ver, "sil"		=> \$var_sil,
																						"pat=s"		=> \$var_pat, "cor=i"	=> \$var_cor,
																						"inp=i"		=> \$var_inp, "are"		=> \$var_are,
																						"lib"			=> \$var_lib, "cus"		=> \$var_cus,
																						"exp=i"		=> \$var_exp, "dup"		=> \$var_dup,
																						"sli=i"		=> \$var_sli, "osi=i"	=> \$var_osi,
																						"exc=s"		=> \$var_exc, "vir=s"		=> \$var_vir,
																						) or pod2usage(2);
	pod2usage																	(1) if $var_hlp;
	pod2usage																	(-exitstatus => 1, -verbose => 2) if $var_man;
	unless($var_ver eq 0){										print PHAGE_version(); exit; } 
	require 																	"PATH/bk_basic.pl"; #replace with direct path to script
	require 																	"PATH/bk_biology.pl"; #replace with path to script
	require 																	"PATH/bk_wrappers.pl"; #replace with path to script
	$var_arg																	= bk_basic_argCapture(@args_out);
	$var_tim																	= bk_basic_time($var_tim);
	open $out_log,														">", "$var_pat/$var_nam-$var_tim.log" 
			or die																"Could not open $var_pat/$var_nam-$var_tim.log - $!";
	bk_basic_operationLog											( Path			=> $var_pat,
																							Time			=> $var_tim,
																							Name			=> $var_nam,
																							ProcessID => $$,
																							Args			=> $var_arg,
																							Version		=> PHAGE_version(),
																							ThirdPary => "BLAST",		
																							File			=>	$out_log							);
	$log_nam																	= "cmdOUT-$var_tim.logf";
	unless($var_sil eq 1){										bk_basic_errorLog
																							( Path => $var_pat,
																								File => $log_nam	); 
	}
}

my $fork_manager														= new Parallel::ForkManager($var_cor);

unless ($var_dev eq 1){											#B<SECTION: File Structure and Validation> 
	bk_basic_directoryStructure								( Path => $var_pat, 
																							Folders => ["_images","_logs","_library",
																													"_metrics","_reports","_sequence"]	);
																						`cp -Rf //data/_scripts/jrk/_resources/_html/_images $var_pat/_reports`;
	print $out_log														bk_basic_writeLog("Completed file structure and file check.");
}

unless ($var_dev eq 1) {										#B<SECTION: Initialize globals.> 
	@arr_sam																	= <$var_pat/*/>;
			foreach my $folder (@arr_sam) 				{ $folder = bk_basic_returnFile($folder); }
			@arr_sam															= grep { substr($_,0,1) !~/_/ } @arr_sam;
	@arr_ftq																	= <$var_pat/*/*_R*fastq.gz>;
	@arr_ind																	= <$var_pat/*/*_I*fastq.gz>;
	if																				(@arr_ftq eq 0){ print "Failed to detect sequence files."; }
	foreach my $sample 												(@arr_sam){
		my @tmp_ftq															= grep { $_ =~ /\/$sample\// } @arr_ftq;
		my @tmp_ind															= grep { $_ =~ /\/$sample\// } @arr_ind;
		if																			(@arr_ftq < 2){ print "Failed to detect sequence files for $sample.\n"; }
		if																			(@arr_ind eq 2){ print "Dual index files detected for $sample.\n"; }
		elsif																		(@arr_ind eq 1){ print "Single index file detected for $sample.\n"; }
		elsif																		(@arr_ind eq 0){ print "No index files detected for $sample. Turning off q30 for this sample.\n"; }
	}
	print 																		"------\n";
	print $out_log														bk_basic_writeLog("Program Initialized."); 
}

unless ($var_are eq 1){											#B<SECTION: Array Description.> 
	print $out_log														bk_basic_writeLog("Starting Array Creation."); 
	my ($tmp_time)														= (time);
	open $out_met,														">", "$var_pat/_metrics/06_Array.met" 
		or die																	"Could not open 06_Array.met - $!";		
	unless 																		($var_cus eq 1){
		my @files 															= <$var_pat/_array/*.fasta>;
				bk_biology_PhageArray								(	Path				=>	$var_pat,
																							Core				=>	$var_cor,
																							Window			=>	$var_osi,
																							Slide				=>	$var_sli,
																							Metrics			=>	$out_met,
																							Exclude			=>  $var_exc,
																							Alignments 	=>	\@files		);
	}else{
		my @files 															= <$var_pat/_array/*.oli>;
	}	
	close 																		$out_met;
	my $tmp_time_out													= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
	print $out_log														bk_basic_writeLog("Array created in ".$tmp_time_out." mins.");	
}
if(-e "$var_pat/_metrics/06_Array.met"){		print read_file($var_pat."/_metrics/06_Array.met"); }
print $out_log															bk_basic_writeLog("-----------------Array Description Completed--------------------");

unless ($var_lib eq 1){											#B<SECTION: Expression.> 
	my $tmp_time_out													= 0;
	print $out_log														bk_basic_writeLog("Starting Library Evaluation."); 
	my (@arr_fasta)														= <$var_pat/_sequence/*.fasta*>;
			foreach 															(@arr_fasta) { unlink $_; }
	my ($tmp_time)														= (time);
	open $out_exp,														">", "$var_pat/_metrics/01_Expression.met" 
		or die																	"Could not open 01_Expression.met - $!";	
	open $out_lib,														">", "$var_pat/_metrics/02_Hits.met" 
		or die																	"Could not open 02_Hits.met - $!";	
	open $out_cor,														">", "$var_pat/_metrics/04_Correction.met" 
		or die																	"Could not open 03_Correction.met - $!";	
	open $out_cle,														">", "$var_pat/_metrics/03_Cleaning.met" 
		or die																	"Could not open 04_Cleaning.met - $!";	
	open $out_dup,														">", "$var_pat/_metrics/05_Duplicates.met" 
		or die																	"Could not open 05_Duplicates.met - $!";	
	open $out_met,														">", "$var_pat/_metrics/seq.met" 
		or die																	"Could not open seq.met - $!";	
	my @gunzip_all														= <$var_pat/*/*fastq.gz>;
			bk_basic_paraGunzip										( Core	=> $var_cor, 
																							Files => \@gunzip_all	);	
			bk_biology_randInp										( Metrics_File	=> $out_met,
																							Path					=> $var_pat,
																							Core					=> $var_cor,
																							Files					=> \@arr_sam,
																							Input_Number	=> $var_inp		);	
	my @rem_uz																= <$var_pat/*/*uz>;
			foreach my $file											(@rem_uz){ unlink $file; }
	$tmp_time_out															= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
	print $out_log														bk_basic_writeLog("Files unzipped and inputs created in ".$tmp_time_out." mins.");
	$tmp_time																	= time;
	unless($var_dup eq 1){
		bk_biology_removeDuplicates							( Metrics_File		=> $out_dup,
																							Path						=> $var_pat,
																							Core						=> $var_cor,	
																							Method					=> "",
																							Samples					=> \@arr_sam );
	print $out_log														bk_basic_writeLog("Duplicate Removal Complete."); 																							
		$tmp_time_out														= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
		print $out_log													bk_basic_writeLog("Duplicates removed in ".$tmp_time_out." mins.");
	}
	$tmp_time																	= time;
	print $out_cle														"CLEANING\nSAMPLE\tTOTAL\tINSERT\t%TOT\tNO INSERT\t%TOT\tCHIMERIC\t%TOT\tREPLICATE\t%TOT\tBAD PRIM\t%TOT\n";
	print $out_cor														"CORRECTION\nSAMPLE\tCORRECTED\tNOMATCH\tINDEL\tAVG ERR\tDEV\tAVG LEN\tDEV\n";
	print $out_lib														"EXPRESSION\nSAMPLE\tGOOD\t%TOTG\tSHORT\t%TOT\tHIT\t%TOT\tTOTAL\n";
	print $out_exp														"CLUSTERING\nSAMPLE\tHITS\tSIG\tREM_REDUC\tREM_STAT\tCLUSTERS\n";
	foreach my $sample												(@arr_sam){
		my (@reads)															= <$var_pat/$sample/*_R*inp>;
		bk_biology_readMateCorrection						( Core					=> $var_cor,
																							Path					=> $var_pat,
																							Sample				=> $sample,
																							Metrics				=> $out_cor,
																							Reads					=> \@reads 		);
		$tmp_time_out														= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
		print $out_log													bk_basic_writeLog("$sample read mate correction completed in ".$tmp_time_out." mins.");		
		$tmp_time																= time;
		bk_biology_PhageCleaning								( Core					=> $var_cor,
																							Path					=> $var_pat,
																							Sample				=> $sample,
																							Metrics				=> $out_cle 	);
		$tmp_time_out														= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
		print $out_log													bk_basic_writeLog("$sample phage cleaning completed in ".$tmp_time_out." mins.");		
		$tmp_time																= time;
	}
	my (@arr_insert)													= <$var_pat/_sequence/*_insert*>;
	foreach 																	(@arr_insert){ $_ = bk_basic_removeFileType(bk_basic_returnFile($_)); }
	foreach my $sample												(@arr_insert){
		bk_biology_PhageLibrary									(	Core					=> $var_cor,
																							Path					=> $var_pat,
																							Sample				=> $sample,
																							Metrics				=> $out_lib		);
		$tmp_time_out														= bk_basic_divide
																							( numerator						=> (time-$tmp_time),
																								denominator					=> 60,
																								significant_digits 	=> 2 								);
		print $out_log													bk_basic_writeLog("$sample hit analysis completed in ".$tmp_time_out." mins.");		
		$tmp_time																= time;
	}
	bk_biology_PhageAnalysis									(	Core					=> $var_cor,
																							Path					=> $var_pat,
																							Sample				=> \@arr_sam,
																							Metrics				=> $out_exp,
																							Exp						=> $var_exp);
	my (@arr_txt)															= <$var_pat/*/*.txt*>;
			foreach 															(@arr_txt) { unlink $_; }
	my (@arr_ind)															= <$var_pat/*/*.ind*>;
			foreach 															(@arr_ind) { unlink $_; }
	my (@arr_inp)															= <$var_pat/*/*.inp*>;
			foreach 															(@arr_inp) { unlink $_; }

}
if(-e "$var_pat/_metrics/01_Expression.met"){	print read_file($var_pat."/_metrics/01_Expression.met"); }
print $out_log															bk_basic_writeLog("-----------------Library Evaluation Completed--------------------");

unless ($var_dev eq 1){											#B<SECTION: Reporting.> 
	unless($var_dev eq 1){
		my @sample_detail												= <$var_pat/_reports/*fin>;
		foreach my $file												(@sample_detail){
			my ($sample)													= bk_basic_removeFileType(bk_basic_returnFile($file));
			my (@in_fin)													= read_file($file);
					foreach														(@in_fin){ chomp; }
					my (@headers)											= split /\t/, $in_fin[0];
					shift															@in_fin;
			my ($out_html)												= "<a name='TOP'></a><table><tr><td align='left'><h1><b>Sample Cluster Enrichment</b></h1></td><td align='right'>".
																								"<a href='#Alignments'>Alignments</a> | <a href='#Diversity'>Diversity</a> | <a href='#Expression'>Expression</a>".
																								" | <a href='#Inserts'>Inserts</a> | <a href='#TOP'>Enrichment</a></td><tr><td colspan=2><hr><br></td></tr></table>";
					$out_html													= $out_html."<table><tr>";
					foreach my $value									(@headers){ $out_html = $out_html."<td class='header' align='center'>".$value."</td>";	}
					$out_html													= $out_html."</tr>";
			my ($odd)															= 0;
			foreach my $line											(@in_fin){
				my ($class)													= "data1";
						if															( 0 == $odd % 2 ) { $class = "data2"; }
				my (@values) 												= split /\t/, $line;
				$out_html														= $out_html."<tr>";
				$values[0]													=~ s/ |\:|\_/\-/g;
				$values[0]													=~ s/\>//g;
				$out_html														= $out_html."<td class='$class' align='left'><b><a href='#$values[0]'>$values[0]</a></b></td>";	
				shift																@values;
				foreach my $value										(@values){ $out_html = $out_html."<td class='$class' align='center'>".$value."</td>";	}
				$out_html														= $out_html."</tr>";
				$odd++;
			}
			$out_html															= $out_html."</table><br><br>";
			$out_html															= $out_html."<table><tr><td align='left'><a name='Alignments'></a><h1><b>Alignments</b></h1></td><td align='right'>".
																								"<a href='#Alignments'>Alignments</a> | <a href='#Diversity'>Diversity</a> | <a href='#Expression'>Expression</a>".
																								" | <a href='#Inserts'>Inserts</a> | <a href='#TOP'>Enrichment</a></td><tr><td colspan=2><hr><br></td></tr></table>";		
			my @alignments												= <$var_pat/_sequence/$sample*aln>;
			my @all_align;
			foreach my $file											(@alignments){
				my (@in_fin)												= read_file($file);
						foreach													(@in_fin){ chomp; }
						push														(@all_align, @in_fin);
						my ($align)											= shift @in_fin;
						my (@headers)										= split /\t/, $in_fin[0];
						shift														@in_fin;
						$out_html												= $out_html."<table><tr><td align='left'><a name='$align'></a><h1><b><br>CLUSTER: $align</b></h1>".
																								"</td><td align='right'><a href='#TOP'>TOP</a></td></tr></table><table><tr>";
						foreach my $value								(@headers){ $out_html = $out_html."<td class='header' align='center'>".$value."</td>";	}
						$out_html												= $out_html."</tr>";
				my ($odd)														= 0;
				foreach my $line										(@in_fin){
					my ($class)												= "data3";
							if														( 0 == $odd % 2 ) { $class = "data4"; }
					my (@values) 											= split /\t/, $line;
					$values[2]												=~ s/ /\&nbsp\;/g;
					$out_html													= $out_html."<tr><td class='$class' align='left' width='250'><b>$values[0]</b></td><td class='$class' width='60' align='center'>$values[1]</td><td class='$class' align='center'>$values[2]</td></td><td class='$class' align='center'>$values[3]</td></tr>"; 
					$odd++;
				}
			}	
			foreach 															(@all_align){ $_ = $_."\n"; }
			write_file														($var_pat."/_reports/$sample.all.aln",@all_align);
			$out_html															= $out_html."</table><br>";	
			$out_html															= $out_html."<table><tr><td align='left'><a name='Diversity'></a><h1><b>Diversity</b></h1></td><td align='right'>".
																								"<a href='#Alignments'>Alignments</a> | <a href='#Diversity'>Diversity</a> | <a href='#Expression'>Expression</a>".
																								" | <a href='#Inserts'>Inserts</a> | <a href='#TOP'>Enrichment</a></td><tr><td colspan=2><hr><br></td></tr></table>";		
			$out_html															= $out_html."<table><tr>";
			my @div																= <$var_pat/_images/$sample*_lib.png>;
			my ($o_div)														= 0;
			foreach my $file											(@div){
				my ($name)													= bk_basic_removeFileType(bk_basic_returnFile($file));
				if																	( 0 == $o_div % 2 ) { 
					$out_html													= $out_html."<tr><td colspan=2><br><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td>";
				}else{
					$out_html													= $out_html."<td colspan=2><br><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td></tr>";
				}
				$o_div++;
			}
			if																		( 0 == $o_div % 2 ) { 
				$out_html														= $out_html."<td colspan=2>&nbsp;</td></tr>";		
			}
			$out_html															= $out_html."</table><br>";	
			$out_html															= $out_html."<table><tr><td align='left'><a name='Expression'></a><h1><b>Expression</b></h1></td><td align='right'>".
																								"<a href='#Alignments'>Alignments</a> | <a href='#Diversity'>Diversity</a> | <a href='#Expression'>Expression</a>".
																								" | <a href='#Inserts'>Inserts</a> | <a href='#TOP'>Enrichment</a></td><tr><td colspan=2><hr><br></td></tr></table>";		
			$out_html															= $out_html."<table><tr>";
			my @exp																= <$var_pat/_images/$sample*_exp.png>;
			my ($o_exp)														= 0;
			foreach my $file											(@exp){
				my ($name)													= bk_basic_removeFileType(bk_basic_returnFile($file));
				if																	( 0 == $o_exp % 2 ) { 
					$out_html													= $out_html."<tr><td colspan=2><br><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td>";
				}else{
					$out_html													= $out_html."<td colspan=2><br><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td></tr>";
				}
				$o_exp++;
			}
			if																		( 0 != $o_exp % 2 ) { 
				$out_html														= $out_html."<td colspan=2>&nbsp;</td></tr>";		
			}
			$out_html															= $out_html."</table><br>";		
			$out_html															= $out_html."<table><tr><td align='left'><a name='Inserts'></a><h1><b>Inserts</b></h1></td><td align='right'>".
																								"<a href='#Alignments'>Alignments</a> | <a href='#Diversity'>Diversity</a> | <a href='#Expression'>Expression</a>".
																								" | <a href='#Inserts'>Inserts</a> | <a href='#TOP'>Enrichment</a></td><tr><td colspan=2><hr><br></td></tr></table>";		
			$out_html															= $out_html."<table><tr>";
			my (@insert)													= <$var_pat/_images/$sample*.png>;
					@insert														= grep { $_ !~ /\_exp|\_lib/ } @insert; 
			my ($o_ins)														= 0;
			foreach my $file											(@insert){
				my ($name)													= bk_basic_removeFileType(bk_basic_returnFile($file));
				if																	( 0 == $o_ins % 2 ) { 
					$out_html													= $out_html."<tr><td colspan=2><br><br><b><h1>$name</h1></b><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td>";
				}else{
					$out_html													= $out_html."<td colspan=2><br><br><b><h1>$name</h1></b><br><img src='../_images/$name.png' alt='$name' width='420' height='160'></td></tr>";
				}
				$o_ins++;
			}
			if																		( 0 != $o_ins % 2 ) { 
				$out_html														= $out_html."<td colspan=2>&nbsp;</td></tr>";		
			}
			$out_html															= $out_html."</table><br>";		
			my @html_log													= read_file("//data/_scripts/jrk/_resources/_html/template.html");
			foreach my $line											(@html_log){
				$line																=~ s/%WIDTH/870/;
				$line																=~ s/%PAGENAME/SAMPLE DETAIL: $sample/;
				$line																=~ s/%CONTENT/$out_html/;
			}
			write_file														($var_pat."/_reports/$sample.html",@html_log);
		}
		my (@align_excel)											= <$var_pat/_reports/*.aln>;
			foreach																(@align_excel){ 
				$_																= "_reports/".bk_basic_returnFile($_);
			}
		bk_basic_ExcelfromArray_Phage						( Worksheets 	=> \@align_excel,
																						Path 				=> $var_pat,
																						VIR					=> $var_vir,
																						Name				=> "projectSummary.xlsx"	);
	}

	unless ($var_dev eq 1){	
		my @sample_detail												= <$var_pat/_reports/*sum>;
		foreach my $file												(@sample_detail){
			my ($sample)													= bk_basic_removeFileType(bk_basic_returnFile($file));
			my (@in_fin)													= read_file($file);
					foreach														(@in_fin){ chomp; }
					my (@headers)											= split /\t/, $in_fin[0];
					pop 															@headers;
					pop 															@headers;
					shift															@in_fin;
			my ($out_html)												= "<a name='TOP'></a><table><tr><td align='left'><h1><b>Sample Replicate Summary</b></h1></td><td align='right'></td></tr></table>";
					$out_html													= $out_html."<table><tr>";
					foreach my $value									(@headers){ $out_html = $out_html."<td class='header' align='center'>".$value."</td>";	}
					$out_html													= $out_html."</tr>";
			my ($odd)															= 0;
			foreach my $line											(@in_fin){
				my ($class)													= "data1";
						if															( 0 == $odd % 2 ) { $class = "data2"; }
				my (@values) 												= split /\t/, $line;
						pop															@values;
						pop															@values;
				$out_html														= $out_html."<tr>";
				$values[0]													=~ s/ |\:|\_/\-/g;
				$values[0]													=~ s/\>//g;
				$out_html														= $out_html."<td class='$class' align='left'>$values[0]</td>";	
				shift																@values;
				foreach my $value										(@values){ $out_html = $out_html."<td class='$class' align='center'>".$value."</td>";	}
				$out_html														= $out_html."</tr>";
				$odd++;
			}
			$out_html															= $out_html."</table><br><br>";	
			my @html_log													= read_file("//data/_scripts/jrk/_resources/_html/template.html");
			foreach my $line											(@html_log){
				$line																=~ s/%WIDTH/870/;
				$line																=~ s/%PAGENAME/SAMPLE REPLICATE SUMMARY: $sample/;
				$line																=~ s/%CONTENT/$out_html/;
			}
			write_file														($var_pat."/_reports/$sample\_summary.html",@html_log);
		}
	}
	
	unless ($var_dev eq 1){	
		my @sample_detail												= <$var_pat/*PHAGE*>;
		my ($out_html)													= "";
		foreach my $file												(@sample_detail){
			my ($sample)													= bk_basic_removeFileType(bk_basic_returnFile($file));
			my (@in_fin)													= read_file($file);
					foreach														(@in_fin){ chomp; }
			$out_html															= $out_html."<a name='TOP'></a><table><tr><td align='left'><h1><b>$sample</b></h1></td><td align='right'></td></tr></table>";
			$out_html															= $out_html."<table>";
			foreach																(@in_fin){ 
				$out_html														= $out_html."<tr><td>$_</td></tr>";
			}
			$out_html															= $out_html."</table><br><br>";
		}
		my @html_log														= read_file("//data/_scripts/jrk/_resources/_html/template.html");
		foreach my $line												(@html_log){
			$line																	=~ s/%WIDTH/870/;
			$line																	=~ s/%PAGENAME/LOGS/;
			$line																	=~ s/%CONTENT/$out_html/;
		}
		write_file															($var_pat."/_reports/logs.html",@html_log);
	}
	
	unless ($var_dev eq 1){	
		if																			(-e "$var_pat/_metrics/seq.met"){
																						`mv $var_pat/_metrics/seq.met $var_pat/_metrics/07_Sequence.met`;
		}
		my @metrics															= <$var_pat/_metrics/*>;
		if																			($var_dup ne 0){ @metrics = grep { $_ !~ /05_Duplicates/ } @metrics; }
		my ($out_html)													= "<a name=TOP></a>";
			$out_html															= $out_html."";
		foreach my $file												(@metrics){
			my ($sample)													= bk_basic_removeFileType(bk_basic_returnFile($file));
			my (@in_fin)													= read_file($file);
					foreach														(@in_fin){ chomp; }
			$out_html															= $out_html."<table><tr><td align='left'><h1><b>$in_fin[0]</b><h1></td><td align='right'><a href='#TOP'>TOP</a></td></tr></table>";
			shift 																@in_fin;
			$out_html															= $out_html."<table><tr>";
					my (@headers)											= split /\t/, $in_fin[0];
					shift															@in_fin;
					foreach my $value									(@headers){ $out_html = $out_html."<td class='header' align='center'>".$value."</td>";	}
					$out_html													= $out_html."</tr>";
			my ($odd)															= 0;
			foreach my $line											(@in_fin){
				my ($class)													= "data1";
						if															( 0 == $odd % 2 ) { $class = "data2"; }
				my (@values) 												= split /\t/, $line;
				$out_html														= $out_html."<tr>";
				$values[0]													=~ s/ |\:|\_/\-/g;
				$values[0]													=~ s/\>//g;
				if 																	($sample eq '01_Expression'){
					$values[0]												=~ s/\-/\_/g;
					$out_html													= $out_html."<td class='$class' align='left'><b>$values[0]</b> (<a href='$values[0].html'>Cluster Analysis</a> | <i><a href='$values[0]_summary.html'>Hit Summary</a></i>)</td>";	
				}else{
					$out_html													= $out_html."<td class='$class' align='left'><b>$values[0]</b></td>";	
				}
				shift																@values;
				foreach my $value										(@values){ $out_html = $out_html."<td class='$class' align='center'>".$value."</td>";	}
				$out_html														= $out_html."</tr>";
				$odd++;
			}
			$out_html															= $out_html."</table><br><br>";		
		}
		my @html_log														= read_file("//data/_scripts/jrk/_resources/_html/template.html");
		foreach my $line												(@html_log){
			$line																	=~ s/%WIDTH/870/;
			$line																	=~ s/%PAGENAME/Samples/;
			$line																	=~ s/%CONTENT/$out_html/;
		}
		write_file															($var_pat."/_reports/samples.html",@html_log);
	}
}
print $out_log															bk_basic_writeLog("-----------------$var_nam.pl Pipeline Completed-----------------");
close $out_log;

sub bk_version{ 														#METHOD: B<bk_version> Tested:20131217. Time:0 secs.
	return 																		bk_basic_version()."\n".bk_biology_version()."\n".bk_wrappers_version()."\n";
}																						# B<SECTION: Methods>
sub PHAGE_version{													#METHOD: B<bk_version> Tested:20130207. Time:0 secs.
	return																		"$var_nam\tv$var_rev\tUpdated:\t$var_upd\n";
}



__END__

=head1 NAME

Using <name>

=head1 SYNOPSIS

sample <name>.pl [options]

Options:

	-help	Brief help message
	-man	Full documentation
	-ver	Version and last revision date
	-sil	Write to terminal (default silent).
	-pat	Path to working directory (default current).
	-cor	Number of cores to use for parallelization (default 1).
	-inp	Input sequence (default all sequences).
	-sli	Slide (default 7)
	-osi	Oligo window (default 36)
	-are	Turn off array creation (default on, 0)
	-cus	Turn on custom array (default off, 0)
	-lib	Turn off library evaluation (default on, 0)
	-exp	Expression threshold (default 1)
	-dup	Duplicate removal (default on, 0)
	-vir	Virus to highlight from array design (default EBOV).
 
=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-ver>

Prints the version and last revision date of the program and exits.

=item B<-sil>

Write to terminal (default to file). The program is set to report to the error
log. cmdOUT-<date>-$var_name.logf. This command redirects traffic 
back to the terminal. 

=item B<-pat>

This is the path to the working directory and defaults to the current 
directory. Specify from the partition ie. //data/<folder>/<folder>. Defualts to
the current directory.

=item B<-inp>

If you wish to restrict the number of sequences, generally when testing, you can set 
this variable to the number of sequences you wish to use. It will use the first number 
of sequence you specify. The default is to use all of the sequence.

=back

=head1 DESCRIPTION

<long description>.

=head1 REVISION HISTORY
B<3.2.0> 20140626: Continued work -jrk
B<3.2.1> 20180530: Continued work -kyg

=cut