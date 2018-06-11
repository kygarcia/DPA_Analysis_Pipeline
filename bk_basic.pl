#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

use Parallel::ForkManager;														#USE:
use POSIX qw/strftime/;																#USE:
use Excel::Writer::XLSX;															#USE:

																											#B<SCRIPT INFO: "bk_basic", "0.0.13", "20180611">
my ($varNam, $varRev, $varUpd, $varDev)								= ("bk_basic", "0.0.13", "20180611", 0);

unless($varDev eq 1){																	#B<GROUP: Miscellaneous>
	sub bk_basic_version{																#METHOD: B<bk_basic_version> Tested:20131217. Time:0 secs.
		return 																						"$varNam\tv$varRev\tUpdated:\t$varUpd";
		
																											#ARGS: None
																											#RETURN: Scalar with Script Information.
																											#DESCRIPTION: Provides versioning information for the script.
	}

	sub bk_basic_uniqueArray {													#METHOD: B<bk_basic_uniqueArray> Tested:20131217. Time:1.7M records/sec.
		my %h;
		return 																						grep { !$h{$_}[0]++ } @_;
		
																											#ARGS: Array of values.
																											#RETURN: Array with unique entries.
																											#DESCRIPTION: Reduces the redundant values of an array and returns.
	}
	
	sub bk_basic_unique2DArray {												#METHOD: B<bk_basic_unique2DArray> #needs args
		my $element																				= shift;
		my %h;
    return 																						grep { !$h{$_->[$element]}++ } @_
	}                            
	
	sub bk_basic_divide {				              					#METHOD: B<bk_basic_divide> Tested:20131217. Time:0 sec.
		my	(%args) 					                  					= (@_);		
		my	($outNum) 					                					= (0) x 1; 
				if																						($args{denominator} and $args{numerator}){
					if						                      				($args{denominator} ne 0 and $args{numerator} ne 0)
																												{ $outNum = $args{numerator}/$args{denominator}; }
				}
		$outNum 						                    					= bk_basic_round
																												( significant_digits	=> $args{significant_digits},
																													number				      => $outNum 						        );
		return							                    					$outNum;

																											#ARGS:	
																											# numerator 					=> sum(@tmpCounts),
																											# denominator					=> $line[6]
																											# significant_digits	=> $args{significant_digits}
																											#RETURN: Quotient of operation.
																											#DESCRIPTION: 	Checks for null denominators and returns a 0 if found.
	}
	
	sub bk_basic_argCapture {														#METHOD: B<bk_basic_argCapture>
		my	(@args) 																			= (@_);		
		my	($outArgs) 																		= "";
		foreach my $ARG																		(@args){															
			$outArgs 																				= $outArgs." ".$ARG;
		}
		return																						$outArgs;
			
																											#ARGS: None.
																											#RETURN: String with Arguements.
																											#DESCRIPTION: 	Accepts the command line arguements and formats
																											#them for reporting.
	}

	sub bk_basic_operationLog {													#METHOD: B<bk_basic_operationLog>
		my	(%args) 																			= (@_);		
		my $oLOG																					= $args{File};
		print $oLOG 																			bk_basic_writeLog("------------------------$args{Name} Started------------------------");
		print $oLOG 																			bk_basic_writeLog("Parent Process: $args{ProcessID}");
		print $oLOG 																			bk_basic_writeLog("Arguments Used: $args{Args}");
		print $oLOG																				bk_basic_writeLog("$args{Name} Version:");
		print $oLOG																				$args{Version};
		print $oLOG																				bk_basic_writeLog("Core Components:");
		print $oLOG																				bk_version();
		print $oLOG 																			bk_basic_writeLog("3rd party software.");
			
																											#ARGS: 
																											#	Path 		=> $varPat,
																											#	Time 		=> $varTim,
																											#	Name 		=> $varNam,
																											#	ProcessID 	=> getppid,
																											# 	Args		=> $varArg,
																											# 	Version		=> VSALIGN_version(),
																											# 	ThirdPary 	=> "BLAST, NGEN",
																											#	File		=> $oLOG
																											#RETURN: None.
																											#DESCRIPTION: 	Accepts script information and initiates the log. 
																											#A list of used 3rd party software are parsed in bk_wrappers to 
																											#provide versioning information.
	}

	sub bk_basic_errorLog {															#METHOD: B<bk_basic_errorLog>
		my	(%args) 																			= (@_);		
		
		open STDOUT, 																			'>>', "$args{Path}/$args{File}" or die "Can't redirect STDOUT: $!";
		open STDERR, 																			">&STDOUT" or die "Can't redirect STDERR: $!";
			
																											#ARGS: 
																											#	Path => "//data/../../",
																											#	File => name
																											#RETURN: None.
																											#DESCRIPTION: 	Reroutes STDERROR and STDOUT to an error log file
																											#provided.
	}

	sub bk_basic_directoryStructure {										#METHOD: B<bk_basic_directoryStructure>
		my	(%args) 																			= (@_);	
		
		foreach my $folder																(@{$args{Folders}}){
			unless 																					(-d $args{Path}."/".$folder) { mkdir $args{Path}."/".$folder; }
		}
					
																											#ARGS: 
																											#	Path 	=> "//data/../../",
																											#	Folders => @Folders
																											#RETURN: None.
																											#DESCRIPTION: 	Creates the directory structure for the project.
																											#Error trapped to not try to recreated said folder.
	}

	sub bk_basic_ExcelfromArray_Phage {									#METHOD: B<bk_basic_ExcelfromArray>
		my	(%args) 																			= (@_);	
	
		my ($workbook)																		= Excel::Writer::XLSX -> new($args{Path}."/".$args{Name});
						my $format																= $workbook -> add_format();
								$format																-> set_num_format('@');
								$format																-> set_align('left');
								$format																-> set_font('Courier New');
						my	$format2															= $workbook -> add_format();
								$format2															-> set_num_format('@');
								$format2															-> set_align('left');
								$format2															-> set_font('Courier New');
								$format2															-> set_bold(1);
								if																		($args{VIR} eq "EBOV"){
								$format2															-> set_color('navy');
								}
						my $format3																= $workbook -> add_format();
								$format3															-> set_num_format('@');
								$format3															-> set_align('left'); 
								$format3															-> set_font('Courier New');
								$format3															-> set_color('gray');
						my $format4																= $workbook -> add_format(); 
								$format4															-> set_num_format('@');
								$format4															-> set_align('left');
								$format4															-> set_font('Courier New');
								$format4															-> set_color('navy');
		my (@sheet_array)																	= @{$args{Worksheets}};

		foreach my $sheet																	(@sheet_array){
			my @lines																				= read_file($args{Path}."/".$sheet);
			my $worksheet																		= $workbook -> add_worksheet(bk_basic_removeFileType(bk_basic_returnFile($sheet)));
			for(my $row=0; $row<@lines; $row++){
				my @values																		= split /\t/, $lines[$row];
				for(my $col=0; $col<@values; $col++){
					chomp $values[$col];
					if($lines[$row] =~ /EBOV/){
					$worksheet																	-> write($row, $col, $values[$col], $format2);}
					elsif($lines[$row] =~ /RESTV/){
					$worksheet																	-> write($row, $col, $values[$col], $format3);}
					elsif($args{VIR} ne "" && $lines[$row] =~ /$args{VIR}/){
					$worksheet																	-> write($row, $col, $values[$col], $format4);}
					else{
					$worksheet																	-> write($row, $col, $values[$col], $format);}
				}
			}
		}
	
																											#ARGS: 
																											#	Worksheets 		=> /@worksheets,
																											#	Path 					=> $var_Pat,
																											#	Workbook_Name	=> "projSummary.xls"
																											#RETURN: None.
																											#DESCRIPTION: 		
	}
 	
	sub bk_basic_log10 {																#METHOD: B<bk_basic_log10>
		my $in_num 																				=	shift;
		return 																						log($in_num)/log(10);
	}
	
}

unless($varDev eq 1){																	#B<GROUP: Formating Methods>
	sub bk_basic_writeLog {															#METHOD: B<bk_basic_writeLog> Tested:20131217. Time:0 sec.
		my ($string)																			= (@_);
		return 																						strftime("%Y%m%d %H:%M:%S", localtime(time())).": ".$string."\n";
		
																											#ARGS:	String to be formatted for log.
																											#RETURN: Formatted string.
																											#DESCRIPTION: Formats string for log output. 
	}

	sub bk_basic_round {																#METHOD: B<bk_basic_round> Tested:20131217. Time:0 sec.
		my																								(%args) = (@_);
		return 																						sprintf("%.".$args{significant_digits}."f", $args{number});
		
																											#ARGS:	
																											#significant_digits	=> 5
																											#number				=> $minAlleles
																											#RETURN: Formatted number.
																											#DESCRIPTION: Rounds up to the specified significant digit. 
	}

	sub bk_basic_time {																	#METHOD: B<bk_basic_time> Tested:20131217. Time:0 sec.
		my ($time)																				= (@_);
		return 																						strftime("%Y%m%d-%H_%M_%S", localtime($time));	
		
																											#ARGS:	Time as a scalar.
																											#RETURN: Formatted time string.
																											#DESCRIPTION: Formats time YYYYMMDD-HR_MIN_SEC. 
	}
	
}

unless($varDev eq 1){																	#B<GROUP: Path operations>
	sub bk_basic_returnFile {														#METHOD: B<bk_basic_returnFile> Tested:20131217. Time:0 sec.
		my ($path)																				= (@_);
		my (@elements) 																		= split "/", $path;
		return 																						$elements[@elements-1];
	
																											#ARGS:	String path to file or Directory.
																											#RETURN: String containing last element of the path.
																											#DESCRIPTION: This method returns the file or directory specified by the 
																											#path.
	}

	sub bk_basic_removeFileType {												#METHOD: B<bk_basic_removeFileType> Tested:20131217. Time:0 sec.
		my ($path) 																				= (@_);
		my (@outFile)																			= split /\./,$path;
				pop																						(@outFile);
		return 																						join ".", @outFile;
		
																											#ARGS: File name.
																											#RETURN: String removing file type.
																											#DESCRIPTION: This method removes file type and returns the name.
	}
	
	sub bk_basic_returnPath {														#METHOD: B<bk_basic_returnPath> Tested:20131217. Time:0 sec.
		my ($path) 																				= (@_);
		my (@outPath) 																		= split /\//,$path;
				pop 																					@outPath;
		return 																						join "/",@outPath;
		
																											#ARGS: String path to file or Directory.
																											#RETURN: Returns the path to the parent directory.
																											#DESCRIPTION: This method strips the file or directory from the 
																											#path and returns the path to the parent directory.
	}
}

unless($varDev eq 1){																	#B<GROUP: File operations>
	sub bk_basic_paraGunzip {														#METHOD: B<bk_basic_paraGunzip> Tested:20131217. Time:66.7K reads(fastq)/sec.
		my (%args)																				= (@_);
		my $forkManager																		= new Parallel::ForkManager($args{Core});
		foreach my $file 																	(@{$args{Files}}){ 	
			$forkManager																		-> start and next;
			my $tmp																					= bk_basic_removeFileType($file);
			`gunzip -c $file > $tmp.uz`;
			$forkManager																		-> finish; 
		}
		$forkManager																			-> wait_all_children;
		
																											#ARGS: 
																											#   Core  => $core, 
																											#   Files => \@sequence
																											#RETURN: None.
																											#DESCRIPTION: This method accepts an array of file paths and unzips them,  
																											#in parallel by the nubmer of core specified. The files are unzipped and 
																											#a .uz extension is appended. The original file remains after the creation  
																											#of the unzipped version.
	}
	
	sub bk_basic_lineCount {														#METHOD: B<bk_basic_lineCount> Tested:20131217. Time:>1M reads(fastq)/sec.
		my (%args)																				= (@_);
		my ($lines) 																			= (`wc -l $args{File}`); 
		my (@outLines)																		= split " ", $lines;
				if 																						(@outLines ne 0){ return $outLines[0]/$args{Type_Number}; }
				else																					{ return 0;	}
		
																											#ARGS: 
																											#   File 			=> $file_path, 
																											#   Type_Number 	=> 2
																											#RETURN: Number of lines or a 0 as a scalar.
																											#DESCRIPTION: This method accepts a file path and a type. Types are for
																											#splitting a file in discreet units. Fastq uses 4 lines per sequence and 
																											#Fasta uses 2. Enter the numeric value for the appropriate type to make 
																											#sure the split does not split a unit.
	}
	
	sub bk_basic_paraFileSplit {												#METHOD: B<bk_basic_paraFileSplit> Tested:20131217. Time:200K reads(fastq)/sec.
		my (%args)																				= (@_);
		my ($path, $outFile)															= (bk_basic_returnPath($args{File}),bk_basic_returnFile($args{File}));
		unless																						(-d $path."/_tmp"){ mkdir $path."/_tmp"; }
		my ($lineNumber)																	= bk_basic_lineCount
																												( File 			=> $args{File},	
																													Type_Number	=> $args{Type_Number} );
				$lineNumber																		= ($lineNumber/($args{Core}-1))+1;
				$lineNumber																		= bk_basic_round
																												( significant_digits => 0,
																													number			 => $lineNumber );
				$lineNumber																		= $lineNumber * $args{Type_Number};		
																											`split -d -l $lineNumber $args{File} $path/_tmp/$outFile.`;
		
																											#ARGS: 
																											#   File 					=> $file_path, 
																											#   Type_Number 	=> 2
																											#   Core					=> 10
																											#RETURN: None.
																											#DESCRIPTION: This method accepts a file and a file type. Types are for
																											#splitting a file in discreet units. Fastq uses 4 lines per sequence and 
																											#Fasta uses 2. Enter the numeric value for the appropriate type to make 
																											#sure the split does not split a unit. The file is broken up into n files
																											#where n is the number of core - 1 to account for the parent process. The
																											#files are created in an _tmp directory in the same location as the file.
										#																	Files are number 00-99.
	}
}

1;