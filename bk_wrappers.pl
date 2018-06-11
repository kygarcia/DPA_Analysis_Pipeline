#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

use File::Slurp;														#USE:
use Math::NumberCruncher;										#USE:

																						#B<SCRIPT INFO: "bk_wrappers", "0.0.17", "20180611">
my ($varNam, $varRev, $varUpd, $varDev)			= ("bk_wrappers", "0.0.17", "20180611", 0);
my ($patBla)																= ("//data/_pipeline_software/ncbi-blast-2.2.28+/bin"); #replace with direct path to script

unless($varDev eq 1){												#B<GROUP: Miscellaneous>
	sub bk_wrappers_version {									#METHOD: B<bk_wrappers_version> Tested:20131218. Time:0 sec.
		return 																	"$varNam\tv$varRev\tUpdated:\t$varUpd";
																						#ARGS: None
																						#RETURN: Scalar with Script Information.
																						#DESCRIPTION: Provides versioning information for the script.
	}  
}

unless($varDev eq 1){												#B<GROUP: Alignment>
	sub bk_wrappers_paraTabBlast {					  #METHOD: B<bk_wrappers_paraTabBlast> Tested:20131219. Time:537 reads(fasta)/sec.
		my (%args)														= (@_);
    
		my ($path) 														= bk_basic_returnPath($args{Query_Fasta}); 
		my ($file)														= bk_basic_returnFile($args{Query_Fasta});
		my $forkmanager  											= new Parallel::ForkManager($args{Core});
		bk_basic_paraFileSplit								( File 			  => $args{Query_Fasta},	
																					Type_Number	=> 2,
                                          Core			  => $args{Core},	      );
		my (@files) 													= <$path/_tmp/*>;
			foreach my $file (@files) { 
				$forkmanager											-> start and next;
																					`$patBla/$args{Blast_Program} -db $args{Database}$args{GI_List} -out $file.out -query $file$args{Task} -max_target_seqs=$args{Number_Hits} -outfmt "7 $args{Output_Format}"`;
				$forkmanager											-> finish; 
			}
			$forkmanager												-> wait_all_children;
																					`cat $path/_tmp/*.out > $args{Query_Fasta}.out`;
																					`rm -Rf $path/_tmp`;

										                      #ARGS: *Caution: Make sure to include a space and the tag for 
										                      #" -gilist" and " -task" all other variables just use the value.
										                      #Core 			    => 50,
										                      #Blast_Program  => "blastn",
										                      #Database 		  => $path."/".$reference.".fasta",
										                      #GI_List 		    => "",
										                      #Query_Fasta 	  => $path."/".$name.".fa",
										                      #Output_Format  => "qseqid sstart send length pident qcovhsp",
										                      #Number_Hits 	  => 1,
										                      #Task 			    => " -task blastn-short"
										                      #RETURN: None.
										                      #DESCRIPTION: This method accepts a number of BLAST parameters
										                      #an indexed database file and a query fasta file. It then splits
										                      #the file based on the number of core-1 requested, performs a 
										                      #a tabular blast ouyput as requested, consolidates the results 
										                      #and cleans up after itself.
	}
  
	sub bk_wrappers_formatBlast {						  #METHOD: B<bk_wrappers_formatBlast> Tested:20131218. Time:32.2K reads(fasta)/sec.
		my 	(%args) 													= (@_);
																					`$patBla/makeblastdb -in $args{Path_Fasta} -dbtype $args{Database_Type}`;
		
																					#ARGS: 
																					#Path_Fasta => $file
																					#Database_Type => "nucl"
																					#RETURN: None.
																					#DESCRIPTION: This method accepts a fasta file and and a database type:
																					#nucl for nucleotide and prot for protein. It then calls makeblastdb to
																					#create the appropriate blast database.
	}
	
}
	
1;