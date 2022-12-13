#NOTE: The first few steps must be completed within a virtual machine with SBanalyzer from ShorelineBiome. 
#We are using SBanalyzer version 3.1-3.

#Demultiplex samples using the SBanalyzer GUI. Use the "Demux-NoTrim.txt" pipeline.
#The output folder "1_phase_find_samples" will contain fasta files for each sample.
#Our ASV-calling involves pooling the samples together, which makes it more sensitive to rare ASVs, but is computationally intensive.
#To avoid crashes and speed the process up, split the samples up into groups by moving the fasta files into new folders.
#	Target between 30 and 40 samples per folder. Note that computation time increases exponentially, so <35 samples per folder is ideal.

#Run sb-dada2 within the ShorelineBiome virtual machine. 
sb-dada2 -n Run1 -i SequencingRun1 -s StrainID -c 1e-120

#Assign taxonomy with Athena database within the ShorelineBiome virtual machine.
/opt/sbanalyzer/bin/sbsearch \
--mode Search --seqs Run1_ASVs_only_nochiR_md5.fasta \
--db /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/index.bin \
--tax /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/athena_v2_2.tax \
-op /home/shoreline/Documents/Run1

#Repeat for other groups. 
sb-dada2 -n Run2 -i SequencingRun2 -s StrainID -c 1e-120

/opt/sbanalyzer/bin/sbsearch \
--mode Search --seqs Run2_ASVs_only_nochiR_md5.fasta \
--db /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/index.bin \
--tax /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/athena_v2_2.tax \
-op /home/shoreline/Documents/Run2

sb-dada2 -n Run3 -i SequencingRun3 -s StrainID -c 1e-120

/opt/sbanalyzer/bin/sbsearch \
--mode Search --seqs Run3_ASVs_only_nochiR_md5.fasta \
--db /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/index.bin \
--tax /home/shoreline/Documents/SBanalyzer_Master/.user_lib/lib/athena_v2_2/athena_v2_2.tax \
-op /home/shoreline/Documents/Run3

#For each run, create a directory containing the *nochiR.csv, *nochiR_md5.fasta, and taxonomy.tax files.
#Place each of these directories into a shared directory. This new directory can have any name.
mkdir /home/shoreline/Documents/merge_runs

#Use post_dada2_merge.R to merge the output files together files together
Rscript post_dada2_merge.R -n merged -i /home/shoreline/Documents/merge_runs