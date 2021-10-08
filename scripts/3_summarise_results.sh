cd loh_pipeline_results_icgc
cat */sequenza/purity_ploidy.txt | grep -v ^Ploidy | less
for i in *; do cat $i/polysolver/results/winners.hla.txt | perl -ne 'BEGIN {$pt = shift @ARGV}; chomp; print "$pt\t$_\t"' $i; echo ""; done | less # polysolver results all in one line
for i in *; do cat $i/lohhla/results/*HLAlossPrediction_CI.xls | grep -v BAF | sort | uniq | perl -ne 'BEGIN {$pt = shift @ARGV}; print "$pt\t$_"' $i; done | less # lohhla results on separate lines
