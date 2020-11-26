##peak overlap
./intervene/intervene venn -i intervene/example_data/ENCODE_hESC/*.bed
./intervene/intervene upset -i intervene/example_data/ENCODE_hESC/*.bed
./intervene/intervene pairwise -i intervene/example_data/dbSUPER_mm9/*.bed

## random peak 
bedtools shuffle -seed 1 -i ../../macs_out_2_peaks.narrowPeak -g ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.size  > 2_random_1.bed
bedtools shuffle -seed 2 -i ../../macs_out_2_peaks.narrowPeak -g ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.size  > 2_random_2.bed 
bedtools shuffle -seed 3 -i ../../macs_out_2_peaks.narrowPeak -g ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.size  > 2_random_3.bed

##get seq.fa
narrowPeak -fo 5_peak.fa
bedtools getfasta -fi ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.fa -bed 2_random_1.bed -fo 2_random_1.bed.fa
bedtools getfasta -fi ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.fa -bed 2_random_2.bed -fo 2_random_2.bed.fa
bedtools getfasta -fi ~/project/project/00.DATABASE/hg38/bowtie2_index/hg38.fa -bed 2_random_3.bed -fo 2_random_3.bed.fa

## motif analysis

Rscript sequence_hits_analysis_ChIP_C.R 5_peak.fa 5_random_1.bed.fa 5_random_2.bed.fa 5_random_3.bed.fa

## plot motif
Rscript plot.motif.R
