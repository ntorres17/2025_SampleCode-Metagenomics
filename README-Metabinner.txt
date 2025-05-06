Creating coverage profile directly from sequencing reads: 
cd path_to_MetaBinner
cd scripts

Files required: 

contig_file=/Absolute/Path/to/File/final_contigs_f1k.fa 

output_dir=/Absolute/Path/to/File/output_test

coverage_profiles=/Absolute/Path/to/File/coverage_profile_f1k.tsv  
 
kmer_profile=/Absolute/Path/to/File/kmer_4_f1000.csv

bash gen_coverage_file.sh -a contig_file \
-o output_dir_of_coveragefile \
path_to_sequencing_reads/*fastq

bash gen_coverage_file.sh -a /Absolute/Path/to/File/Guinea_1852.fa \ 
      -o /Absolute/Path/to/File/ \
      /Absolute/Path/to/File/*.fastq


bash run_metabinner.sh -a /Absolute/Path/to/File/guinea6922.fa \
                       -o /Absolute/Path/to/File/output \
                       -d /Absolute/Path/to/File/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))


bash run_metabinner.sh -a /Absolute/Path/to/File/Garbon_12988.fa \
                       -o /Absolute/Path/to/File/output \
                       -d /Absolute/Path/to/File/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))
