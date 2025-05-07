Creating coverage profile directly from sequencing reads: 
cd path_to_MetaBinner
cd scripts

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory (to save the coverage files)
	   -b STR          directory for the bam files (optional)
        -t INT          number of threads (default=1)
        -m INT          amount of RAM available (default=4)
        -l INT          minimum contig length to bin (default=1000bp).
        --single-end    non-paired reads mode (provide *.fastq files)
        --interleaved   input read files contain interleaved paired-end reads
        -f              Forward read suffix for paired reads (default="_1.fastq")
	   -r              Reverse read suffix for paired reads (default="_2.fastq")

Files required: 

contig_file=/Absolute/Path/to/File/final_contigs_f1k.fa 

output_dir=/Absolute/Path/to/File/output_test

coverage_profiles=/Absolute/Path/to/File/coverage_profile_f1k.tsv  
 
kmer_profile=/Absolute/Path/to/File/kmer_4_f1000.csv

____________________________
Project folder - it works!: 

4 = Kmer 
1000 = contig length threshold 

gen_kmer.py /Absolute/Path/to/File/GA12988/Garbon_12988.fa 1000 4
resulting file: Garbon_12988_kmer_4_f1000.csv

gen_kmer.py /Absolute/Path/to/File/GU6922/guinea6922.fa 1000 4
resulting file: guinea6922_kmer_4_f1000.csv

gen_kmer.py /Absolute/Path/to/File/Lab/nbm_MAGs.fa 1000 4

gen_kmer.py /Absolute/Path/to/File/Lab/pbm_MAGs.fa 1000 4

gen_kmer.py /Absolute/Path/to/File/Lab/sm_MAGs.fa 1000 4

____________________________
bash gen_coverage_file.sh -a contig_file \
-o output_dir_of_coveragefile \
path_to_sequencing_reads/*fastq

bash gen_coverage_file.sh -a /Absolute/Path/to/File/GU6922/guinea6922.fa \
-o /Absolute/Path/to/File/GU6922/Coverage_profile_GU \
--interleaved /Absolute/Path/to/File/GU6922/fq_files/guinea_reads.fastq 

bash gen_coverage_file.sh -a /Absolute/Path/to/File/GA12988/Garbon_12988.fa \
-o /Absolute/Path/to/File/GA12988/Coverage_profile_GA \
--interleaved /Absolute/Path/to/File/GA12988/fq_files/gabon_reads.fastq 

bash gen_coverage_file.sh -a /Absolute/Path/to/File/Lab/nbm_MAGs.fa \
-o /Absolute/Path/to/File/Lab/Coverage_profile_nbm \
--interleaved /Absolute/Path/to/File/Lab/fq_files/nbm.fastq 

bash gen_coverage_file.sh -a /Absolute/Path/to/File/Lab/pbm_MAGs.fa \
-o /Absolute/Path/to/File/Lab/Coverage_profile_pbm \
--interleaved /Absolute/Path/to/File/Lab/fq_files/pbm.fastq 

bash gen_coverage_file.sh -a /Absolute/Path/to/File/Lab/sm_MAGs.fa \
-o /Absolute/Path/to/File/Lab/Coverage_profile_sm \
--interleaved /Absolute/Path/to/File/Lab/fq_files/sm.fastq 

____________________________
bash run_metabinner.sh -a /Absolute/Path/to/File/GU6922/guinea6922.fa \
                       -o /Absolute/Path/to/File/GU6922/Output_MetaBinner_All \
                       -d /Absolute/Path/to/File/GU6922/Coverage_profile_All/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/GU6922/guinea6922_kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))

bash run_metabinner.sh -a /Absolute/Path/to/File/GA12988/Garbon_12988.fa \
                       -o /Absolute/Path/to/File/GA12988/Output_MetaBinner_GA \
                       -d /Absolute/Path/to/File/GA12988/Coverage_profile_GA/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/GA12988/Garbon_12988_kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))

bash run_metabinner.sh -a /Absolute/Path/to/File/Lab/nbm_MAGs.fa \
                       -o /Absolute/Path/to/File/Lab/Output_MetaBinner_nbm \
                       -d /Absolute/Path/to/File/Lab/Coverage_profile_nbm/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/Lab/nmb_MAGs_kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))

bash run_metabinner.sh -a /Absolute/Path/to/File/Lab/pbm_MAGs.fa \
                       -o /Absolute/Path/to/File/Lab/Output_MetaBinner_pbm \
                       -d /Absolute/Path/to/File/Lab/Coverage_profile_pbm/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/Lab/pmb_MAGs_kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))

bash run_metabinner.sh -a /Absolute/Path/to/File/Lab/sm_MAGs.fa \
                       -o /Absolute/Path/to/File/Output_MetaBinner_sm \
                       -d /Absolute/Path/to/File/Lab/Coverage_profile_nbm/coverage_profile_f1k.tsv \
                       -k /Absolute/Path/to/File/Lab/sm_MAGs_kmer_4_f1000.csv \
                       -p $(dirname $(which run_metabinner.sh))