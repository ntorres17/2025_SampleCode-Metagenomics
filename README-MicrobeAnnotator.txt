MicrobeAnnotator 

conda install -c bioconda diamond
conda install pandas
conda install matplotlib
conda install seaborn
conda install psutil
conda install attrs
pip install wget
pip install biopython

# Building database

nohup microbeannotator_db_builder -d MicrobeAnnotator_DB -m diamond -t 7 --light &

nohup microbeannotator_db_builder -d MicrobeAnnotator_DB -m diamond -t 7 & # Full version 

Vanilla template code: 

microbeannotator -i [fasta_1.fa fasta_2.fa] \
                 -d [microbeannotator_db_dir] \
                 -o [output folder] \
                 -m [blast,diamond,sword] \
                 -p [# processes] \
                 -t [# threads]

Multi-input Full 

nohup microbeannotator -i /Absolute/Path/to/File/aa_seqs/*.fa \
-d /Absolute/Path/to/File/MicrobeAnnotator_DB \
-o /Absolute/Path/to/File/Refined_Input_Full \
-m diamond \
-p 5 \
-t 2 \
--refine
&

Multi-input Light

nohup microbeannotator -i /Absolute/Path/to/File/aa_seqs/*.fa \
-d /Absolute/Path/to/File/MicrobeAnnotator_DB \
-o /Absolute/Path/to/File/Refined_Input_Light \
-m diamond \
-p 5 \
-t 2 \
--light \
--refine
&