# dsbra
Double Strand Break Analysis Program

Analyze .fastq files for insertions, deletions, and microhomologous end joining at specified break site within reference sequence.

Use:
1. Must install samtools, bowtie2 or bwa aligner and edit approprite paths within sam_dsbra_interface.py if not in your PATH variable.
2. sam_dsbra_interface.py requires pysam, numpy, scipy installation prior to use
3. sam_dsbra_grapher.py additionally requires matplotlib for graph generation

sam_dsbra_interface.py -r <reference.fa> -q <sample.fastq|.fq> -b <break_index>,<margin>,<last_margin> -o <output_filename.txt> [-v] [-c]

Example use:
python sam_dsbra_interface.py -r ./reference.fa -q ./my_sample.fq -b 100,10,5 -o my_output.txt

Original Reference can be found here:

Soong, C. P., Breuer, G. A., Hannon, R. A., Kim, S. D., Salem, A. F., Wang, G., . . . Bindra, R. S. (2015). Development of a novel method to create double-strand break repair fingerprints using next-generation sequencing. DNA Repair (Amst), 26, 44-53. doi:10.1016/j.dnarep.2014.12.002
