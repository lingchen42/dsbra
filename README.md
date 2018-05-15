# dsbra
Double Strand Break Analysis Program
Analyze .fastq files for insertions, deletions, and microhomologous end joining at specified break site within reference sequence.

## Run On ACCRE
1. Load necessary packages with LMOD.
```
ml Intel git Anaconda2 Bowtie2/2.3.2
```
2. Clone the git repo to ACCRE and cd into the directory
```
git clone https://github.com/lingchen42/dsbra.git
cd dsbra
```
3. Install dependencies within a conda virtual enviroment. (This step only needs to be done the first time.)
```
conda env create -n dsbra -f dsbra_conda.yml
```
4. Activate the virtual enviroment.
```
source activate dsbra
```
5. Run the scripts. 
```
python sam_dsbra_interface.py -r <reference.fa> -q <sample.fastq|.fq> -b <break_index>,<margin>,<last_margin> -o <output_filename.txt> [-v] [-c]
```
example use:
```
python sam_dsbra_interface.py -r example_data/nfr_ref_seq.fa -q example_data/RIF1_NFR.fastq -b 100,10,5 -o temp_out.txt
```
The results will be in the `temp_out.txt` file. <br>
This python script calls the bowtie2 and samtools, therefore will also generate some other files. Those files, together with running
information will be stored in `temp_outputs/` directory. <br>

6. Finally, exit the virtual environment
```
source deactivate
```

## Original Documentation
Use:
1. Must install samtools AND bowtie2 OR bwa aligner and edit approprite paths within sam_dsbra_interface.py if not in your PATH variable.
2. sam_dsbra_interface.py requires pysam, numpy, scipy installation prior to use
3. sam_dsbra_grapher.py additionally requires matplotlib for graph generation

sam_dsbra_interface.py -r <reference.fa> -q <sample.fastq|.fq> -b <break_index>,<margin>,<last_margin> -o <output_filename.txt> [-v] [-c]

Example use:
python sam_dsbra_interface.py -r ./reference.fa -q ./my_sample.fq -b 100,10,5 -o my_output.txt

Original Reference can be found here:

Soong, C. P., Breuer, G. A., Hannon, R. A., Kim, S. D., Salem, A. F., Wang, G., . . . Bindra, R. S. (2015). Development of a novel method to create double-strand break repair fingerprints using next-generation sequencing. DNA Repair (Amst), 26, 44-53. doi:10.1016/j.dnarep.2014.12.002
