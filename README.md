# dsbra
Double Strand Break Analysis Program
Analyze .fastq files for insertions, deletions, and microhomologous end joining at specified break site within reference sequence.

## Run On ACCRE
1. Load necessary packages with LMOD.
```
ml Intel git Anaconda2 Bowtie2/2.3.2
```
2. Clone the git repo to ACCRE and cd into the directory (This step only needs to be done the first time.)
```
git clone https://github.com/lingchen42/dsbra.git
cd dsbra
```
3. Install dependencies within a conda virtual enviroment. (This step only needs to be done the first time. If the dsbra virtual enviroment already exists, please remove the old one or give the new one a different name, using the `-n` flag.)
```
conda env create -n dsbra -f dsbra_conda.yml
```
4. Activate the virtual enviroment.
```
source activate dsbra
```
5. Run the `sam_dsbra_interface.py` to infer the mutation events. 
```
python sam_dsbra_interface_ling.py -r <reference.fa> -q <sample.fastq|.fq> -b <break_index>,<margin>,<last_margin> -o <output_filename.txt> [-v] [-c]
```
Example use:
```
python sam_dsbra_interface_ling.py -r example_data/nfr_ref_seq.fa -q example_data/RIF1_NFR.fastq -b 84,10,5 -o temp_out.txt
```
The results will be in the `temp_out.txt` file. <br>

6. To plot the results, run the `sam_dsbra_grapher_ling.py`.
```
./sam_dsbra_grapher_ling.py -h
usage: sam_dsbra_grapher_ling.py [-h] [-R R] [-o OUT_DIR] -i INPUT_FILES
                                 [INPUT_FILES ...] [--all] [--mut_type]
                                 [--del_len] [--del_seq] [--ins_len]
                                 [--ins_seq]

Visualize Mutation Events

optional arguments:
  -h, --help            show this help message and exit
  -R R                  R plotting script path
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory, default is plots_2018-05-30/
  -i INPUT_FILES [INPUT_FILES ...], --input_files INPUT_FILES
[INPUT_FILES ...]
                        the output txt files from sam_dsbra_interface.py
  --all                 plot all types of plots
  --mut_type            plot Mutation Event Frequency by Type
  --del_len             plot Frequency of Deletions by Length
  --del_seq             plot Sequences with Deletion Event
  --ins_len             plot Frequency of Insertions by Length
  --ins_seq             plot Sequences with Insertion Event

```
Example use:
If we want to plot Mutation Event Frequency by Type, use: <br>
```
./sam_dsbra_grapher_ling.py -i example_outputs/RIF1_NFR_output.txt -o example_outputs/ --mut_type
```
The Mutation Event Frequency by Type plot will be `example_outputs/RIF1_NFR_output_mutation_event_frequency_by_type.png`. <br>
<br>
If we want to plot all types of plots: <br>
```
./sam_dsbra_grapher_ling.py -i example_outputs/RIF1_NFR_output.txt -o example_outputs/ --all
```
The plots will be in `example_outputs` directory. <br>

7. Finally, exit the virtual environment
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
