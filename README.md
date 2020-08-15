# dsbra
Double Strand Break Analysis Program <br>
Analyze .fastq files for insertions, deletions, and microhomologous end joining at specified break site within reference sequence. Example output:
<p align="center">
  <img src="https://github.com/lingchen42/dsbra/blob/master/example_plots/RIF1_NFR_b84_m5_repair_pattern_summary_table_aligned_mutation_events.png"/>
</p>
<p align="center">
  <img src="https://github.com/lingchen42/dsbra/blob/master/example_plots/RIF1_NFR_b84_m5_repair_pattern_summary_table_sequences_with_mutation_events.png"/>
</p>

## Usage
**1. Load necessary programs** (Here I use LMOD, which is applied to running this on vanderbilt ACCRE cluster.)
```
ml Intel git Anaconda2
```
**2. Clone the git repo and cd into the directory** (This step only needs to be done the first time.)
```
git clone https://github.com/lingchen42/dsbra.git
cd dsbra
```
**3. Install dependencies within a conda virtual enviroment** (This step only needs to be done the first time. If the dsbra virtual enviroment already exists, please remove the old one or give the new one a different name, using the `-n` flag.)
```
conda env create -n dsbra -f dsbra_conda.yml
```
**4. Activate the virtual enviroment**
```
source activate dsbra
```
**5. Run the `sam_dsbra_interface_ling.py` to infer the repair patterns**
```
usage: sam_dsbra_interface_ling.py [-h] -o OUT_DIR
                                   [--align_mode {global,local}] --pcr_primer1
                                   PCR_PRIMER1 --pcr_primer2 PCR_PRIMER2
                                   --pcr_tail_seq PCR_TAIL_SEQ
                                   [--min_len MIN_LEN] [--not_valid_reads]
                                   [-r REF_FA] [-q FASTQ] [-b BREAK_INDEX]
                                   [-m MARGIN] [-qm Q_REPAIR_MARGIN]
                                   [--last_margin LAST_MARGIN]
                                   [--count_cutoff COUNT_CUTOFF COUNT_CUTOFF]
                                   [--full_table]

Double strand break repair analyzer (DSBRA)

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory
  --align_mode {global,local}
                        choose the method for alignment
  --pcr_primer1 PCR_PRIMER1
                        primer seq, at the begining of the ref seq; will
                        require the read to starts with this primer seq
  --pcr_primer2 PCR_PRIMER2
                        primer seq, at the end of the ref seq; will require
                        last 4 nucleotide of the read is part of this primer
                        seq
  --pcr_tail_seq PCR_TAIL_SEQ
                        added sequence for PCR;default
                        ATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  --min_len MIN_LEN     miminal valid length (reads - pcr_tail_seq); default
                        44bp
  --not_valid_reads     output not valid reads
  -r REF_FA, --ref_fa REF_FA
                        path of the reference sequence fasta file
  -q FASTQ, --fastq FASTQ
                        path of the reads to analyze in fastq format
  -b BREAK_INDEX, --break_index BREAK_INDEX
                        break index. 0-indexed, the position ofnucleotide
                        right before the break site.
  -m MARGIN, --margin MARGIN
                        the margin before and after the break indexthat should
                        be considered for analysis;default 5bp.
  -qm Q_REPAIR_MARGIN, --q_repair_margin Q_REPAIR_MARGIN
                        the margin before and after the last mutationthat
                        should be reportedin the repaired_read_sequence;
                        default 5bp.
  --last_margin LAST_MARGIN
                        the margin outside of the last mutation eventthat
                        should be included in the repair sequencepattern for
                        those compound mutation events; default 5bp
  --count_cutoff COUNT_CUTOFF COUNT_CUTOFF
                        take probability cutoff, number of colonies,apply
                        count cutoff to the summary tablebased on the poisson
                        distribution;default 0.001, 300
  --full_table          Will generate a large table withrepair pattern for
                        each analyzed

```
_Example use:_
```
ref_seq=NFR
break_idx=85
pcr_primer1=CGTAGGAAGTAGATTGTGTTAG
pcr_primer2=-1
pcr_tail_seq=ATCGGAAGAGCAACGTCTGAACTCCAGTCAC
margin=10
n_colony=300
poisson_cutoff=0.001
python sam_dsbra_interface_ling.py -o rif1_nfr_wt/ -q example_data/RIF1_NFR.fastq -b $break_idx -m $margin --ref example_data/nfr_ref_seq.fa --pcr_primer1 $pcr_primer1 --pcr_primer2 $pcr_primer2 --pcr_tail_seq $pcr_tail_seq --not_valid_reads --full_table
```

The results will be in the `rif1_nfr_wt` directory as specified by `-o`,  including:
- `RIF1_NFR_b84_m5_repair_pattern_summary_table.txt` This is the **_main_** result table, summarizing the dsbra repair pattern observed in the reads. 
- `RIF1_NFR_b84_m5_repair_pattern_full_table.txt` This has the repair pattern for each read. If we want to confirm how the read is converted to corresponding repair pattern, we can check this table.
- `RIF1_NFR_b84_m5_run_info.txt` This is a text file documenting the setting of the script.
- `RIF1_NFR.sam*` These are the intermediate files generated by samtools when align the reads in the FASTQ file to the reference.
- `RIF1_NFR_failed_alignment.sam` This is the SAM file for failed alignments.
- `RIF1_NFR_b84_m5_ref.fa` This is the intermediate file generated by the script for storing the reference sequences.
- `RIF1_NFR_b84_m5_mut_reads.sam*`, `RIF1_NFR_b84_m5_non_mut_reads.sam` These are the intermediate filse generated by the script to store the reads with mutation within range of interest.
<br>

**6. To plot the results, run the `sam_dsbra_grapher_ling.py`.**
```
usage: sam_dsbra_grapher_ling.py [-h] [-R R] [-o OUT_DIR] -i INPUT_FILES
                                 [INPUT_FILES ...] [--run_info RUN_INFO]
                                 [--align_stats] [--all] [--mut_type]
                                 [--del_len] [--del_seq] [--ins_len]
                                 [--ins_seq] [--repair_seq]
                                 [--aligned_mutations]
                                 [--ref_bottom REF_BOTTOM] [--ref_top REF_TOP]
                                 [--break_index BREAK_INDEX] [--fts FTS]
                                 [--count_cutoff COUNT_CUTOFF COUNT_CUTOFF]

Visualize Mutation Events

optional arguments:
  -h, --help            show this help message and exit
  -R R                  R plotting script path
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory, default is plots_2020-08-14/
  -i INPUT_FILES [INPUT_FILES ...], --input_files INPUT_FILES [INPUT_FILES ...]
                        the output txt files from sam_dsbra_interface.py
  --run_info RUN_INFO   the run_info csv files from sam_dsbra_interface.py
  --align_stats         plot run_info.csv, plot alignment stats by Type
  --all                 plot all types of plots
  --mut_type            plot Mutation Event Frequency by Type
  --del_len             plot Frequency of Deletions by Length
  --del_seq             plot Sequences with Deletion Event
  --ins_len             plot Frequency of Insertions by Length
  --ins_seq             plot Sequences with Insertion Event
  --repair_seq          plot top most common repair sequences
  --aligned_mutations   plot frequency and aligned mutations
  --ref_bottom REF_BOTTOM
                        the start reference base to show aligned mutation
                        events
  --ref_top REF_TOP     the start reference base to show aligned mutation
                        events
  --break_index BREAK_INDEX
                        the index of break site
  --fts FTS             axis font size; default 6
  --count_cutoff COUNT_CUTOFF COUNT_CUTOFF
                        take probability cutoff, number of colonies,apply
                        count cutoff to the summary tablebased on the poisson
                        distribution;default 0.001, 300
```

_Example use:_<br>
If we want to plot all types of plots: <br>
```
python sam_dsbra_grapher_ling.py --all -i example_outputs/RIF1_NFR_output.txt --break_index ${break_idx} --count_cutoff $poisson_cutoff $n_colony -o  example_outputs/ --run_info example_outputs/RIF1_NFR_first100k_global_aln_b${break_idx}_m${margin}_run_info.csv --fts 5
```
The plots will be in `example_outputs` directory. <br>
More options can be found in `./sam_dsbra_grapher_ling.py -h`, such as setting a count cutoff.

**7. Finally, exit the virtual environment**
```
source deactivate
```

## Original Reference
Soong, C. P., Breuer, G. A., Hannon, R. A., Kim, S. D., Salem, A. F., Wang, G., . . . Bindra, R. S. (2015). Development of a novel method to create double-strand break repair fingerprints using next-generation sequencing. DNA Repair (Amst), 26, 44-53. doi:10.1016/j.dnarep.2014.12.002
