# Repetitive_Elements_Pipeline
This pipeline is implemented in the Python programming language. The purpose of the pipeline is to insert repetitive elements into genomes, break them apart into short reads, and then reassemble them to test the ability of the assemblers to deal with repetitive elements. This pipeline was created by Hillary Dapsauski, Kimia Boreiri, and Joshua Melnick.


## Required Dependencies 
To utilize and run the pipeline, install the following dependencies and ensure that these dependencies are within the directory where you will be running the pipeline: 

- **Python** - [Download](https://www.python.org/downloads/)
- **NCBI Datasets** - [Download](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
- **ART** - [Download](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) - use the MountRainier version 
- **SPAdes** - [GitHub](https://github.com/ablab/spades)
- **Unicycler** - [GitHub](https://github.com/rrwick/Unicycler)
- **Conda** - [Download](https://anaconda.org/anaconda/conda)
- **QUAST** - [GitHub](https://github.com/ablab/quast)
- **Matplotlib** - [Download](https://matplotlib.org/stable/install/index.html)

The pipeline currently uses 4 specific genomes, but these can be changed to any other genomes by modifying the `accession` variable. 


## Running the pipeline 

To run the pipline use the following lines of code within the terminal. Make sure to run the pipeline within the directory called Repetitive_Elements_Pipeline or else this will not function properly. 

```bash
git clone https://github.com/h-dapsauski/Repetitive_Elements_Pipeline.git
cd Repetitive_Elements_Pipeline
python3 pipeline.py
```

Be mindful that this pipeline will skip certain steps depending on if it detects that certain folders already exist, indicating a prior run. This is by design to allow for replicability and to avoid rerunning long processes, but it may cause the program to behave in undesired ways. The SPAdes and Unicycler steps will skip entirely if it detects the existence of folders called `Spades_Output` or `Unicycler_Output`, respectively, so if you want to rerun those parts either delete them (not just their contents!) or move them outside of the `Repetitive_Elements_Pipeline` folder. The `run_art()` function can still run if the `artgens` folder exists, but it won't redo any files that would have the same names as any already in it.

To reuse motifs without them being lost to randomization in subsequent runs, the program will skip the generation step of `generate_and_insert_repeats()` if it detects a `Motifs` folder with the files `motif1.txt`, `motif2.txt`, and `ip.txt`. 

### Sample data (IMPORTANT!)
Sample data has been provided using a truncated version of one of the genomes and a pre-generated set of motifs. If you run the sample data, **remove the names of the other 3 genomes** from the `accession` variable on line 12, leaving just `GCF_004379335.1`. (Note: the QUAST output produced by the sample run will be empty for the iterations with coverage values of 10, likely due to the short sequence)

If you are **NOT** running the sample data and you don't plan to reuse anything from previous runs, the only things that should be in your `Repetitive_Elements_Pipleline` folder prior to beginning the run should be `pipeline.py` and this README. **DELETE** the `Genomes` and `Motifs` folders that come with the pipeline as the sample data will conflict with your run.

## Pipeline steps
### Downloading the genomes
The pipeline begins by downloading the genomes from RefSeq using the `datasets download` command if they have not already been downloaded. 
The 4 genomes we use by default are:
* Spiroplasma gladiatoris (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004379335.1/)
* Cruoricaptor ignavus (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_014961145.1/)
* Pseudomonas sivasensis (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_021391435.1/)
* Iamia majanohamensis (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_028532485.1/)

However, the pipeline can be modified to use any genome. Simply update the `accession` variable in the `download_genomes` function with your desired accession codes.

### Creation and insertion of repetitive elements
Next, repetitive elements are inserted into the pipeline to create artificial genomes. From a "pool" of nucleotides of "A", "C", "T", and "G", the repetitive element is randomly generated at the specified k-value or base pair length (in our case we had a k-value of 100 and 500). These repetitive elements are stored in a text file called 'motif.txt'. Then, a list of insertion points is randomly generated, and then the repetitive element is inserted either 2 times, 3 times, 4 times, or 5 times. These newly modified genomes are stored in a directory called 'Modified_Genomes'. Each file is named `accession#`_`motif#`_`repeat#`. 

### Short read simulation
The artificial genomes are broken apart by ART as it simulates 150 bp paired end short reads from an Illumina sequencer. It does this at a depth of 10 and 100, but these values can be modified by changing the `depth` variable in the for loop of the `run_art` function. An example of a modified for loop that increments depth by 10 from 10 to 100 is available in the comments above the current for loop. 

The call to ART uses the following flags: 
* `-p` indicates paired ends
* `-na` omits the alignment files so ART just outputs the short reads
* `-i` is the input genome, which in this case is each modified genome
* `-l` is the read length, currently set to 151
* `-f` is the fold coverage, which means depth
* `-o` is the output prefix used for naming
* `-m` is the mean fragment length, set to 200
* `-s` is the standard deviation, set to 10

### Assembly
The short reads are then reassembled by SPAdes and Unicycler to compare their ability to handle repetitive elements. Both SPAdes and Unicycler take as input the paired-end short reads simulated by ART. For each genome and each sequencing depth (e.g., 10x and 100X) and both tools are run on the same sets of reads to allow a direct comparision.

## SPAdes
In this pipeline, SPAdes reconstructs the original genome from the paired-end short reads simulated by ART. The assembled genome is saved in a folder specific to each genome and coverage depth.

The call SPAdes uses the following flags:
* `-1` is the forward read file (R1)
* `-2` is the reverse read file (R2)
* `-t` sets the number of threads used
* `-o` is the output directory

## Unicycler
In this pipeline, unicycler is run in short-read mode using paired-end reads simulated by ART and assembled genome is saved in a folder specific to each genome and coverage depth.

The call to Unicycler uses the following flags:
* `-1` is the forward read file (R1)
* `-2` is the reverse read file (R2)
* `-t` sets the number of threads used
* `-o` is the output directory

## Install Conda and QUAST (unused)
The `install_conda_and_quast()` function is used to automatically install Miniconda and create a Conda environment for running QUAST. When this function is run, it first checks whether Miniconda has already been installed. After Miniconda is installed, the function ensures that the Conda environment exists. If it does not, important channels are added.To install QUAST inside the conda environment, the script uses the command `-y` to automatically approves the installation without asking for confirmation.

## QUAST
The `run_quast` function runs QUAST to compare genome assemblies generated by SPAdes and Unicycler. It  uses the modified genome (with inserted repetitive elements) as reference to evaluate the accuracy of each assembler.

The call to QUAST uses the following flags:
* `{spades_dir} (unicycler_dir}` input folders for SPAdes and Unicycler assemblies
* `-r` is the modified reference genome
* `-l` labels for the assemblies
* `-o` output directory for quast results 

### Interpretation
An additional, separate script has been provided to create a table from the QUAST data. However, this script implements hard-coding and is unlikely to work with different data. However, it has been provided for testing and reference purposes.
