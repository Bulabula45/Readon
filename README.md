## Readon

Readon: a novel algorithm to identify readthrough transcripts with long-read sequencing data




## Build Instructions and Examples

1. Clone this repository locally using the command `git clone https://github.com/Bulabula45/Readon.git && cd Readon`. Or `wget https://github.com/Bulabula45/Readon/archive/refs/heads/main.zip && unzip main.zip && cd Readon-main`

3. Compile the project using  `make`.

4. Uncompress reference files using the command `tar -xzvf hg38REF.tar.gz`. Upon extraction, the directory will include transcriptome files grouped by chromosome and strand, along with a geneName file generated to record each gene and its corresponding multiple transcripts. These files serve as the reference for readon's alignment operations(There is also an instruction for building references files of other organisms).

5. You can use the following command to test whether readon can run without any problems and produce the expected output:

   `./readon -i test/example.fa -o test/out`

6. If you want to use readon globally, you can add the directory containing the executable file readon to the environment variable, and then source it. Similarly, you can use `readon -i test/example.fa -o test/test_out` to run the test.

   ```sh
   export PATH=$PATH:/path/to/readon_directory
   source ~/.bashrc  # if you are using macOS or if your system is configured to use ~/.bash_profile
   ```



   

## Usage and options

```sh
Usage: ./readon -i <reads_file> -o <output_file> [-m <mode>] [-k <kmer>] [-s <step>] [-w <window_size>] [-h]
```

- `-i <reads_file>`: Specify the input file containing reads(required).
- `-o <output_file>`: Specify the output file to store results(required).
- `-m <mode>`: Specify the mode. If set to 'error', the program will operate in error mode, considering sequencing data errors. Otherwise, a larger k-mer size will be used by default.
- `-k <kmer>` : Specify the k-mer size (default is 25), maximum is 32
- `-s <step> `: Specify the interval value (default is 10. If error mode is error, the default value is 6)
- `-w <window_size>` : Specify the window_s value (default is 15)
- `-h` : Display this help message





## Downstream Analysis 

We provide two downstream analysis tools. First, you need to uncompress the `analysis_tools.tar.gz` file by `tar -xzvf analysis_tools.tar.gz `, which will yield two files: `NMD_or_Protein_prediction.ipynb` and `visualization.ipynb`.

`./analysis_tools/NMD_or_Protein_prediction.ipynb` is used for predicting whether a read-through transcript is likely to undergo nonsense-mediated decay or encodes a protein. The notebook contains detailed instructions guiding you on what reference files to download, how to process these files, and the necessary bash commands to execute.

 `./analysis_tools/visualization.ipynb` is used for visualizing splicing patterns.

Dependencies include: subprocess, matplotlib, and biopython, which can be easily installed via pip.





## Construct References for Other Organisms

As an example, we use the mouse. We need to prepare its transcriptome reference sequences, including gene names, chromosome names, strands, transcription start and end positions. The specific file processing depends on the format of your data.

For example, we can download the fasta file using martview provided by Emsembl(https://www.ensembl.org/biomart/martview/), and we choose Dataset - Mouse genes (GRCm39), and attributes including "Transcript stable ID", "cDNA sequences", "Gene name", "Transcript start (bp)", "Transcript end (bp)", "Chromosome/scaffold name", "Strand"(ordered).

Then, these sequences are grouped by chromosome and strand, and reordered to generate new files according to the order of each gene on the chromosome. Meanwhile, a geneName file is generated to record each gene and its corresponding multiple transcripts. These commands are included in `./analysis_tools/preprocess_for_other_organisms.ipynb`. If your data format differs, the part that needs modification is section 1.2, where sequence descriptions are processed.







