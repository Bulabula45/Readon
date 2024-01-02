## Readon

Readon: a novel algorithm to identify readthrough transcripts with long-read sequencing data



## Build Instructions

1. Clone this repository locally using the command `git clone https://github.com/Bulabula45/Readon.git && cd Readon`.
2. Compile the project using  `make`.



## Usage

```sh
readon -i <reads_file> -o <output_file> [-m <mode>]
```



### Options

- `-i <reads_file>`: Specify the input file containing reads.
- `-o <output_file>`: Specify the output file to store results.
- `-m <mode>`: Specify the mode. If set to 'error', the program will operate in error mode, considering sequencing data errors. Otherwise, a larger k-mer size will be used by default.
