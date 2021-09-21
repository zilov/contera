# contera
Contera - a tool to detect and delete non-target organisms and adapters contamination in your genome assembly. Contera Uses BLAST against NT and custom adapters databases to find contaminations.

## Dependencies
Contera is written and tested on python >= 3.5. To run it you should install BLAST NT database on your system (~ 110GB)

## Install

To install contera with Conda use following command:

`conda create -c contera_env -c zilov contera`

And than activate the environment: `conda activate contera_env`

## Download or set BLAST NT database

If you do not have database in your system, you can download it with following command (**Please make sure that you have 120 Gb of free space!**):

`contera -m download_db -o /path/to/outdir`

If you already have database in your system, please add CONTERA_DB in your PATH with the following command (Make sure tha database was build with makeblastdb!):

`export CONTERA_BLAST=/path/to/blast_db_folder/nt`

## Usage

Contera has three modes to run:

1) Contamination - to check for contamination with BLAST NT database and delete it from your genome assembly, to run it use:

    `contera -m cont -a /path/to/assembly.fasta -db /path/to/blast_db -o /path/to/outdir`

2) Adapters - to check for contamination with adapters database, to run it use:

    `contera -m adapter -a /path/to/assembly.fasta -o /path/to/outdir`

3) All-in-one (default) - to check for both adapters and other organisms contaminations, to run it use:

    `contera -a /path/to/assembly.fasta -db /path/to/blast_db -o /path/to/outdir`

## Options
```
  -h, --help            show this help message and exit
  -m {aio, cont, adapter, download_db}, --mode {aio, cont, adapter, download_db}
                        mode to use [default = aio]
  -a ASSEMBLY, --assembly ASSEMBLY
                        path to genome assembly in FASTA format
  -p PREFIX, --prefix PREFIX
                        prefix for output files
  -l LENGTH, --length LENGTH
                        length threshold to delete short contigs
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -db BLAST_DB, --blast_db BLAST_DB
                        path to blast nucleotide (nt) database
  -t THREADS, --threads THREADS
                        number of threads [default == 8]
  -d, --debug           do not run, print commands only
  ```
