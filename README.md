# ecc_workflow

  * This analysis employs a modified version of Pierre Joubert's tool, [ecc_caller](https://github.com/pierrj/ecc_caller) to call eccDNAs in a small circle-seq experiment of *Oryza sativa* cv. Nipponbare.
    - Untreated Plants (labeled RC, 3 biological reps)
      + 3 technical reps of each biological rep
    - Infected plants (labeled IF, 3 biological reps)
      + 3 technical reps of each biological rep

### Preparing the Data
  * First, you will need to move all fastq.gz files into `raw_data/`.
    - `bin/units.tsv` explains names and conditions of the files needed.
    - *You are welcome to run this workflow with other samples, but you will need to edit `bin/units.tsv` to match your inpt data.*
    - *If you want to change the histone mark datasets used in the analysis, you'll need to re-write the necessary rules in `Snakefile` as this part of the analysis is currently hard-coded in an effort to ensure reproducibility of this analysis.*

### Set-Up Snakemake
  * This project is designed to be cloned to Savio, the UC Berkeley HPC, and run using the SLURM scheduler via jobscripts automatically written by [snakemake](https://snakemake.readthedocs.io/en/stable/) and submitted to savio's SLURM scheduler.
  * the user will need to have `conda` available and initialized for their profile on savio.
    - I'd recommend using this one-liner modify your `~/.bashrc` so that `conda` can be called by batch scripts submitted to savio.
      ```
      $ echo 'module load python/3.6' >> ~/.bashrc
      ```
    - personally, i've chosen to use the conda installation available to all savio users via the python3.6 module:
      ```
      $ which conda
      /global/software/sl-7.x86_64/modules/langs/python/3.6/condabin/conda
      ```
    - Once you are able to call the `conda` executable in your interactive sessions, you should initialize `conda`. This modifies your `~/.bashrc`, enabling your batch jobs to call `conda` to activate environments. To initialize:
      ```
      $ conda init
      ```
    - Now you need to create an environment and install `snakemake`. This one-liner creates an environment called `snakemake`, looks in the [bioconda](https://anaconda.org/bioconda) channel for programs and installs [snakemake (v5.25.0)](https://snakemake.readthedocs.io/en/v5.25.0/)
      ```
      $ conda create -n snakemake -c conda-forge  snakemake==5.25.0
      ```
    - This is essential so that when you submit the workflow to savio, it is able to source your `~/.bashrc` and activate a conda environment called `snakemake` where snakemake>=5.14.0 is installed. this program will be the brain organizing the analysis by installing all required software into discrete environments for each job and activating them as required by each analysis step.
    - i used the python 3.6 module managed by savio sysadmins here:`/global/software/sl-7.x86_64/modules/langs/python/3.6/bin/conda`
  * Snakemake identifies the jobs to be run by designing a diacyclic graph from the bottom-up using the desire output requested by the user in `/Snakefile`.
    - the first *rule* is a special one called `rule all:`. Here, the user specifies all the desired outputs under the directive called `input:`.
    - in my workflow, the outputs from each of the rules seen below `rule all:` the output from each intermediate rule leading to the final output have been listed for clarity, but commented out as they are not explicitly requested. Snakemake is able to automatically determine the dependencies between rules by matching file names. read more about it in the documentation [example workflow](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html).

### Running the workflow
    - To see which jobs will be queued, activate your snakemake environment perform a dry-run at the base of the repo:
      ```
      # activate snakemake conda environment
      $ conda activate snakemake

      # perform dry run
      (snakemake)$ snakemake -np
      ```
    - You should see a long list of all the jobs to be run.
    - If you'd like a visualization of the jobs, produce a .png of the DAG created by snakemake for job scheduling with the following command:
      ```
      $ snakemake --dag | dot -Tpng > my_DAG.png
      ```

### Analyzing the results

#### Epigenetic Marks
  * `H3K4ac`
    - **activation**
  * `H3K9me1`
    - **repression**
    - TSS of active genes
  * `H3K9me3`
    - **repression**
    - Satellite repeats, telomeres, pericentromeres
    - permanent signal for heterochromatin formation in gene-poor chromosomal regions with tandem repeat structures, such as satellite repeats, telomeres, and pericentromeres. It also marks retrotransposons and specific families of zinc finger genes (KRAB-ZFPs)
  * `H3K9ac`
    - active
  * `H3K27me3`
    - **repression**
    - gene-rich regions, promoters
    - temporary signal at promoter regions that controls development regulators in embryonic stem cells, including Hox and Sox genes
  * `H3K27ac`
    - **activation**
  * `Acetylation`
    - Acetylation adds a negative charge to lysine residues on the N-terminal histone tails that extend out from the nucleosome.
    - These negative charges repel negatively charged DNA, which results in a relaxed chromatin structure.
    - The open chromatin conformation allows transcription factor binding and significantly increases gene expression (Roth et al., 2001)
  * `Methylation`
    - lysine methylation is implicated in both transcriptional activation and repression depending on the methylation site
    - this flexibility may be explained by the fact that that methylation does not alter histone charge or directly impact histone-DNA interactions, unlike acetylation.
