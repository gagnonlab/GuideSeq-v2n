# GuideSeq-v2n: Guide-Seq Analysis Pipeline for Nanopore Reads

A modified version of the Guide-Seq-2 (https://github.com/tsailabSJ/guideseq) pipeline that accomodates nanopore reads (or other single-end reads) lacking UMIs.

## Table of Contents
- [Pipeline](#pipeline)
  - [Running Pipeline Steps Individually](#steps)
- [Dependencies](#dependencies)
- [Getting Set Up](#setup)
	- [Installation](#Installation)
- [Running the Full Analysis Pipeline](#full_pipeline)
	- [Quickstart](#quickstart)
	- [Writing A Manifest File](#write_manifest)
	- [A Full Manifest File Example](manifest_example)
	- [Pipeline Outputs](#pipeline_output)



## Pipeline<a name="pipeline"></a>

The package implements a pipeline consisting of a series of modules to take a raw series of raw Nanopore fastq files from GuideSeq experiments and process them to derive potential CRISPR off-target sites. 

The individual pipeline steps are:

1. **Sample filtering**: Pre-demultiplexed nanopore fastq files are filtered by a series of criteria. Reads are checked for the full dsODN sequence, as a partial sequence could be a sign of nonspecific binding during PCR. If the dsODN is on the 3' end of the read, the read is replaced with its reverse-complement to ensure the dsODN is always on the 5' end of the read. Reads too short to align are also discarded. 

2. **Read Alignment**: The filtered reads are aligned to a reference genome using minimap2 with default map-ont parameters. After alignment, a UMI-free deduplication process is run to attempt to remove PCR duplicates. Because nanopore reads can be any length, and the Nanopore RAPID transposase inserting randomly, the odds that two reads will be identical without PCR duplication are low. Therefore, reads that align to the same location and have the same aligned length are assumed to be PCR duplicates.

3. **Candidate Site Identification**: The mapping positions of the deduplicated reads containing the dsODN sequence are tabulated on a genome-wide basis. Mapping positions are consolidated using a 10-bp sliding window. Windows with multiple reads are flagged as sites of potential DSBs. 25 bp of reference sequence (`window_size`) is retrieved on either side of the most frequently occuring start-mapping position in each flagged window. The retrieved sequence is aligned to the guide sequence using a Smith-Waterman local-alignment algorithm. Mismatches are counted. Off-target cleavage sites with more than 7 mismatches (`max_score`, mismatch counted as 1, indel counted as 3) to the intended target sequence are filtered out.

4. **Visualization**: Identified off-targets, sorted by GUIDE-Seq deduplicated read count, are annotated in a final output table. The GUIDE-Seq read count is expected to scale approximately linearly with cleavage rates (Tsai et al., *Nat Biotechnol.* 2015). Alignment of detected off-target sites is visualized via a color-coded sequence grid, as seen below:

![guideseq_flowchart](vis_example.PNG)

5. **Tablization**: Results of individual samples are combined into a single tab-seperated-value (TSV) table for ease of comparison.

### Running Pipeline Steps Individuallye<a name="steps"></a>

In addition to end-to-end pipeline analysis functionality, the guideseq package also allows for every step of the analysis to be run individually. There are 5 steps: `filter`, `align`, `identify`, `visualize`, and 'tablize', specified using `--step filter+align+identify+visualize+tablize`. You can run individual steps just by giving the single step name, or a combination of steps, using `+` to connect them. What you can't do is to skip a step, for example skipping the alignment step, e.g. `filter+identify`, will cause an error. This can be useful for troubleshooting failed steps.

## Dependencies<a name="dependencies"></a>
* Python 3
* Reference genome fasta file ([Example](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/))
* [`minimap2`](<https://github.com/lh3/minimap2>) alignment tool
* [`bedtools`](<https://bedtools.readthedocs.io/en/latest/>) genome arithmetic utility
* [`samtools`](<http://www.htslib.org/>) genome arithmetic utility
* [`homer`](<http://homer.ucsd.edu/homer/>) off-target annotation
* [`seqkit`](<https://bioinf.shenwei.me/seqkit/>) tool used to filter and process reads
* R `ggplot` and `bsgenome` for off-target manhattan plot, optional


## Getting Set Up<a name="setup"></a>

### Installation<a name="Installation"></a>

The simplest way to install the guideseq-v2n pipeline is via conda.

```
# Create a conda environment with dependencies:
conda create -n guideseq-v2n minimap2 bedtools samtools biopython homer seqkit pandas pyfaidx pysam pyyaml regex tabulate Levenshtein dill r-base r-essentials r-ggplot2 bioconda::bioconductor-bsgenome bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 svgwrite

# Homer requires additional installation steps
# Replace with your own conda environment path
perl ~/.conda/envs/guideseq-v2n/share/homer/configureHomer.pl -install hg38

# Download guideseq-v2n code:
git clone https://github.com/gagnonlab/guideseq-v2n.git

```



### Quickstart <a name="Quickstart"></a>

```

git clone https://github.com/tsailabSJ/guideseq

cd guideseq/test

python ../guideseq/guideseq.py main -m test_manifest.yaml

```

If you are using an HPC with LSF job system, you can parallelize the guide-seq analysis by `guideseq.py parallel -m test_manifest.yaml`

## Running the Full Analysis Pipeline<a name="full_pipeline"></a>

### Quickstart<a name="quickstart"></a>

To run the full guideseq analysis pipeline, you must first create a manifest YAML file that describes all pipeline inputs. Once you have done so, you can simply run

```
guideseq main -m /path/to/manifest.yaml
```

to run the entire pipeline. Below are specific instructions detailing how to write the manifest file.

### Writing A Manifest File<a name="write_manifest"></a>
When running the end-to-end analysis functionality of the guideseq package, a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest formatted as a YAML file. YAML files allow easy-to-read specification of key-value pairs. This allows us to easily specify our parameters. The following fields are required in the manifest:

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `blacklist`: The absolute or relative path to the blacklist file (optional).
- `genome`: hg19 or hg38 for homer off-target annotation.
- `analysis_folder`: The absolute or relative path to the folder in which all pipeline outputs will be saved.
- `bwa`: The absolute path to the `bwa` executable
- `bedtools`: The absolute path to the `bedtools` executable
- `PAM`: PAM sequence (optional), default is NGG.
- `demultiplex`: The absolute paths to the undemultiplexed paired end sequencing files. The required parameters are:
	- `forward`: The absolute path to the FASTQ file containing the forward reads.
	- `reverse`: The absolute path to the FASTQ file containing the reverse reads.
	- `index1`: The absolute path to the FASTQ file containing the forward index reads.
	- `index2`: The absolute path to the FASTQ file containing the reverse index reads.

An example `undemultiplexed` field:

```
demultiplex:
    min_reads: 10000
    mismatch: 1
    forward: test_fastq_R1_001.fastq
    reverse: test_fastq_R2_001.fastq
    index1: test_fastq_I1_001.fastq
    index2: test_fastq_I2_001.fastq
    out_dir: demultiplexed # relative to analysis_folder
    input_dir: ./ # relative to analysis_folder
```

- `samples`: A nested field containing the details of each sample. At least two samples must be specified: a "control" sample (to be used to filter out background off-target sites) and at least one treatment sample. The required parameters are:
	- `target`: The sample targetsites
	- `barcode1`: The forward barcode
	- `barcode2`: The reverse barcode
	- `controlbarcode1`: The forward controlbarcode1
	- `controlbarcode2`: The reverse controlbarcode2
	- `description`: A description of the sample

An example `samples` field:

```
samples:
    test:
        target: GGGGCCACTAGGGACAGGATNGG
        barcode1: CCTGAGGA
        barcode2: GTAAGGAG
        controlbarcode1: GGGGGGGG
        controlbarcode2: CGGCGACC
        description:  test_description
```

### A Full Manifest File Example<a name="manifest_example"></a>

Below is an example of a full manifest file. Feel free to copy it and replace the parameters with your own experiment data. Remember that you can input more than just one treatment sample (e.g. the "AAVS1_site_14" data below).

```
reference_genome: chr19.fa
blacklist: hg38-blacklist.v2.bed.gz
genome: hg38
# if not specified, default is NGG
PAM: NGG

# output folder directory (absolute or relative path), if not specified, default is current dir
analysis_folder: ./


# PATH to analysis tools
# If not specified, default is just the tool name, i.e., accessible in the $PATH
bwa: bwa
samtools: samtools
bedtools: bedtools
umi_tools: umi_tools

# optional, to generate the off-target Manhattan plot
Rscript: Rscript


window_size: 25
max_score: 7
mapq_threshold: 50


# DEMULTIPLEX 
# ------------
demultiplex:
    min_reads: 10000
    mismatch: 1
    forward: test_fastq_R1_001.fastq
    reverse: test_fastq_R2_001.fastq
    index1: test_fastq_I1_001.fastq
    index2: test_fastq_I2_001.fastq
    out_dir: demultiplexed # relative to analysis_folder
    input_dir: ./ # relative to analysis_folder
    subsample_reads: -1

samples:
    test:
        target: GGGGCCACTAGGGACAGGATNGG
        barcode1: CCTGAGGA
        barcode2: GTAAGGAG
        controlbarcode1: GGGGGGGG
        controlbarcode2: CGGCGACC
        description:  test_description




```

### Pipeline Output<a name="pipeline_output"></a>

When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:


#### Output Folders
- `output_folder/demultiplexed`: Contains the four undemultiplexed reads files (forward, reverse, index1, index2) for each sample.
- `output_folder/aligned`: Contains an alignment `.sam` file for each sample.
- `output_folder/identified`: Contains a tab-delimited `.txt` file for each sample with an identified off-target in each row.
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.

The final detected off-target sites are placed in the `output_folder/identified` folder, with one `.rmblck.txt` file for each sample specified in the manifest. The fields that are populated in each row of these off-target files are specified below:

####Output Off-Targets `.txt` Fields:

- `BED Chromosome`: Window chromosome
- `BED Min.Position`: Window 0-based start position
- `BED Max.Position`: Window 0-based end position
- `BED Name`: Name of window 
- `Filename`: The name of the current `.SAM` file used in analysis.
- `WindowIndex`: Index number of window
- `Chromosome`: Chromosome corresponding to position with maximum reads in window (matches `BED Chromosome`)
- `Position`: Position with maximum number of reads in window
- `Sequence`: The window sequence, starting 25 bp upstream and ending 25 bp downstream of `Chromosome:Position`
- `+.mi`: Number of forward reads with distinct molecular indices
- `-.mi`: Number of reverse reads with distinct molecular indices
- `bi.sum.mi`: Sum of the `+.mi` and `-.mi` fields (GUIDE-seq Read Count)
- `bi.geometric_mean.mi`: Geometric mean of the `+.mi` and `-.mi` fields
- `+.total`: Total number of forward mapping reads 
- `-.total`: Total number of reverse mapping reads 
- `total.sum`: Sum of `+.total` and `-.total` fields
- `total.geometric_mean`: Geometric mean of the `+.total` and `-.total` fields
- `primer1.mi`: Number of reads amplified by forward primer with distinct molecular indices
- `primer2.mi`: Number of reads amplified by reverse primer with distinct molecular indices
- `primer.geometric_mean`: Geometric mean of the `primer1.mi` and `primer2.mi` fields
- `position.stdev`: Standard deviation of positions within genomic window
- `Off-Target Sequence`: Off-target sequence derived from genome reference
- `Mismatches`: Number of mismatches between the intended target sequence and the off-target sequence
- `Length`: Length of the target sequence
- `Cells`: Cell type
- `Target site`: Targetsite name
- `Target Sequence`: Intended target site sequence (including PAM)
- The following information are in column names prefix with `Site_SubstitutionsOnly` or `Site_GapsAllowed`.
- `BED off-target Chromosome`: Off-target chromosome
- `BED off-target start`: Off-target 0-based start position
- `BED off-target end`: Off-target 0-based end position
- `BED off-target name`: Off-target name
- `BED Score`: Field to conform to standard BED format
- `Strand`: Indicates the strand of detected off-target site. `+` for forward strand and `-` for reverse strand
The key fields for interpreting this output and identifying off-target sites are: `BED off-target Chromosome`, `BED off-target start`, `BED off-target end`, `BED off-target name`, `BED off-target strand`, `Off-Target Sequence`, `bi.sum.mi`

#### Output Visualizations

The outputted visualizations are in the `.svg` vector format, which is an open image standard that can be viewed in any modern web browser (e.g. Google Chrome, Apple Safari, Mozilla Firefox), and can be viewed and edited in any vector editing application (e.g. Adobe Illustrator). Because the output visualizations are vector images, they can be scaled up or down infinitely without a loss in quality, and can also be edited as shapes with ease. This makes the images produced by the guideseq package ideal for posters, presentations, and papers.
