# Nextflow Package `seq-data-to-lane-fastq`

A simple wrapper written in `nextflow` for the sequencing processing tool to convert all input sequencing data into unaligned and lane level fastq files.  
The tool support both aligned bam or unaligned fastq formats with paired or single end reads.

## Package development

The initial version of this package was created by the WorkFlow Package Manager CLI tool, please refer to
the [documentation](https://wfpm.readthedocs.io) for details on the development procedure including
versioning, updating, CI testing and releasing.


## Inputs
### Required
- `metadata_json`: JSON file contains donor/sample/specimen/experiment/read_groups/files metadata
- `seq_files`: Sequencing reads in aligned BAM or unaligned FASTQ formats. Supported input format: {BAM, *.fq.gz, *.fastq.gz, *.fq.bz2, *.fastq.bz2}

### Optional
- `reads_max_discard_fraction`: Max fraction of reads allowed to be discarded when reverting aligned BAM to unaligned
- `tempdir`: Specify directory for temporary files
- `cpus`: Set cpu number for running the tool
- `mem`: Set memory(G) for running the tool
- `publish_dir`: Specify directory for getting output files

## Outputs
- `lane_fastq`: All fastq files 
- `file_pair_map_csv`: CSV file contains the 3 columns per lane: `read_group_id`, `file_r1`, `file_r2` 

## Usage

### Run the package directly

With inputs prepared, you should be able to run the package directly using the following command.
Please replace the params file with a real one (with all required parameters and input files). Example
params file(s) can be found in the `tests` folder.

```
nextflow run icgc-argo-workflows/dna-seq-processing-tools/seq-data-to-lane-fastq/main.nf -r seq-data-to-lane-fastq.v0.1.0 -params-file <your-params-json-file>
```

### Import the package as a dependency

To import this package into another package as a dependency, please follow these steps at the
importing package side:

1. add this package's URI `github.com/icgc-argo-workflows/dna-seq-processing-tools/seq-data-to-lane-fastq@0.1.0` in the `dependencies` list of the `pkg.json` file
2. run `wfpm install` to install the dependency
3. add the `include` statement in the main Nextflow script to import the dependent package from this path: `./wfpr_modules/github.com/icgc-argo-workflows/dna-seq-processing-tools/seq-data-to-lane-fastq@0.1.0/main.nf`
