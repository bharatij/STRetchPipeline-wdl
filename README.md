# STRetchPipeline-wdl
Workflows for detecting short tandem repeat expansions using STRetch (https://github.com/Oshlack/STRetch) on Terra Platform.
These workflows are designed for Whole Genome Sequencing starting from mapped bam/cram files.
The whole pipeline is divided into three steps as follows:

### STRetchPipeline-Steps
1. Generate STR decoy reference genome (optional)
2. Read realignment and generating STR read counts (for each sample in parallel)
3. Estimate allele size and One vs all samples outlier test (for each chromosome in parallel)

  
