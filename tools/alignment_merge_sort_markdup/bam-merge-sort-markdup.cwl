class: CommandLineTool
cwlVersion: v1.0
id: bwa-mem-aligner
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/lindaxiang/dna-seq-processing:latest'

baseCommand: [ 'bam-merge-sort-markdup.py' ]

inputs:
  aligned_lane_bams:
    type: File[]
    inputBinding:
      prefix: -i
  ref_genome:
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles: [.fai]
  cpus:
    type: int?
    inputBinding:
      prefix: -n
  aligned_basename:
    type: string
    inputBinding:
      prefix: -o
  markdup:
    type: boolean
    inputBinding:
      prefix: -d
  cram:
    type: boolean
    inputBinding:
      prefix: -c

outputs:
  aligned_bam:
    type: File
    outputBinding:
      glob: $(inputs.aligned_basename).bam
    secondaryFiles: [.bai]
  aligned_bam_duplicate_metrics:
    type: File
    outputBinding:
      glob: $(inputs.aligned_basename).bam.duplicates-metrics.txt
  aligned_cram:
    type: File
    outputBinding:
      glob: $(inputs.aligned_basename).cram
    secondaryFiles: [.crai]





