class: CommandLineTool
cwlVersion: v1.0
id: bwa-mem-aligner
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/dna-seq-processing-tools:0.1.1'

baseCommand: [ 'bwa-mem-aligner.py' ]

inputs:
  input_bam:
    type: File
    inputBinding:
      prefix: -i
  ref_genome_gz:
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles: [.amb, .sa, .pac, .ann, .bwt, .fai, .alt]
  cpus:
    type: int?
    inputBinding:
      prefix: -n
  aligned_lane_prefix:
    type: string
    inputBinding:
      prefix: -o

outputs:
  aligned_lane_bam:
    type: File
    outputBinding:
      glob: $(inputs.aligned_lane_prefix).$(inputs.input_bam.basename)

