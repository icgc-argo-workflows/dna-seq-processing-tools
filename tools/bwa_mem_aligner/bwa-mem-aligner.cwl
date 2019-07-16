class: CommandLineTool
cwlVersion: v1.0
id: bwa-mem-aligner
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/pancancer/dna-seq-processing:latest'

baseCommand: [ 'bwa-mem-aligner.py' ]

inputs:
  input_bam:
    type: File
    inputBinding:
      prefix: -i
  ref_genome:
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles: [.amb, .sa, .pac, .ann, .bwt, .fai, .alt]
  cpus:
    type: int?
    inputBinding:
      prefix: -n
  aligned_seq_output:
    type: string
    inputBinding:
      prefix: -o

outputs:
  aligned_seq:
    type: File
    outputBinding:
      glob: $(inputs.aligned_seq_output)
  output_meta:
    type: File
    outputBinding:
      glob: output.json
