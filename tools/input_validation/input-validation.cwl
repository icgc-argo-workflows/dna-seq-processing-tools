class: CommandLineTool
cwlVersion: v1.0
id: bwa-mem-aligner
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/pancancer/dna-seq-processing:latest'

baseCommand: [ 'input-validation.py' ]

inputs:
  seqexp_tsv:
    type: File
    inputBinding:
      prefix: -e
  rg_tsv:
    type: File
    inputBinding:
      prefix: -r
  seq_tsv:
    type: File
    inputBinding:
      prefix: -s
  seq_files
    type: Directory
    inputBinding:
      prefix: -f
  seqexp_json_name:
    type: string
    inputBinding:
      prefix: -o
  seq_rg_json_name:
    type: string
    inputBinding:
      prefix: -p

outputs:
  seqexp_json:
    type: File
    outputBinding:
      glob: $(inputs.seqexp_json_name)
  seq_rg_json:
    type: File
    outputBinding:
      glob: $(inputs.seq_rg_json_name)
  output_meta:
    type: File
    outputBinding:
      glob: output.json
