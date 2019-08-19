class: CommandLineTool
cwlVersion: v1.0
id: metadata-validation
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/dna-seq-processing-tools:0.1.1'

baseCommand: [ 'metadata-validation.py' ]

inputs:
  meta_format:
    type: string?
    inputBinding:
      prefix: -m
  exp_json:
    type: string?
    inputBinding:
      prefix: -j
  exp_tsv:
    type: File?
    inputBinding:
      prefix: -e
  rg_tsv:
    type: File?
    inputBinding:
      prefix: -r
  file_tsv:
    type: File?
    inputBinding:
      prefix: -f
  seq_exp_json_name:
    type: string
    inputBinding:
      prefix: -o
  seq_rg_json_name:
    type: string
    inputBinding:
      prefix: -p

outputs:
  seq_exp_json:
    type: File
    outputBinding:
      glob: $(inputs.seq_exp_json_name)
  seq_rg_json:
    type: File
    outputBinding:
      glob: $(inputs.seq_rg_json_name)
