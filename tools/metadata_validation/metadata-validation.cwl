class: CommandLineTool
cwlVersion: v1.0
id: bwa-mem-aligner
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/pancancer/dna-seq-processing:latest'

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
  input_format:
    type: string
    outputBinding:
      glob: input_validation.json
      loadContents: true
      outputEval: |
        ${
           var data = JSON.parse(self[0].contents)["input_format"];
           return data;
         }

stdout: input_validation.json
