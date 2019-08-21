class: CommandLineTool
cwlVersion: v1.0
id: seq-validation
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/dna-seq-processing-tools:seq-validation.update'

baseCommand: [ 'seq-validation.py' ]

inputs:
  seq_rg_json:
    type: File
    inputBinding:
      position: 1
      prefix: -p
  seq_files:
    type: File[]
    inputBinding:
      position: 2
      prefix: -d


outputs:
  seq_status:
    type: string
    outputBinding:
      glob: seq_validation.json
      loadContents: true
      outputEval: |
        ${
           var data = JSON.parse(self[0].contents)["valid"];
           return data;
         }

stdout: seq_validation.json

