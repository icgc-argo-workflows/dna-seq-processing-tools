class: CommandLineTool
cwlVersion: v1.0
id: seq-validation
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: InitialWorkDirRequirement
  listing: $(inputs.seq_files_dir.listing)
- class: DockerRequirement
  dockerPull: 'quay.io/pancancer/dna-seq-processing:0.1.0'

baseCommand: [ 'seq-validation.py' ]

inputs:
  seq_format:
    type: string
    inputBinding:
      position: 1
      prefix: -i
  seq_rg_json:
    type: File
    inputBinding:
      position: 2
      prefix: -p
  seq_files_dir:
    type: Directory
    inputBinding:
      position: 3
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

