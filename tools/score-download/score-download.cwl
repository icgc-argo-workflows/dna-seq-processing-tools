cwlVersion: v1.1
class: CommandLineTool
id: score-download
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: NetworkAccess
  networkAccess: true
- class: InitialWorkDirRequirement
  listing: $(inputs.seq_files)
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/score-download:score-download.0.1.3'

baseCommand: [ 'score-download.py' ]

inputs:
  seq_files:
    type: File[]?
    default: []
    inputBinding:
      prefix: -s
  file_tsv:
    type: File?
    inputBinding:
      prefix: -f
  repository:
    type: string?
    inputBinding:
      prefix: -r
  token_file:
    type: File?
    inputBinding:
      prefix: -t

outputs:
  seq_files:
    type: File[]
    outputBinding:
      glob: "*"
