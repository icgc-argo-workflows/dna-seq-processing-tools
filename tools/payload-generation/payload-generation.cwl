class: CommandLineTool
cwlVersion: v1.1
id: payload-generation
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/dna-seq-processing-tools:payload-generation.initial'

baseCommand: [ 'payload-generation.py' ]

inputs:
  input_seq_format:
    type: string?
    inputBinding:
      prefix: -s
  bundle_type:
    type: string
    inputBinding:
      prefix: -t
  payload_schema_version:
    type: string
    inputBinding:
      prefix: -p
  metadata_lane_seq:
    type: File?
    inputBinding:
      prefix: -m
  file_to_upload:
    type: File
    inputBinding:
      prefix: -f
    secondaryFiles: [.bai?, .crai?]
  lane_seq_analysis:
    type: File[]?
    inputBinding:
      prefix: -a



outputs:
  payload:
    type: File
    outputBinding:
      glob: '$(inputs.bundle_type).*.json'


