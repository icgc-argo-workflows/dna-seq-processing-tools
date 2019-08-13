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
      prefix: -sf
  payload_type:
    type: string
    inputBinding:
      prefix: -pt
  payload_schema_zip:
    type: File
    inputBinding:
      prefix: -ps
  metadata_lane_seq:
    type: File?
    inputBinding:
      prefix: -mls
  file_to_upload:
    type: File
    inputBinding:
      prefix: -fu
    secondaryFiles: [.bai?, .crai?]
  lane_seq_analysis:
    type: File[]?
    inputBinding:
      prefix: -lp



outputs:
  payload:
    type: File
    outputBinding:
      glob: '$(inputs.payload_type).*.json'


