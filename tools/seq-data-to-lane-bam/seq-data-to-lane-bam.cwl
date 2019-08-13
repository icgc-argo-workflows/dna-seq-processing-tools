class: CommandLineTool
cwlVersion: v1.0
id: seq-data-to-lane-bam
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/dna-seq-processing-tools:seq-data-to-lane-bam.update'

baseCommand: [ 'seq-data-to-lane-bam.py' ]

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
  picard_jar:
    type: File?
    inputBinding:
      position: 4
      prefix: -j


outputs:
  lane_bams:
    type: File[]
    outputBinding:
      glob: '*.lane.bam'
  aligned_basename:
    type: string
    outputBinding:
      glob: preprocess.json
      loadContents: true
      outputEval: |
        ${
           var data = JSON.parse(self[0].contents)["aligned_basename"];
           return data;
         }
  payload_type:
    type: string
    outputBinding:
      glob: preprocess.json
      loadContents: true
      outputEval: |
        ${
           var data = JSON.parse(self[0].contents)["payload_type"];
           return data;
         }


stdout: preprocess.json

