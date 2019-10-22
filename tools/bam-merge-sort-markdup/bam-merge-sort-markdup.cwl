class: CommandLineTool
cwlVersion: v1.0
id: bam-merge-sort-markdup
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/bam-merge-sort-markdup:bam-merge-sort-markdup.0.1.3'

baseCommand: [ 'bam-merge-sort-markdup.py' ]

inputs:
  aligned_lane_bams:
    type: File[]
    inputBinding:
      prefix: -i
  ref_genome:
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles: [.fai]
  cpus:
    type: int?
    inputBinding:
      prefix: -n
  aligned_basename:
    type: string
    inputBinding:
      prefix: -b
  markdup:
    type: boolean?
    default: true
    inputBinding:
      prefix: -d
  output_format:
    type: string[]?
    inputBinding:
      prefix: -o
  lossy:
    type: boolean?
    default: false
    inputBinding:
      prefix: -l

outputs:
  aligned_seq:
    type:
      - "null"
      - type: array
        items: File
    outputBinding:
      glob:
        - '$(inputs.aligned_basename).bam'
        - '$(inputs.aligned_basename).cram'
    secondaryFiles: ['.bai', '.crai']
  aligned_duplicate_metrics:
    type: ["null", File]
    outputBinding:
      glob: $(inputs.aligned_basename).duplicates-metrics.txt
  bundle_type:
    type: string
    outputBinding:
      glob: stdout.json
      loadContents: true
      outputEval: |
        ${
           var data = JSON.parse(self[0].contents)["bundle_type"];
           return data;
         }

stdout: stdout.json



