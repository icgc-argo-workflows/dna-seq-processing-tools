class: CommandLineTool
cwlVersion: v1.0
id: s3-upload
requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/s3-upload:s3-upload.init'

baseCommand: [ 's3-upload.py' ]

inputs:
  endpoint_url:
    type: string
    inputBinding:
      prefix: -s
  bucket_name:
    type: string
    inputBinding:
      prefix: -b
  bundle_type:
    type: string
    inputBinding:
      prefix: -t
  payload_json:
    type: File
    inputBinding:
      prefix: -p
  s3_credential_file:
    type: File
    inputBinding:
      prefix: -c
  upload_files:
    type: File[]
    inputBinding:
      prefix: -f

outputs: []
