class: CommandLineTool
cwlVersion: v1.1
id: payload-generation
requirements:
- class: NetworkAccess
  networkAccess: true
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: 'quay.io/icgc-argo/payload-ceph-submission:payload-ceph-submission.initial'
- class: EnvVarRequirement
  envDef:
    - envName: "AWS_SHARED_CREDENTIALS_FILE"
      envValue: $(inputs.credentials_file.path)

baseCommand: [ 'payload-ceph-submission.py' ]

inputs:
  credentials_file:
    type: File
  payload:
    type: File
    inputBinding:
      prefix: -p
  metadata:
    type: File
    inputBinding:
      prefix: -m
  endpoint_url:
    type: string
    inputBinding:
      prefix: -s
  bucket_name:
    type: string
    inputBinding:
      prefix: -b

outputs:
  payload:
    type: File
    outputBinding:
      glob: '*.json'

