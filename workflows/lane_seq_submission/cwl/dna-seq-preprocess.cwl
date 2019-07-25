#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
id: dna-seq-preprocess

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:
  exp_tsv: File
  rg_tsv: File
  file_tsv: File
  seq_exp_json_name: string
  seq_rg_json_name: string
  seq_files_dir: Directory
  picard_jar: File


outputs:
  lane_bams:
    type: File[]
    outputSource: preprocess/lane_bams
  aligned_basename:
    type: string
    outputSource: preprocess/aligned_basename


steps:
  metadata_validation:
    run: ../../../tools/metadata_validation/metadata-validation.cwl
    in:
      exp_tsv: exp_tsv
      rg_tsv: rg_tsv
      file_tsv: file_tsv
      seq_exp_json_name: seq_exp_json_name
      seq_rg_json_name: seq_rg_json_name
    out:
      [ seq_exp_json, seq_rg_json, input_format]

  sequence_validation:
    run: ../../../tools/sequence_validation/seq-validation.cwl
    in:
      seq_rg_json: metadata_validation/seq_rg_json
      seq_files_dir: seq_files_dir
      seq_format: metadata_validation/input_format
    out:
      [  ]

  preprocess:
    run: ../../../tools/seq_data_to_lane_bam/preprocess.cwl
    in:
      picard_jar: picard_jar
      seq_rg_json: metadata_validation/seq_rg_json
      seq_format: metadata_validation/input_format
      seq_files_dir: seq_files_dir
    out:
      [ lane_bams, aligned_basename ]
