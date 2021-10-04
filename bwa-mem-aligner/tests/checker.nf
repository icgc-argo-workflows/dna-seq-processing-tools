#!/usr/bin/env nextflow

/*
  Copyright (C) 2021,  icgc-argo

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Junjun Zhang
    Linda Xiang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo/dna-seq-processing-tools.bwa-mem-aligner'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.input_bam = "NO_FILE1"
params.aligned_lane_prefix = 'grch38-aligned'
params.ref_genome_gz = "NO_FILE2"
params.sequencing_experiment_analysis = "NO_FILE3"
params.tempdir = "NO_DIR"
params.publish_dir = ""

params.expected_output = ""

include { bwaMemAligner } from '../main'
include { getBwaSecondaryFiles; getSecondaryFiles as getSec } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.1/main'


def getBase(main_file){

  parts = main_file.split("\\.").toList()
  parts.remove(0)
  baseFile = parts.join(".")

  return baseFile
}

process file_smart_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    tuple val(baseName), file(grouped_files)

  output:
    stdout()

  script:
    """
    diff <(samtools view --no-PG ${grouped_files[0]} | sort) <(samtools view --no-PG ${grouped_files[1]} | sort) \
      && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, bam files mismatch." && exit 1 )
    """
}


workflow checker {
  take:
    input_bam
    ref_genome_gz
    ref_genome_gz_secondary_files
    sequencing_experiment_analysis
    tempdir
    expected_output

  main:
    bwaMemAligner(
      input_bam,
      ref_genome_gz,
      ref_genome_gz_secondary_files,
      sequencing_experiment_analysis,
      tempdir,
      true
    )

    expected_output.concat(bwaMemAligner.out.aligned_bam)
    .map { file -> tuple(getBase(file.name), file) }
    .groupTuple()
    .set { grouped_files }

    file_smart_diff(
      grouped_files
    )
}

workflow {

  checker(
    Channel.fromPath(params.input_bam, checkIfExists: true),
    file(params.ref_genome_gz),
    Channel.fromPath(getBwaSecondaryFiles(params.ref_genome_gz), checkIfExists: true).collect(),
    file(params.sequencing_experiment_analysis),
    params.tempdir,
    Channel.fromPath(params.expected_output, checkIfExists: true)
  )
}
