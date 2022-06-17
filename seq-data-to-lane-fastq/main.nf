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
    Linda Xiang
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/dna-seq-processing-tools.seq-data-to-lane-fastq'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
params.metadata_json = ""
params.seq_files = ""
params.reads_max_discard_fraction = 0.05
params.tempdir = "NO_DIR"


process seqDataToLaneFastq {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir ? true : false

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    path metadata_json
    path seq

  output:  // output, make update as needed
    path "out/*{fq,fastq,fq.gz,fastq.gz}", emit: lane_fastq
    path "out/rgs_file_pair_map.csv", emit: file_pair_map_csv

  script:
    // add and initialize variables here as needed

    arg_tempdir = params.tempdir != 'NO_DIR' ? "-t ${params.tempdir}" : ""

    """
    mkdir -p out

    main.py \
      -p ${metadata_json} \
      -s ${seq} \
      -d ${params.reads_max_discard_fraction} \
      -n ${params.cpus} \
      -o out ${arg_tempdir}
    
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  seqDataToLaneFastq(
    file(params.metadata_json),
    Channel.fromPath(params.seq_files).collect()
  )
}
