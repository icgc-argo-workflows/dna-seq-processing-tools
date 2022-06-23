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

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
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

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.metadata_json = ""
params.seq_files = ""
params.reads_max_discard_fraction = 0.05
params.tempdir = "NO_DIR"

include { seqDataToLaneFastq } from '../main'

workflow {

  main:
    seqDataToLaneFastq(
      file(params.metadata_json),
      Channel.fromPath(params.seq_files).collect()
    )

  emit:
    lane_fastq = seqDataToLaneFastq.out.lane_fastq
    rg_fastq_map = seqDataToLaneFastq.out.file_pair_map_csv
}
