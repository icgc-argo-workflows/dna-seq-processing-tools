#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).
 *                                                                                                               
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

nextflow.preview.dsl=2

params.metadata_json = "input/1a1fbac3-00bf-4606-9fba-c300bf46068d.bam.sequencing_experiment.song-analysis.json"
params.seq_files = "input/test_rg_3.v2.bam"
params.reads_max_discard_fraction = -1
params.tool = ""

include '../seq-data-to-lane-bam' params(params)

workflow {
  main:
    seqDataToLaneBam(
      file(params.metadata_json),
      Channel.fromPath(params.seq_files).collect(),
      params.reads_max_discard_fraction,
      params.tool
    )
    // seqDataToLaneBam.out.lane_bams.view()

  publish:
    seqDataToLaneBam.out.lane_bams to: 'outdir', overwrite: true
}
