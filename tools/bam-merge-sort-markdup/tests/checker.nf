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

params.aligned_lane_bams = ""
params.ref_genome_gz = ""

include '../bam-merge-sort-markdup.nf' params(params)


Channel
  .fromPath(params.aligned_lane_bams, checkIfExists: true)
  .set { aligned_lane_bams_ch }

Channel
  .fromPath(getMdupSecondaryFile(params.ref_genome_gz), checkIfExists: true)
  .set { ref_genome_gz_ch }

// will not run when import as module
workflow {
  main:
    bamMergeSortMarkdup(
      aligned_lane_bams_ch.collect(),
      file(params.ref_genome_gz),
      ref_genome_gz_ch.collect()
    )

  publish:
    bamMergeSortMarkdup.out to: "outdir", overwrite: true
}
