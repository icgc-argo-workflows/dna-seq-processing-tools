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

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.2.0.1'

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/dna-seq-processing-tools.bam-merge-sort-markdup'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.aligned_lane_bams = ""
params.ref_genome_gz = ""
params.aligned_basename = "grch38-aligned.merged"
params.markdup = true
params.output_format = "cram"
params.lossy = false
params.tempdir = "NO_DIR"
params.expected_output = ""
params.expected_metrics_output = "NO_FILE"

include { bamMergeSortMarkdup } from '../main'
include { getSecondaryFiles } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.1/main'

process output_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path output_file
    path expected_file
    path ref_genome_gz
    path ref_genome_gz_secondary_file

  output:
    stdout()

  script:
    """
    diff <(samtools view -T ${ref_genome_gz} --no-PG ${output_file} | sort) <(samtools view -T ${ref_genome_gz} --no-PG ${expected_file} | sort) \
      && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, bam files mismatch." && exit 1 )
    """
}

process metrics_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path output_file
    path expected_file

  output:
    stdout()

  script:
    """
    mkdir output expected
    tar xzf ${output_file} -C output
    tar xzf ${expected_file} -C expected
    
    cd output
    # only compare txt file
    for f in *; do 
      if [ ! -f "../expected/\$f" ]
      then
        echo "Test FAILED, found unexpected file: \$f in the output tarball" && exit 1
      fi
      echo diff \$f ../expected/\$f
      EFFECTIVE_DIFF=`diff \$f ../expected/\$f | egrep '<|>' || true`
      if [ ! -z "\$EFFECTIVE_DIFF" ]
      then
        echo -e "Test FAILED, output file \$f mismatch:\n\$EFFECTIVE_DIFF" && exit 1
      fi
    done
    echo "All files match, test PASSED" && exit 0
    """
}

workflow checker {
  take:
    aligned_lane_bams
    ref_genome_gz
    ref_genome_gz_secondary_file
    tempdir
    expected_output
    expected_metrics_output

  main:
    bamMergeSortMarkdup(
      aligned_lane_bams,
      ref_genome_gz,
      ref_genome_gz_secondary_file,
      tempdir
    )

    output_diff(
      bamMergeSortMarkdup.out.merged_seq,
      expected_output,
      ref_genome_gz,
      ref_genome_gz_secondary_file
    )

    if (params.markdup) {
      metrics_diff(
        bamMergeSortMarkdup.out.duplicates_metrics,
        expected_metrics_output
      )
    }

}


workflow {
  checker(
    Channel.fromPath(params.aligned_lane_bams, checkIfExists: true).collect(),
    file(params.ref_genome_gz),
    Channel.fromPath(getSecondaryFiles(params.ref_genome_gz, ['fai', 'gzi']), checkIfExists: true).collect(),
    params.tempdir,
    file(params.expected_output),
    file(params.expected_metrics_output)
  )
}
