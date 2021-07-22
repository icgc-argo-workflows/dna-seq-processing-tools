#!/usr/bin/env nextflow

/*
  Copyright (c) 2019-2021, Ontario Institute for Cancer Research (OICR).

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
*/

nextflow.enable.dsl = 2

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

bwaSecondaryExts = ['fai', 'sa', 'bwt', 'ann', 'amb', 'pac', 'alt']

params.file_name = null
params.file_size = null

include {
    getSecondaryFiles;
    getBwaSecondaryFiles
} from '../main.nf'

include {
    generateDummyFile as gFile1;
    generateDummyFile as gFile2;
} from './generate-dummy-file.nf'

include {
    filesExist as fExist1;
    filesExist as fExist2;
} from './files-exist.nf'

Channel.from(params.file_name).set{ file_name_ch }
Channel.from(bwaSecondaryExts).set{ bwa_ext_ch }


workflow {
    // generate the main file
    gFile1(
        file_name_ch.flatten(),
        params.file_size
    )

    // generate the BWA secondary files
    gFile2(
        file_name_ch.combine(bwa_ext_ch),
        params.file_size
    )

    // test 'getSecondaryFiles' for expected 'fai' file exists
    fExist1(
        getSecondaryFiles(params.file_name, ['fai']),
        'exist',
        gFile2.out.file.collect(),
        true  // no need to wait
    )

    // test 'getBwaSecondaryFiles' for all expected bwa secondary files exist
    fExist2(
        getBwaSecondaryFiles(params.file_name).collect(),
        'exist',
        gFile2.out.file.collect(),
        true  // no need to wait
    )

}
