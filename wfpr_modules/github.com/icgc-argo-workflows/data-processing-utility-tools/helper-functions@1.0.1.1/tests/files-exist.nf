#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019-2021, Ontario Institute for Cancer Research (OICR).
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
 * Contributors:
 *   Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '1.0.1.1'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/data-processing-utility-tools.helper-functions'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB

process filesExist {
    container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

    cpus params.cpus
    memory "${params.mem} GB"

    input:
        val file_names  // file name shall not have spaces
        val expect  // 'exist' for files expected to exist; 'not_exist' for files expected not exist
        path files
        val dependency_flag  // any output from process(es) you'd like to make this process depend on

    script:
        file_name_arg = file_names instanceof List ? file_names.join(" ") : file_names
        """
        if [[ "${expect}" = "exist"  ]]; then
            for f in \$(echo "${file_name_arg}"); do
                if [[ ! -f \$f ]]; then
                    echo "Expected \$f not exists."
                    exit 1
                fi
            done
        elif [[ "${expect}" = "not_exist"  ]]; then
            for f in \$(echo "${file_name_arg}"); do
                if [[ -f \$f ]]; then
                    echo "Unexpected \$f exists."
                    exit 1
                fi
            done
        else
            echo "Second argument must be either 'exist' or 'not_exist'. '${expect}' is supplied."
            exit 1
        fi
        """
}
