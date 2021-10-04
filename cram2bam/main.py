#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
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
"""

import os
import sys
import argparse
import subprocess


def main():
    """
    Python implementation of tool: cram2bam

    This is auto-generated Python code, please update as needed!
    """

    parser = argparse.ArgumentParser(description='Tool: cram2bam')
    parser.add_argument('-i', '--input-file', dest='input_file', type=str,
                        help='Input file', required=True)
    parser.add_argument('-o', '--output-dir', dest='output_dir', type=str,
                        help='Output directory', required=True)
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.input_file)

    if not os.path.isdir(args.output_dir):
        sys.exit('Error: specified output dir %s does not exist or is not accessible!' % args.output_dir)

    subprocess.run(f"cp {args.input_file} {args.output_dir}/", shell=True, check=True)


if __name__ == "__main__":
    main()

