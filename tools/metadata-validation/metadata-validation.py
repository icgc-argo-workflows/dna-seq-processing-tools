#!/usr/bin/env python3

import sys
import csv
from argparse import ArgumentParser
import json
import datetime
import copy

def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set) and len(obj) > 1:
        return list(obj)
    if isinstance(obj, set) and len(obj) == 1:
        return obj.pop()
    raise TypeError

def reshape_metadata(input_metadata):
    output_metadata = {}
    output_files = {}
    for key, value in input_metadata.items():
        if not key == 'read_groups':
            output_metadata[key] = value
            continue
        # process readGroups
        for rg in input_metadata['read_groups']:
            output_rg = {}
            for k, v in rg.items():
                if not k == 'files':
                    output_rg[k] = v
                    continue
                for fn in rg['files']:
                    if not fn['name'] in output_files:
                        output_files[fn['name']] = {}
                        output_files[fn['name']]['read_groups'] = []
                    for fk, fv in fn.items():
                        if not fk in output_files[fn['name']]: output_files[fn['name']][fk] = set()
                        output_files[fn['name']][fk].add(fv)
            output_files[fn['name']]['read_groups'].append(output_rg)
    output_metadata['files'] = []
    for key, value in output_files.items():
        output_metadata['files'].append(value)
    return output_metadata


def build_entity_map(entity_tsv, pointer_id):
    with open(entity_tsv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        entity_dict = {}
        for l in reader:
            if not entity_dict.get(l[pointer_id]): entity_dict[l[pointer_id]] = []
            sub_entity_dict = {}
            for field in reader.fieldnames:
                if field in ['read_group_count', 'read_length_r1', 'read_length_r2', 'insert_size', 'size']:
                    sub_entity_dict[field] = int(l.get(field))
                elif l.get(field) in ['True', 'true', 'TRUE']:
                    sub_entity_dict[field] = True
                elif l.get(field) in ['False', 'false', 'FALSE']:
                    sub_entity_dict[field] = False
                else:
                    sub_entity_dict[field] = l.get(field, None)
            entity_dict[l[pointer_id]].append(sub_entity_dict)
    return entity_dict

def get_input_format(file_tsv):
    with open(file_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        input_format = set()
        for l in reader:
            input_format.add(l.get('format').upper())

    if not len(input_format) == 1: sys.exit('\nError: The input files should have the same format.')

    return input_format.pop()

def generate_metadata(args):

    # build {submitter_read_group_id: [files]} map
    rg_file_map = build_entity_map(args.file_tsv, "submitter_read_group_id")

    # build {submitter_sequencing_experiment_id: [readgroups]} map
    exp_rg_map = build_entity_map (args.rg_tsv, "submitter_sequencing_experiment_id")

    # build metadata based on schema
    exp_json_map = build_entity_map(args.exp_tsv, "submitter_sequencing_experiment_id")

    # only permit one experiment input
    if not len(exp_json_map) == 1: sys.exit('\nError: The input should only contain one experiment!')

    exp_id = set()
    exp_json = {}
    for key, val in exp_json_map.items():
        exp_id.add(key)
        exp_json = val[0]
        if not exp_rg_map.get(exp_json.get('submitter_sequencing_experiment_id')):
            sys.exit('\nError: The input experiment.tsv and read_group.tsv have mismatch experiment IDs!')
        exp_json['read_groups'] = exp_rg_map.get(exp_json.get('submitter_sequencing_experiment_id'))
        rg_id = set()
        for rg in exp_json['read_groups']:
            rg_id.add(rg.get('submitter_read_group_id'))
            if not rg_file_map.get(rg.get('submitter_read_group_id')):
                sys.exit('\nError: The input read_group.tsv and file.tsv have mismatch read_group IDs!')
            rg['files'] = rg_file_map.get(rg.get('submitter_read_group_id'))

        # validate read_group ids match across read_group.tsv and file.tsv
        if not rg_id == set(rg_file_map.keys()): sys.exit('\nError: The input read_group.tsv and file.tsv have mismatch read_group IDs!')

    # validate experiment ids match across experiment.tsv and read_group.tsv
    if not exp_id == set(exp_rg_map.keys()): sys.exit('\nError: The input experiment.tsv and read_group.tsv have mismatch experiment IDs!')

    return exp_json

def run_validation(args):

    # generate the exp_json with experiment, read_group, file info
    if args.meta_format == 'tsv':
        input_seq_format = get_input_format(args.file_tsv)
        exp_json = generate_metadata(args)
    else:
        exp_json = args.exp_json

    # generate metadata from exp_json by reshaping the exp_json in condition of the input_format
    # cross check the read_group id for BAM
    if input_seq_format == "BAM":
        metadata = reshape_metadata(exp_json)
        metadata["input_seq_format"] = "BAM"

    elif input_seq_format == "FASTQ":
        metadata = copy.deepcopy(exp_json)
        metadata["input_seq_format"] = "FASTQ"

    else:
        sys.exit('\nError: The input files should have format in BAM or FASTQ.')


    with open(args.metadata_json, 'w') as f:
        f.write(json.dumps(metadata, default=set_default, indent=2))


    # remove files from the exp_json
    # validate the json with the given schema and the data get from the server
    for rg in exp_json['read_groups']:
        rg.pop('files')

    # write the exp_rg.json file for later submission
    with open(args.exp_rg_json, 'w') as f:
        f.write(json.dumps(exp_json, default=set_default, indent=2))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--meta", dest="meta_format", type=str, help="metadata file format",
                        choices=["tsv", "json"],
                        default="tsv")
    parser.add_argument("-j", "--exp_json", dest="exp_json", help="json format file containing readgroup and input BAM/FASTQ file information")
    parser.add_argument("-e", "--exp_tsv", dest="exp_tsv", help="tsv format file containing experiment information")
    parser.add_argument("-r", "--rg_tsv", dest="rg_tsv", help="tsv format file containing readgroup information")
    parser.add_argument("-f", "--file_tsv", dest="file_tsv", help="tsv format file containing BAM/FASTQ input file information")
    parser.add_argument("-o", "--exp_rg_json", dest="exp_rg_json", help="json file containing experiment and read_group information", required=True)
    parser.add_argument("-p", "--metadata_json", dest="metadata_json",
                        help="json file containing experiment, read_group and file information for sequence preprocessing", required=True)
    args = parser.parse_args()

    if (args.meta_format == 'json' and not args.exp_json):
        parser.error('The json input argument requires the --exp-json')
    if (args.meta_format == 'tsv' and (not args.exp_tsv or not args.rg_tsv or not args.file_tsv)):
        parser.error('The tsv input argument requires both --exp-tsv, --rg-tsv and --file-tsv')

    run_validation(args)