#!/usr/bin/env python3

import sys
import csv
from argparse import ArgumentParser
import json
import datetime


def set_default(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    if isinstance(obj, set) and len(obj) > 1:
        return list(obj)
    if isinstance(obj, set) and len(obj) == 1:
        return obj.pop()
    raise TypeError


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

    # generate the metadata with experiment, read_group, file tsv info
    metadata = generate_metadata(args)

    # additional checks



    # write the metadata.json as output
    with open('metadata.json', 'w') as f:
        f.write(json.dumps(metadata, default=set_default, indent=2))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--exp_tsv", dest="exp_tsv", help="tsv format file containing experiment information")
    parser.add_argument("-r", "--rg_tsv", dest="rg_tsv", help="tsv format file containing readgroup information")
    parser.add_argument("-f", "--file_tsv", dest="file_tsv", help="tsv format file containing BAM/FASTQ input file information")
    args = parser.parse_args()

    run_validation(args)