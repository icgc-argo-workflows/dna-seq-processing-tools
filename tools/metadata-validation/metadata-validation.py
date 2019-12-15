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


def check_experiment(exp_tsv):
    exp_id = set()
    with open(exp_tsv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # only permit one experiment input
        for l in reader:
            experiment_dict = {}
            for field in reader.fieldnames:
                if field not in ['submitter_matched_normal_sample_id', 'sequencing_center', 'sequencing_date', 'type'] and l.get(field) is None:
                    sys.exit('\nError: Missing required field: %s in experiment.tsv!' % field)
                elif field in ['read_group_count']:
                    experiment_dict[field] = int(l.get(field))
                else:
                    experiment_dict[field] = l.get(field)
            exp_id.add(l.get('submitter_sequencing_experiment_id'))

        if not len(exp_id) == 1: sys.exit('\nError: The input should only contain one experiment!')

    return (experiment_dict, exp_id)

def check_files(file_tsv):
    rg_file = {}
    files_dict = {}
    with open(file_tsv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for l in reader:
            for field in reader.fieldnames:
                if l.get(field) is None: sys.exit('\nError: Missing required field: %s in file.tsv!' % field)
            if not files_dict.get(l["name"]): files_dict[l["name"]] = {}
            for field in ['name', 'size', 'md5sum', 'path', 'format']:
                if not files_dict[l["name"]].get(field): files_dict[l["name"]][field] = set()
                value = int(l.get(field)) if field == 'size' else l.get(field)
                files_dict[l["name"]][field].add(value)

            if not rg_file.get(l.get('submitter_read_group_id')):
                rg_file[l.get('submitter_read_group_id')] = {}
            if not l.get('r1_r2') in ['r1/r2', 'r1', 'r2']:
                sys.exit('\nError: Invalid value of r1_r2 in file.tsv!')
            elif l.get('r1_r2') == 'r1/r2':
                rg_file[l.get('submitter_read_group_id')]['file_r1'] = l.get('name')
                rg_file[l.get('submitter_read_group_id')]['file_r2'] = l.get('name')
            elif l.get('r1_r2') == 'r1':
                rg_file[l.get('submitter_read_group_id')]['file_r1'] = l.get('name')
            else:
                rg_file[l.get('submitter_read_group_id')]['file_r2'] = l.get('name')

    files = []
    for value in files_dict.values():
        for key in value.keys():
            if isinstance(value[key], set) and len(value[key]) > 1:
                sys.exit('\nError: Inconsistent values of field: %s in file.tsv!' % key)
        files.append(value)

    return (files, rg_file)

def check_read_group(rg_tsv, rg_file):
    read_group_dict = {}
    experiment_id = set()
    read_group_id = set()
    with open(rg_tsv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for l in reader:
            for field in reader.fieldnames:
                if field not in ['read_length_r2', 'insert_size', 'sample_barcode'] and l.get(field) is None:
                    sys.exit('\nError: Missing value of field: %s in read_group.tsv!' % field)
            experiment_id.add(l.get('submitter_sequencing_experiment_id'))
            read_group_id.add(l.get('submitter_read_group_id'))

            if not read_group_dict.get(l['submitter_read_group_id']): read_group_dict[l['submitter_read_group_id']] = {}
            read_group_dict.get(l['submitter_read_group_id']).update(rg_file.get(l['submitter_read_group_id']))
            for field in reader.fieldnames:
                if field in ['type', 'submitter_sequencing_experiment_id']: continue
                if not read_group_dict.get(l['submitter_read_group_id']).get(field):
                    read_group_dict.get(l['submitter_read_group_id'])[field] = set()
                if field in ['read_length_r1', 'read_length_r2', 'insert_size']:
                    value = int(l.get(field))
                elif l.get(field) in ['True', 'true', 'TRUE']:
                    value = True
                elif l.get(field) in ['False', 'false', 'FALSE']:
                    value = False
                else:
                    value = l.get(field)
                read_group_dict.get(l['submitter_read_group_id'])[field].add(value)

    read_groups = []
    for value in read_group_dict.values():
        for key in value.keys():
            if isinstance(value[key], set) and len(value[key]) > 1:
                sys.exit('\nError: Inconsistent values of field: %s in read_group.tsv!' % key)
        read_groups.append(value)

    return (read_groups, experiment_id, read_group_id)

def run_validation(args):

    # check file.tsv
    # input file.tsv
    # output files dict, rg_file_map={rg_id: file_name_r1, file_name_r2}
    files, rg_file = check_files(args.file_tsv)

    # check experiment.tsv
    # input experiment.tsv
    # output experiment dict, exp_id
    experiment_dict, exp_id = check_experiment(args.exp_tsv)

    # check read_group.tsv
    # input read_group.tsv, rg_file_map
    # output read_groups dict
    read_groups, experiment_id, read_group_id = check_read_group(args.rg_tsv, rg_file)

    # additional checks
    if not set(rg_file.keys()) == read_group_id:
        sys.exit('\nError: The input read_group.tsv and file.tsv have mismatch read_group IDs!')

    if not exp_id == experiment_id:
        sys.exit('\nError: The input experiment.tsv and read_group.tsv have mismatch experiment IDs!')

    if not experiment_dict['read_group_count'] == len(read_groups):
        sys.exit('\nError: The input read_group.tsv and experiment.tsv has mismatch read_group_count')

    # generate the metadata with experiment, read_groups, files dict
    metadata = {}
    metadata.update(experiment_dict)
    metadata['files'] = files
    metadata['read_groups'] = read_groups

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