#!/usr/bin/env python3

import sys
import csv
import re
from argparse import ArgumentParser
import json
import datetime


tsv_fields = {
    'experiment': [
        'type', 'program_id', 'submitter_sequencing_experiment_id', 'submitter_donor_id', 'gender',
        'submitter_specimen_id', 'tumour_normal_designation', 'specimen_tissue_source', 'submitter_sample_id',
        'sample_type', 'submitter_matched_normal_sample_id', 'sequencing_center', 'platform', 'platform_model',
        'library_strategy', 'sequencing_date', 'read_group_count'
    ],
    'read_group': [
        'type', 'submitter_read_group_id', 'submitter_sequencing_experiment_id', 'platform_unit',
        'is_paired_end', 'read_length_r1', 'read_length_r2', 'insert_size', 'sample_barcode', 'library_name'
    ],
    'file': [
        'type', 'name',	'size', 'submitter_read_group_id', 'md5sum', 'path', 'format', 'r1_r2'
    ]
}


def tsv_confomity_check(ftype, tsv):
    expected_fields = tsv_fields[ftype]

    header_processed = False
    with open(tsv, 'r') as t:
        uniq_row = {}
        for l in t:
            l = l.rstrip('\n').rstrip('\r')  # remove trailing newline, remove windows `\r` (just in case)
            if not header_processed:  # it's header
                fields = l.split('\t')
                if len(fields) != len(set(fields)):
                    sys.exit("Error found: Field duplicated in input TSV: %s, offending header: %s\n" % (tsv, l))

                missed_fields = set(expected_fields) - set(fields)
                if missed_fields:  # missing fields
                    sys.exit("Error found: Field missing in input TSV: %s, offending header: %s. Missed field(s): %s\n" % \
                        (tsv, l, ', '.join(missed_fields)))

                unexpected_fields = set(fields) - set(expected_fields)
                if unexpected_fields:  # unexpected fields
                    sys.exit("Error found: Unexpected field in input TSV: %s, offending header: %s. Unexpected field(s): %s\n" % \
                        (tsv, l, ', '.join(unexpected_fields)))

                header_processed = True

            else:  # it's data row
                # at this point we only check whether number of values matches number of expected fields and uniqueness check,
                # later steps will perform more sophisticated content check
                values = l.split('\t')
                if len(expected_fields) != len(values):
                    sys.exit("Error found: number of fields: %s does not match expected: %s, offending data row: %s\n" % \
                        (len(values), len(expected_fields), l))

                if l in uniq_row:
                    sys.exit("Error found: data row repeated in file: %s, offending data row: %s\n" % (tsv, l))
                else:
                    uniq_row[l] = True


def load_all_tsvs(exp_tsv, rg_tsv, file_tsv):
    metadata_dict = {}
    with open(exp_tsv, 'r') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
        if len(rows) != 1:
            sys.exit("Error found: experiment TSV expects exactly one data row, offending file: %s has %s row(s)\n" % \
                (exp_tsv, len(rows)))
        metadata_dict.update(rows[0])

    with open(rg_tsv, 'r') as f:
        metadata_dict['read_groups'] = list(csv.DictReader(f, delimiter='\t'))
        if len(metadata_dict['read_groups']) == 0:
            sys.exit("Error found: read group TSV does not contain any read group information\n")

    with open(file_tsv, 'r') as f:
        metadata_dict['files'] = list(csv.DictReader(f, delimiter='\t'))
        if len(metadata_dict['files']) == 0:
            sys.exit("Error found: file TSV does not contain any file information\n")

    return metadata_dict


def check_relationships(metadata_dict):
    # first let's make sure all entities have PK populated
    # submitter_sequencing_experiment_id must match regex: e.3
    submitter_sequencing_experiment_id = metadata_dict.get('submitter_sequencing_experiment_id')
    if not re.match(r'[a-zA-Z0-9_\.\-]+', submitter_sequencing_experiment_id):
        sys.exit("Error found: invalid submitter_sequencing_experiment_id: '%s' in experiment TSV\n" % \
            submitter_sequencing_experiment_id)

    # submitter_read_group_id must match regex, unique: g.2
    uniq_rg = {}
    for rg in metadata_dict['read_groups']:
        submitter_read_group_id = rg.get('submitter_read_group_id')
        if not re.match(r'[a-zA-Z0-9_\:\.\-]+', submitter_read_group_id):
            sys.exit("Error found: invalid submitter_read_group_id: '%s' in read group TSV\n" % \
                submitter_read_group_id)
        if submitter_read_group_id in uniq_rg:  # read group id not unique
            sys.exit("Error found: submitter_read_group_id: '%s' duplicated in read group TSV\n" % \
                submitter_read_group_id)
        else:
            uniq_rg[submitter_read_group_id] = True

    # submitter_read_group_id in file TSV must exist in read group TSV: f.2
    for f in metadata_dict['files']:
        submitter_read_group_id = f.get('submitter_read_group_id')
        if not re.match(r'[a-zA-Z0-9_\:\.\-]+', submitter_read_group_id):
            sys.exit("Error found: invalid submitter_read_group_id: '%s' in file TSV\n" % \
                submitter_read_group_id)
        if not submitter_read_group_id in uniq_rg:
            sys.exit("Error found: submitter_read_group_id: '%s' in file TSV does not exist in read group TSV\n" % \
                submitter_read_group_id)

    # all read groups must have this field with value matching that in experiment
    for rg in metadata_dict['read_groups']:
        if rg.get('submitter_sequencing_experiment_id') != submitter_sequencing_experiment_id:
            sys.exit("Error found: submitter_sequencing_experiment_id in read group '%s' does not match '%s'\n" % \
                (rg.get('submitter_read_group_id'), submitter_sequencing_experiment_id))

    # read_group_count in experiment TSV must match total number of read groups: g.2
    read_group_count = metadata_dict['read_group_count']
    if not re.match(r'[0-9]+', read_group_count):
            sys.exit("Error found: read_group_count %s in experiment TSV is not a positive integer\n" % \
                read_group_count)

    # now convert from str to int
    metadata_dict['read_group_count'] = int(read_group_count)

    if metadata_dict['read_group_count'] != len(metadata_dict['read_groups']):
        sys.exit("Error found: specified read_group_count %s in experiment TSV does not match number of read groups %s in read group TSV\n" % \
            (metadata_dict['read_group_count'], len(metadata_dict['read_groups'])))


def check_experiment(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    pass


def check_read_groups(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    pass


def check_files(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    pass


def run_validation(args):
    # fistly TSV format conformity check, if not well-formed no point to continue
    tsv_confomity_check('experiment', args.exp_tsv)
    tsv_confomity_check('read_group', args.rg_tsv)
    tsv_confomity_check('file', args.file_tsv)

    # all TSV are well-formed, let's load them
    metadata_dict = load_all_tsvs(args.exp_tsv, args.rg_tsv, args.file_tsv)

    # check relationships
    check_relationships(metadata_dict)

    # check experiment
    check_experiment(metadata_dict)

    # check read groups
    check_read_groups(metadata_dict)

    # check files
    check_files(metadata_dict)

    # write the metadata.json as output
    with open('metadata.json', 'w') as f:
        f.write(json.dumps(metadata_dict, indent=2) + "\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-e", "--exp_tsv", dest="exp_tsv", help="tsv format file containing experiment information")
    parser.add_argument("-r", "--rg_tsv", dest="rg_tsv", help="tsv format file containing readgroup information")
    parser.add_argument("-f", "--file_tsv", dest="file_tsv", help="tsv format file containing BAM/FASTQ input file information")
    args = parser.parse_args()

    run_validation(args)