#!/usr/bin/env python3

import os
import sys
import csv
import re
from argparse import ArgumentParser
import json
import datetime


TSV_FIELDS = {
    'experiment': [
        'type', 'program_id', 'submitter_sequencing_experiment_id', 'submitter_donor_id', 'gender',
        'submitter_specimen_id', 'tumour_normal_designation', 'specimen_tissue_source', 'submitter_sample_id',
        'sample_type', 'submitter_matched_normal_sample_id', 'sequencing_center', 'platform', 'platform_model',
        'library_strategy', 'sequencing_date', 'read_group_count'
    ],
    'read_group': [
        'type', 'submitter_read_group_id', 'submitter_sequencing_experiment_id', 'platform_unit',
        'is_paired_end', 'file_r1', 'file_r2', 'read_length_r1', 'read_length_r2', 'insert_size', 'sample_barcode', 'library_name'
    ],
    'file': [
        'type', 'name', 'size', 'md5sum', 'path', 'format'
    ]
}


def tsv_confomity_check(ftype, tsv):
    expected_fields = TSV_FIELDS[ftype]

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
    # first let's make sure all entities have type and PK populated
    if metadata_dict.get('type') != 'sequencing_experiment':  # check: e.1
        sys.exit("Error found: type field in experiment TSV must be 'sequencing_experiment', %s found\n" % \
            metadata_dict.get('type'))

    # submitter_sequencing_experiment_id must match regex: e.3
    submitter_sequencing_experiment_id = metadata_dict.get('submitter_sequencing_experiment_id')
    if not re.match(r'^[a-zA-Z0-9_\.\-]+$', submitter_sequencing_experiment_id):
        sys.exit("Error found: invalid submitter_sequencing_experiment_id: '%s' in experiment TSV\n" % \
            submitter_sequencing_experiment_id)

    uniq_files = {}
    for f in metadata_dict['files']:
        if f.get('type') != 'file':  # check: f.1
            sys.exit("Error found: type field in file TSV must be 'file', %s found\n" % \
                f.get('type'))

        file_name = f.get('name', '')  # file name serves as PK in file TSV
        if not re.match(r'^[a-zA-Z0-9]{1}[a-zA-Z0-9_\.\-]*\.(bam|fq|fq\.gz|fq\.bz2|fastq|fastq\.gz|fastq\.bz2)$', file_name):
            sys.exit(
                "Error found: invalid file name %s in file TSV."
                " File name must include only characters [a-zA-Z0-9_\\.-] and start with one of [a-zA-Z0-9]."
                " Accepted file extensions include: .bam, .fq, .fq.gz, .fq.bz2, .fastq, .fastq.gz, .fastq.bz2\n" \
                % file_name)

        if file_name in uniq_files:  # file name duplicated
            sys.exit("Error found: file name %s duplicated in file TSV" % file_name)
        else:
            uniq_files[file_name] = True

    uniq_rgs = {}
    file_to_rgs = {}
    for rg in metadata_dict['read_groups']:
        if rg.get('type') != 'read_group':  # check: g.1
            sys.exit("Error found: type field in read group TSV must be 'read_group', %s found\n" % \
                rg.get('type'))

        # submitter_read_group_id must match regex, and be unique: g.2
        submitter_read_group_id = rg.get('submitter_read_group_id')
        if not re.match(r'^[a-zA-Z0-9_\:\.\-]+$', submitter_read_group_id):
            sys.exit("Error found: invalid submitter_read_group_id: '%s' in read group TSV\n" % \
                submitter_read_group_id)
        if submitter_read_group_id in uniq_rgs:  # read group id not unique
            sys.exit("Error found: submitter_read_group_id: '%s' duplicated in read group TSV\n" % \
                submitter_read_group_id)
        else:
            uniq_rgs[submitter_read_group_id] = True

        # check is_paired_end: g.5
        if rg['is_paired_end'] == 'true':
            rg['is_paired_end'] = True
        elif rg['is_paired_end'] == 'false':
            rg['is_paired_end'] = False
        else:
            sys.exit("Error found: is_paired_end must be 'true' or 'false' for read group '%s' in read group TSV\n" % \
                submitter_read_group_id)

        # check fields reference to file: file_r1 and file_r2
        # file_r1 must be populated regardless is_paired_end or not
        if rg['file_r1'] not in uniq_files:  # file_r1 must be populated with a known file
            sys.exit("Error found: unknown file_r1 %s in read group %s" % \
                (rg['file_r1'], submitter_read_group_id))

        if rg['file_r1'] not in file_to_rgs:
            file_to_rgs[rg['file_r1']] = set([submitter_read_group_id])
        else:
            file_to_rgs[rg['file_r1']].add(submitter_read_group_id)

        if rg['is_paired_end']:  # is paired end
            # file_r2 must be populated with a known file
            if rg['file_r2'] not in uniq_files:
                sys.exit("Error found: unknown file_r2 %s in read group %s\n" % \
                    (rg['file_r2'], submitter_read_group_id))

            # for BAM, file_r1 and file_r2 must be the same
            if rg['file_r1'].endswith('.bam') and rg['file_r1'] != rg['file_r2']:
                sys.exit("Error found: file_r1 must be the same as file_r2 for BAM file."
                    " However in read group %s, file_r1 is %s, file_r2 is %s\n" % \
                    (submitter_read_group_id, rg['file_r1'], rg['file_r2']))

            # for FASTQ, file_r1 and file_r2 must not be the same
            if (not rg['file_r1'].endswith('.bam')) and rg['file_r1'] == rg['file_r2']:
                sys.exit("Error found: file_r1 must not be the same as file_r2 for FASTQ files."
                    " However in read group %s, file_r1 and file_r2 are both %s\n" % \
                    (submitter_read_group_id, rg['file_r1']))

            if rg['file_r2'] not in file_to_rgs:
                file_to_rgs[rg['file_r2']] = set([submitter_read_group_id])
            else:
                file_to_rgs[rg['file_r2']].add(submitter_read_group_id)

        else:  # not paired end
            if rg['file_r2'] != "":  # single end file_r2 must be empty
                sys.exit("Error found: for single end read group %s, file_r2 must be empty, however %s is found\n" % \
                    (submitter_read_group_id, rg['file_r2']))

            rg['file_r2'] = None  # convert from "" to null

    for f in file_to_rgs:
        # this is to check FASTQ file can only be associated with one read group
        if not f.endswith('.bam') and len(file_to_rgs[f]) > 1:
            sys.exit("Error found: FASTQ file must be linked to exactly one read group. However %s is linked to: %s\n" % \
                (f, ', '.join(file_to_rgs[f])))

    # check to make sure all files in file TSV are referenced in read group TSV
    unreferenced_files = set(uniq_files.keys()) - set(file_to_rgs.keys())
    if unreferenced_files:
        sys.exit("Error found: sequencing file(s) not linked to any read group: %s" % \
            ', '.join(sorted(unreferenced_files)))

    # all read groups must have submitter_sequencing_experiment_id with value matching that in experiment: g.3
    for rg in metadata_dict['read_groups']:
        if rg.get('submitter_sequencing_experiment_id') != submitter_sequencing_experiment_id:
            sys.exit("Error found: submitter_sequencing_experiment_id in read group '%s' does not match '%s'\n" % \
                (rg.get('submitter_read_group_id'), submitter_sequencing_experiment_id))

    # read_group_count in experiment TSV must match total number of read groups: g.2
    read_group_count = metadata_dict['read_group_count']
    if not re.match(r'^[1-9]{1}[0-9]*$', read_group_count):
            sys.exit("Error found: read_group_count %s in experiment TSV is not a positive integer\n" % \
                read_group_count)

    metadata_dict['read_group_count'] = int(read_group_count)  # now convert from str to int and compare
    if metadata_dict['read_group_count'] != len(metadata_dict['read_groups']):
        sys.exit("Error found: specified read_group_count %s in experiment TSV does not match number of read groups %s in read group TSV\n" % \
            (metadata_dict['read_group_count'], len(metadata_dict['read_groups'])))


def check_experiment(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    if not re.match(r'^[A-Z]{1}[A-Z0-9\-]+$', metadata_dict.get('program_id')):  # check: e.2
        sys.exit("Error found: invalid program_id in experiment TSV: %s\n" % \
            metadata_dict.get('program_id'))

    tumour_normal_designation = metadata_dict.get('tumour_normal_designation')
    # very coarse check here, leave it to SONG to perform more precise check with permissible vaules
    if 'normal' not in tumour_normal_designation.lower() and 'tumour' not in tumour_normal_designation.lower():
        sys.exit("Error found: invalid tumour_normal_designation %s in experiment TSV\n" % 
            tumour_normal_designation)

    submitter_matched_normal_sample_id = metadata_dict.get('submitter_matched_normal_sample_id')
    if 'normal' not in tumour_normal_designation.lower():
        # for non-normal specimen, submitter_matched_normal_sample_id must be provided
        if submitter_matched_normal_sample_id == "":
            sys.exit("Error found: when tumour_normal_designation of the current specimen is not normal: %s, "
                "submitter_matched_normal_sample_id must not be empty in experiment TSV\n" % \
                tumour_normal_designation)

        if submitter_matched_normal_sample_id == metadata_dict.get('submitter_sample_id'):
            sys.exit("Error found: submitter_matched_normal_sample_id %s must not be the same as submitter_sample_id in experiment TSV\n" % \
                submitter_matched_normal_sample_id)

        if not re.match(r'^[a-zA-Z0-9_\.\-]+$', submitter_matched_normal_sample_id):
            sys.exit("Error found: invalid submitter_matched_normal_sample_id: '%s' in experiment TSV\n" % \
                submitter_matched_normal_sample_id)
    else:
        # otherwise, submitter_matched_normal_sample_id must be empty
        if submitter_matched_normal_sample_id != "":
            sys.exit("Error found: when tumour_normal_designation of the current specimen is normal: %s, "
                "submitter_matched_normal_sample_id must be empty in experiment TSV. However, it's populated with %s\n" % \
                (tumour_normal_designation, submitter_matched_normal_sample_id))

        metadata_dict['submitter_matched_normal_sample_id'] = None  # now we can convert "" to null


def check_read_groups(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    platform_units = {}
    for rg in metadata_dict.get('read_groups'):
        read_length_r1 = rg.get('read_length_r1')
        if not re.match(r'^[1-9]{1}[0-9]*$', read_length_r1):
            sys.exit("Error found: read_length_r1 %s in read group %s is not a positive integer\n" % \
                (read_length_r1, rg.get('submitter_read_group_id')))

        rg['read_length_r1'] = int(read_length_r1)

        read_length_r2 = rg.get('read_length_r2')
        insert_size = rg.get('insert_size')
        if rg.get('is_paired_end'):  # paired end
            if not re.match(r'^[1-9]{1}[0-9]*$', read_length_r2):
                sys.exit("Error found: read_length_r2 %s in read group %s is not a positive integer\n" % \
                    (read_length_r2, rg.get('submitter_read_group_id')))

            rg['read_length_r2'] = int(read_length_r2)

            if not re.match(r'^[1-9]{1}[0-9]*$', insert_size):
                sys.exit("Error found: insert_size %s in read group %s is not a positive integer\n" % \
                    (insert_size, rg.get('submitter_read_group_id')))

            rg['insert_size'] = int(insert_size)

        else:  # single end
            if read_length_r2 != "":
                sys.exit("Error found: read_length_r2 in read group %s must not be populated for single end sequencing, however %s is found\n" % \
                    (rg.get('submitter_read_group_id'), read_length_r2))

            rg['read_length_r2'] = None

            if insert_size != "":
                sys.exit("Error found: insert_size in read group %s must not be populated for single end sequencing, however %s is found\n" % \
                    (rg.get('submitter_read_group_id'), insert_size))

            rg['insert_size'] = None

        platform_unit = rg.get('platform_unit')
        if not re.match(r'^[a-zA-Z0-9]{1}[a-zA-Z0-9_\:\.\-]*$', platform_unit):
            sys.exit("Error found: invalid platform_unit: '%s' in experiment TSV\n" % \
                platform_unit)

        if platform_unit in platform_units:  # duplicated
            sys.exit("Error found: duplicated platform_unit %s in read group TSV\n" % platform_unit)
        else:
            platform_units[platform_unit] = True


def check_files(metadata_dict):
    # check individual fields, it's also a good time to convert value to correct type
    for f in metadata_dict.get('files'):
        file_name = f.get('name')
        if file_name.endswith('.bam'):  # bam file
            if f['format'] != 'BAM':
                sys.exit("Error found: for BAM file %s, format field must be BAM. However %s is found\n" % \
                    (file_name, f['format']))
        else:  # fastq file, we have verified all other acceptable files are FASTQ with different file extension
            if f['format'] != 'FASTQ':
                sys.exit("Error found: for FASTQ file %s, format field must be FASTQ. However %s is found\n" % \
                    (file_name, f['format']))

        file_size = f.get('size')
        if not re.match(r'^[1-9]{1}[0-9]*$', file_size):
            sys.exit("Error found: size %s for file %s is not a positive integer\n" % \
                (file_size, file_name))

        f['size'] = int(file_size)

        md5sum = f.get('md5sum')
        if not re.match(r'^[a-f0-9]{32}$', md5sum):
            sys.exit("Error found: md5sum %s for file %s does not seem like an md5sum. Hint md5sum value must be in lower case\n" % \
                (md5sum, file_name))

        path_to_file = f.get('path')
        if path_to_file.startswith('score://'):
            if not re.match(r'^score:\/\/(collab|aws)\/[a-zA-Z0-9\-]+\/[0-9a-f]{8}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{12}$', path_to_file):
                sys.exit("Error found: path %s for file %s looks like a SCORE URL but it's not well-formed\n" % \
                    (path_to_file, file_name))
        else:
            if file_name != os.path.basename(path_to_file):
                sys.exit("Error found: basename of path %s must match file name %s" % \
                    (path_to_file, file_name))


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