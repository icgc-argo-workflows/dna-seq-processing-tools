#!/usr/bin/env python3

import os
import sys
import json
from argparse import ArgumentParser
import subprocess


def run_command(cmd):
    p = subprocess.run(cmd, capture_output=True, shell=True)
    return (p.returncode, p.stdout, p.stderr)


def object_exist(endpoint_url, object_key):
    ret = run_command('aws --endpoint-url %s s3 ls %s' % (endpoint_url, object_key))
    if ret[0] == 0:
        return True
    else:
        return False

def copy_credential_file(credential_file):
    run_command('mkdir ~/.aws')
    run_command('cp %s ~/.aws/credentials' % credential_file)


def main(args):
    filename_to_file = {}
    for f in args.upload_files:
        filename_to_file[os.path.basename(f)] = f

    copy_credential_file(args.s3_credential_file)

    with open(args.metadata_json) as f:
        metadata = json.load(f)

    path_prefix = "PCAWG2/%s/%s/%s/%s.%s" % (
                                                metadata['library_strategy'],
                                                metadata['program'],
                                                metadata['donor_submitter_id'],
                                                metadata['sample_submitter_id'],
                                                'normal' if 'normal' in metadata['specimen_type'].lower() else 'tumour'
                                            )

    with open(args.payload_json) as f:
        payload = json.load(f)

    bundle_id = payload['id']

    if args.bundle_type == 'lane_seq_submission':
        read_group_submitter_id = payload['inputs']['read_group_submitter_id']

        """ don't do this for now
        payload_object_key = "%s/lane_seq_submission/%s/%s.json" % (
            path_prefix,
            read_group_submitter_id,
            bundle_id)
        if not object_exist(args.endpoint_url, 's3://%s/%s' % (args.bucket_name, payload_object_key)):
            sys.exit('Not able to access object store, or payload object does not exist: s3://%s/%s' % (args.bucket_name, payload_object_key))
        """

        for object in payload['files']:
            object_id = payload['files'][object]['object_id']
            filename = payload['files'][object]['name']
            object_key = "%s/lane_seq_submission/%s/%s/%s" % (path_prefix,
                                                                read_group_submitter_id,
                                                                bundle_id,
                                                                object_id)

            file_to_upload = filename_to_file[filename]
            run_command('aws --endpoint-url %s s3 cp %s s3://%s/%s' % (
                    args.endpoint_url,
                    file_to_upload,
                    args.bucket_name,
                    object_key)
                )

    elif args.bundle_type == 'dna_alignment':
        bam_cram = 'bam'
        for f in filename_to_file:
            if f.endswith('.cram'):
                bam_cram = 'cram'
                break

        """
        payload_object_key = "%s/dna_alignment/%s/%s.json" % (
            path_prefix,
            bam_cram,
            bundle_id)
        if not object_exist(args.endpoint_url, 's3://%s/%s' % (args.bucket_name, payload_object_key)):
            sys.exit('Not able to access object store, or payload object does not exist: s3://%s/%s' % (args.bucket_name, payload_object_key))
        """

        for object in payload['files']:
            object_id = payload['files'][object]['object_id']
            filename = payload['files'][object]['name']
            object_key = "%s/dna_alignment/%s/%s" % (path_prefix,
                                                        bundle_id,
                                                        object_id)

            file_to_upload = filename_to_file[filename]
            run_command('aws --endpoint-url %s s3 cp %s s3://%s/%s' % (
                    args.endpoint_url,
                    file_to_upload,
                    args.bucket_name,
                    object_key)
                )

    else:
        sys.exit('Unknown or unimplemented bundle_type: %s' % args.bundle_type)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--s3-endpoint-url", dest="endpoint_url")
    parser.add_argument("-b", "--bucket-name", dest="bucket_name")
    parser.add_argument("-t", "--bundle-type", dest="bundle_type")
    parser.add_argument("-m", "--metadata-json", dest="metadata_json")
    parser.add_argument("-p", "--payload-json", dest="payload_json")
    parser.add_argument("-c", "--s3-credential-file", dest="s3_credential_file")
    parser.add_argument("-f", "--upload-files", dest="upload_files", type=str, nargs='+')
    args = parser.parse_args()

    main(args)
