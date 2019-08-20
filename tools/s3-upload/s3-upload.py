#!/usr/bin/env python3

import os
import sys
import json
from argparse import ArgumentParser
import subprocess


def run_command(cmd):
    return subprocess.run(cmd, capture_output=True, shell=True)


def object_exists(endpoint_url, object_key):
    ret = run_command('aws --endpoint-url %s s3 ls %s' % (endpoint_url, object_key))
    if ret.returncode == 0:
        return True
    else:
        return False


def copy_credential_file(credential_file):
    run_command('mkdir ~/.aws')
    run_command('cp %s ~/.aws/credentials' % credential_file)


def get_payload(filename_to_file, payload_jsons):
    for pf in payload_jsons:
        with open(pf) as f:
            payload = json.load(f)

        match = True
        if len(payload['files']) != len(filename_to_file):
            match = False
            break
        for f in payload['files']:
            fname = payload['files'][f]['name']
            if fname not in filename_to_file:
                match = False
                break
        if match:
            return payload

    return


def main(args):
    filename_to_file = {}
    sfile_filename_to_file = {}
    f = args.upload_file
    filename_to_file[os.path.basename(f)] = f
    # add secondary files if any
    for sf in ('.bai', '.crai'):
        sfile = "%s%s" % (f, sf)
        if os.path.isfile(sfile):
            sfile_filename_to_file[os.path.basename(sfile)] = sfile

    if sfile_filename_to_file:
        filename_to_file = {**filename_to_file, **sfile_filename_to_file}

    copy_credential_file(args.s3_credential_file)

    # find the corresponding payload from provided list of payloads
    payload = get_payload(filename_to_file, args.payload_jsons)
    if not payload:
        sys.exit("Could not find payload matching the files to be uploaded. Files to be uploaded: %s" %
            ", ".join(list(filename_to_file.keys()))
        )

    if args.bundle_type != payload['type']:
        sys.exit("Specified bundle type '%s' does not match what's defined in the payload: '%s'" % \
            (args.bundle_type, payload['type']))

    bundle_id = payload['id']

    path_prefix = "PCAWG2/%s/%s/%s/%s.%s" % (
                                                payload['info']['library_strategy'],
                                                payload['program'],
                                                payload['info']['donor_submitter_id'],
                                                payload['info']['sample_submitter_id'],
                                                'normal' if 'normal' in payload['info']['specimen_type'].lower() else 'tumour'
                                            )

    if args.bundle_type == 'lane_seq_submission':
        read_group_submitter_id = payload['inputs']['read_group_submitter_id']

        """ don't do this for now
        payload_object_key = "%s/lane_seq_submission/%s/%s.json" % (
            path_prefix,
            read_group_submitter_id,
            bundle_id)
        if not object_exists(args.endpoint_url, 's3://%s/%s' % (args.bucket_name, payload_object_key)):
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
            p = run_command('aws --endpoint-url %s s3 cp %s s3://%s/%s' % (
                    args.endpoint_url,
                    file_to_upload,
                    args.bucket_name,
                    object_key)
                )

            if p.returncode != 0:
                sys.exit('Object upload failed: %s; err: %s' % (object_key, p.stderr))

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
        if not object_exists(args.endpoint_url, 's3://%s/%s' % (args.bucket_name, payload_object_key)):
            sys.exit('Not able to access object store, or payload object does not exist: s3://%s/%s' % (args.bucket_name, payload_object_key))
        """

        for object in payload['files']:
            object_id = payload['files'][object]['object_id']
            filename = payload['files'][object]['name']
            object_key = "%s/dna_alignment/%s/%s/%s" % (path_prefix,
                                                        bam_cram,
                                                        bundle_id,
                                                        object_id)

            file_to_upload = filename_to_file[filename]
            p = run_command('aws --endpoint-url %s s3 cp %s s3://%s/%s' % (
                    args.endpoint_url,
                    file_to_upload,
                    args.bucket_name,
                    object_key)
                )

            if p.returncode != 0:
                sys.exit('Object upload failed: %s; err: %s' % (object_key, p.stderr))

    else:
        sys.exit('Unknown or unimplemented bundle_type: %s' % args.bundle_type)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--endpoint-url", dest="endpoint_url")
    parser.add_argument("-b", "--bucket-name", dest="bucket_name")
    parser.add_argument("-t", "--bundle-type", dest="bundle_type")
    parser.add_argument("-p", "--payload-jsons", dest="payload_jsons", type=str, nargs='+')
    parser.add_argument("-c", "--s3-credential-file", dest="s3_credential_file")
    parser.add_argument("-f", "--upload-file", dest="upload_file")
    args = parser.parse_args()

    main(args)
