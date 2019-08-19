#!/usr/bin/env python3

import os
import sys
import json
import subprocess
from argparse import ArgumentParser
import uuid

"""
Major steps:
- convert input Seq to unaligned BAM for each read group
"""


def run_cmd(cmd):
    stdout, stderr, p, success = '', '', None, True
    try:
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
        p.communicate()
    except Exception as e:
        print('Execution failed: %s' % e)
        success = False

    if p and p.returncode != 0:
        print('Execution failed, none zero code returned. %s' % p.returncode)
        success = False

    if not success:
        sys.exit(p.returncode if p.returncode else 1)

    return


def get_uuid5(bid, fid):
    uuid5 = str(uuid.uuid5(uuid.UUID("6ba7b810-9dad-11d1-80b4-00c04fd430c8"), "%s/%s" % (bid, fid)))
    return uuid5

def main(args):

    with open(args.metadata, 'r') as f:
        metadata = json.load(f)

    with open(args.payload, 'r') as f:
        payload = json.load(f)

    # generate bundle_id
    payload['id'] = str(uuid.uuid4())

    # generate info field
    payload['info'] = {
        "library_strategy": metadata.get("library_strategy", None),
        "program_id": metadata.get("program", None),
        "submitter_donor_id": metadata.get("submitter_donor_id", None),
        "sample_submitter_id": metadata.get("sample_submitter_id", None),
        "specimen_type": metadata.get("specimen_type", None)
    }

    # generate object_id
    for key, val in payload['files'].items():
        val['object_id'] = get_uuid5(payload['id'], val['name'])

    payload_fname = ".".join([payload['id'], 'json'])
    with open(payload_fname, 'w') as f:
        f.write(json.dumps(payload))

    # payload bucket basepath
    specimen_type = 'normal' if 'normal' in metadata.get("specimen_type", '').lower() else 'tumour'
    payload_bucket_basepath = os.path.join(args.bucket_name, 'PCAWG2',
                                        payload['info']['library_strategy'],
                                        payload['info']['program_id'],
                                        payload['info']['donor_submitter_id'],
                                        payload['info']['sample_submitter_id']+'.'+specimen_type,
                                        payload['type'])

    if payload['type'] in ['sequencing_experiment', 'dna_alignment_qc']:
        payload_object = os.path.join(payload_bucket_basepath, payload_fname)
    elif payload['type'] in ['lane_seq_submission', 'lane_seq_qc']:
        payload_object = os.path.join(payload_bucket_basepath, payload['inputs']['read_group_submitter_id'], payload_fname)

    elif payload['type'] in ['dna_alignment']:
        alignment_type = "bam" if payload['files']['aligned_seq']['name'].endswith('bam') else "cram"
        payload_object = os.path.join(payload_bucket_basepath, alignment_type, payload_fname)

    else:
        sys.exit('Unknown payload type!')

    cmd = "aws --endpoint-url %s s3 cp %s s3://%s" % (args.endpoint_url, payload_fname, payload_object)
    run_cmd(cmd)



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--endpoint-url", dest="endpoint_url")
    parser.add_argument("-b", "--bucket-name", dest="bucket_name")
    parser.add_argument("-m", "--metadata", dest="metadata", help="Metadata file contains the information user submit")
    parser.add_argument("-p", "--payload", dest="payload", help="Payload file")
    args = parser.parse_args()

    main(args)
