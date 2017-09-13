#!/usr/bin/env python
import os
from os.path import join, abspath, dirname
from subprocess import check_call
import re
import shutil
import argparse
import boto3
from logger import logger


def split_s3_path(path):
    bucket, key = re.match(r's3://([^/]*)/(.*)', path).groups()
    base = key.rsplit('/', 1)[1]
    return bucket, key, base


def input_file(path):
    if not path.startswith('s3://'):
        return path
    else:
        client = boto3.client('s3', 'us-west-2')
        transfer = boto3.s3.transfer.S3Transfer(client)
        bucket, key, base = split_s3_path(path)
        dlpath = join(os.environ.get('TMPDIR'), base)
        transfer.download_file(bucket, key, dlpath)
        return dlpath


def output_file(in_path, out_path):
    if out_path.startswith('s3://'):
        client = boto3.client('s3', 'us-west-2')
        transfer = boto3.s3.transfer.S3Transfer(client)
        bucket, key, _ = split_s3_path(out_path)
        transfer.upload_file(in_path, bucket, key, extra_args={'ServerSideEncryption': 'AES256'})
    else:
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        shutil.copy(in_path, out_path)


if __name__ == '__main__':
    logger.info('Xie Chao\'s HLA typing algorithm')
    parser = argparse.ArgumentParser(description='HLA typing')
    parser.add_argument('--sample_id', type=str, required = True, help='Sample ID/Name')
    parser.add_argument('--input_bam_path', type=str, required = True, help='Input file')
    parser.add_argument('--output_path', type=str, required = True, help='Output directory')
    parser.add_argument('--delete', action="store_true",
                        help='Delete all intermediate files')
    parser.add_argument('--full', action="store_true",
                        help='Run full-digit resolution, default: 4-digit')
    args, _ = parser.parse_known_args()

    logger.info('Sample_id: {} Input file: {}'.format(args.sample_id, args.input_bam_path))
    out_local_path = join('hla-' + args.sample_id, args.sample_id + '.json')
    bin_path = join(dirname(abspath(__file__)), 'typer.sh')
    bin_args = [bin_path, args.input_bam_path, args.sample_id]
    if args.delete:
        bin_args += ['delete']
    if args.full:
        bin_args += ['full']

    check_call(bin_args)

    out_final_path = join(args.output_path, 'report-'+args.sample_id+'-hla.json')
    output_file(out_local_path, out_final_path)
    done_path = join('hla-' + args.sample_id, '_SUCCESS')
    open(done_path, 'a').close()
    done_final_path = join(args.output_path, '_SUCCESS')
    output_file(done_path, done_final_path)

    logger.info('Successfully wrote output file')
    logger.info('HLA typing: shutting down.')
