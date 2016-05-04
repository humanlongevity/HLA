#!/usr/bin/env python
import os
import sys
import json
import glob
import shutil
import filecmp
import argparse
import subprocess

INPUTKEYS = 'sample_id input_bam_path'.split()


def make_cli_args(inpaths):
    return ' '.join('--'+k+"='"+v+"'" for k, v in inpaths)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run integration tests')
    parser.add_argument('--bind-source', action='store_true', help='Whether to bind src code into Docker at runtime')
    parser.add_argument('--image', type=str, help='Docker image name')
    args = parser.parse_args()
    print(args.bind_source)
    print(args.image)
    try:
        shutil.rmtree('_candidates')
    except:
        pass
    os.makedirs('_candidates')
    for injson in glob.glob('inputs/*.json'):
        print("")
        print("### Generating output for input manifest {} ###".format(injson))
        inpaths = json.load(open(injson))
        assert set(inpaths.keys()).issubset(set(INPUTKEYS)), "Your input file must contain keys of {} and nothing else"\
            .format(INPUTKEYS)
        comparison_name = os.path.basename(injson).rsplit('.', 1)[0]
        os.makedirs('_candidates/{}'.format(comparison_name))

        curdir = os.path.dirname(os.path.abspath(__file__))
        bindir = os.path.join(curdir, '..', 'bin')
        datadir = os.path.join(curdir, '..', 'data')
        # From the input paths, we now construct the docker bind paths and cli
        # args
        bind_mounts = []
        cli_args = []
        for key, path in inpaths.iteritems():
            if not path.startswith('s3://') and key.endswith('path'):
                full_path = os.path.join(curdir, 'inputs', path)
                assert os.path.exists(full_path), "Input file doesn't exist at {}".format(full_path)
                bind_mounts.append('-v {}:/in/{}'.format(full_path, path))
                cli_args.append('--' + key + "='" + "/in/{}'".format(path))
            else:
                cli_args.append('--' + key + "='" + path + "'")

        if args.bind_source:
            bind_mounts.append('-v {}:/opt/bin'.format(bindir))
            bind_mounts.append('-v {}:/opt/data'.format(datadir))

        bind_mounts.append('-v {}:/out'.format(os.path.join(curdir, '_candidates', comparison_name)))
        cli_args.append("--output_path=/out")
        cmd = 'docker run --rm -ti \
            -v $HOME/.aws:/root/.aws:ro \
            -v $TMPDIR:/scratch \
            -e TMPDIR=/scratch \
            {} \
            --net=host \
            {} {}'.format(' '.join(bind_mounts), args.image, ' '.join(cli_args))
        print("Calling out to {}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        out_new = json.load(open('_candidates/{}/report-{}-hla.json'.format(comparison_name, inpaths['sample_id'])))
        out_old = json.load(open('outputs/{}/report-{}-hla.json'.format(comparison_name, inpaths['sample_id'])))
        out_new = sorted(out_new['hla']['alleles'])
        out_old = sorted(out_old['hla']['alleles'])
        print("")
        print("### Comparing outputs ###")
        print("\nGold standard:")
        print('\n'.join(out_old))
        print("\nNew output:")
        print('\n'.join(out_new))
        print("")
        if out_new == out_old:
            print("Outputs match.")
        else:
            print("")
            print("ERROR: candidate directory {} doesn't match the gold standard\n".format(comparison_name))
            print("")
            sys.exit(1)
