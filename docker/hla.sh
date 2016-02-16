#!/bin/bash
set -e

DOCKER=$1
BAM_FILE=$2
SAMPLE_ID=$3

KEY=`perl -ne 'print $1 if m/aws_access_key_id = (\\S+)/' ~/.aws/credentials`
SECRET=`perl -ne 'print $1 if m/aws_secret_access_key = (\\S+)/' ~/.aws/credentials`

docker run --rm -u `id -u`:`id -g` -v `pwd`:/work -w /work \
	-e AWS_ACCESS_KEY_ID=$KEY \
	-e AWS_SECRET_ACCESS_KEY=$SECRET \
    $DOCKER typer.sh $BAM_FILE $SAMPLE_ID
