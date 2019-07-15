#!/bin/bash
set -e # exit on error

if [ $# -lt 2 ]; then
  echo "Usage: ./test.sh [WDL] [INPUT_JSON] [DOCKER_IMAGE](optional)"
  echo "Make sure to have caper installed (pip install caper)."
  exit 1
fi

WDL=$1
INPUT=$2
if [ $# -gt 2 ]; then
  DOCKER_IMAGE=$3
else
  DOCKER_IMAGE=quay.io/encode-dcc/hic-pipeline:template
fi
PREFIX=$(basename ${WDL} .wdl)
METADATA=${PREFIX}.metadata.json # metadata

caper run ${WDL} --docker ${DOCKER_IMAGE} -i ${INPUT} -m ${METADATA}
