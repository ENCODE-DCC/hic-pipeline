#!/bin/bash

set -euo pipefail

if [ -z "$(command -v caper)" ]; then
    echo "caper not installed, install with `pip install caper`"
    exit 1
fi

if [ $# -lt 2 ]; then
    echo "Usage: ./caper_run.sh [WDL] [INPUT_JSON]"
    exit 1
fi

if [ -z "${DOCKER_IMAGE}" ]; then
    echo "Must specify HIC_DOCKER_IMAGE via environment variable."
    exit 1
fi

WDL=$1
INPUT=$2

echo "Running caper with WDL ${WDL}, input ${INPUT}, and image ${HIC_DOCKER_IMAGE_TAG}"

caper run "${WDL}" -i "${INPUT}" --docker "${HIC_DOCKER_IMAGE}" -o ./tests/cromwell_options.json
