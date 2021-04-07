#!/bin/bash

set -euo pipefail

if [ -z "$(command -v caper)" ]; then
    echo "caper not installed, install with 'pip install caper'"
    exit 1
fi

if [ $# -lt 2 ]; then
    echo "Usage: ./caper_run.sh [WDL] [INPUT_JSON]"
    exit 1
fi

if [ -z "${HIC_DOCKER_IMAGE_TAG}" ]; then
    echo "Must specify HIC_DOCKER_IMAGE_TAG via environment variable."
    exit 1
fi

WDL=$1
INPUT=$2
WORKFLOW_OPTIONS="tests/pytest_workflow_options.json"
WORKFLOW_OPTIONS_FLAG="--no-relative-output-paths"

if [ $# -gt 2 ]; then
    if [ "$3" != "${WORKFLOW_OPTIONS_FLAG}" ]; then
        echo "Third argument must be ${WORKFLOW_OPTIONS_FLAG}"
        exit 1
    fi
    WORKFLOW_OPTIONS="tests/pytest_workflow_no_relative_output_paths.json"
fi

echo "Running caper with WDL ${WDL}, input ${INPUT}, workflow options ${WORKFLOW_OPTIONS}, and image ${HIC_DOCKER_IMAGE_TAG}"

caper run "${WDL}" -i "${INPUT}" --docker "${HIC_DOCKER_IMAGE_TAG}" -o "./${WORKFLOW_OPTIONS}"

if [[ -f "cromwell.out" ]]; then
    cat cromwell.out
fi
