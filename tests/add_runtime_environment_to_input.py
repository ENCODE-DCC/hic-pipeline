import argparse
import json


def main():
    parser = _get_parser()
    args = parser.parse_args()
    with open(args.input_json_path) as f:
        input_json = json.load(f)
    workflow_name = next(iter(input_json.keys())).split(".")[0]
    input_json["{}.runtime_environment".format(workflow_name)] = {
        "docker": args.docker_image,
        "singularity": "",
    }
    with open(args.input_json_path, "w") as f:
        json.dump(input_json, f, indent=2)


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json_path")
    parser.add_argument("docker_image")
    return parser


if __name__ == "__main__":
    main()
