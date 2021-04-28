import argparse
import json
import os


def main():
    parser = get_parser()
    args = parser.parse_args()
    with open(args.filepath) as f:
        data = f.readlines()
    parsed_data = parse_to_dict(data)
    jsonified = jsonify_stats(parsed_data)
    filename = os.path.splitext(args.filepath)[0] + ".json"
    with open(filename, "w") as f:
        f.write(json.dumps(jsonified, indent=4, sort_keys=True))


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filepath", type=str, help="path to the stats file to convert to json"
    )
    return parser


def parse_to_dict(data):
    """
    Convert a text file of key: value pairs to a dictionary, preprocessing the
    field names such that "Key Name" becomes "key_name" in the resulting
    dictionary
    """
    output = {}
    for line in data:
        key, value = line.split(":")
        if not value.strip():
            continue
        output[clean_key(key)] = value.strip()
    return output


def clean_key(key):
    """
    Rather complex key parsing to make the keys look more normal.
    """
    key = "_".join(key.strip().split(" ")).lower()
    key = key.replace("+", "_and_")
    key = key.replace("hi-c", "hic")
    if "bp" in key or "kb" in key:
        key = key.replace("-", "_to_")
    else:
        key = key.replace("-", "_")
    key = key.replace("(", "")
    key = key.replace(")", "")
    key = key.replace("'", "_prime_")
    key = key.replace("%", "_percent_")
    key = key.replace(">", "_greater_than_")
    key = key.replace("<", "_less_than_")
    key = key.lstrip("_")
    key = key.replace("__", "_")
    return key


def jsonify_stats(parsed_data):
    """
    Does more fine grained parsing on the qc values, including converting some flat
    strings into structured objects.
    """
    value_sep = "-"
    output = {}
    for k, v in parsed_data.items():
        if "%" in v:
            if value_sep not in v:
                split_value = v.split("(")
                values = {
                    "total": split_value[0].strip(),
                    "percent": split_value[1].rstrip("%)"),
                }
            else:
                split_value = v.split(value_sep)
                split_value = [i.strip() for i in split_value]
                if "bias" in k:
                    fields = ["long_range", "short_range"]
                elif "pair" in k:
                    fields = ["left", "inner", "outer", "right"]
                values = dict(zip(fields, split_value))
            output[k] = values
        else:
            output[k] = v
    return output


if __name__ == "__main__":
    main()
