import argparse
import json


def main():
    parser = get_parser()
    args = parser.parse_args()
    with open(args.filepath) as f:
        data = f.readlines()
    processed = process_data(data)
    with open(args.outfile, "w") as f:
        f.write(json.dumps(processed, indent=4, sort_keys=True))


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "infile", help="path to the stats file to convert to json"
    )
    parser.add_argument(
        "outfile", help="Name of file to output stats to"
    )
    return parser


def process_data(data):
    parsed_data = parse_to_dict(data)
    jsonified = jsonify_stats(parsed_data)
    return jsonified


def parse_to_dict(data):
    """
    Convert a text file of key: value pairs to a dictionary, preprocessing the
    field names such that "Key Name" becomes "key_name" in the resulting
    dictionary
    """
    output = {}
    for line in data:
        key, value = line.split(":")
        value = value.strip()
        if not value:
            continue
        output[clean_key(key)] = value
    return output


def parse_to_total_and_percentage(value):
    total, percent = value.split("(")
    values = {
        "total": parse_int_with_commas(total.strip()),
        "pct": float(percent.rstrip("%)")),
    }
    return values


def parse_to_total_and_two_percentages(
    value,
    first_percentage_name="pct_of_unique_reads",
    second_percentage_name="pct_of_sequenced_reads",
):
    """
    For unique and duplicates, they are percentage of alignable / percentage of
    sequenced reads. For the rest, they are percentage of unique reads / percentage of
    sequenced reads.

    The separator of the two percentages can be either a slash and a comma.
    """
    total, percentages = value.split("(")
    try:
        first_percentage, second_percentage = percentages.split("/")
    except ValueError:
        first_percentage, second_percentage = percentages.split(",")
    values = {
        "total": parse_int_with_commas(total.strip()),
        first_percentage_name: float(first_percentage.rstrip("% ")),
        second_percentage_name: float(second_percentage.rstrip("%)")),
    }
    return values


def parse_to_multiple_percentages(value, fields):
    """
    Handles QC with multiple percentage values like 25% - 23% - 27% - 25%
    """
    split_value = [float(i.strip(" %")) for i in value.split("-")]
    values = dict(zip(fields, split_value))
    return values


def parse_to_int_or_float(value):
    """
    Tries to parse `value` as float or int, if not possible will raise a ValueError
    """
    try:
        return parse_int_with_commas(value)
    except ValueError:
        return float(value)


def parse_int_with_commas(value):
    """
    Python `int` builtin does not handle strings with commas and this is more readable
    than using `atoi`
    """
    return int(value.replace(",", ""))


def clean_key(key):
    """
    Rather complex key parsing to make the keys look more normal. Converts special
    characters to English.
    """
    key = "_".join(key.strip().split(" ")).lower()
    key = key.replace("+", "_and_")
    key = key.replace("hi-c", "hic")
    if "bp" in key or "kb" in key:
        key = key.replace("-", "_to_")
        if "long" not in key:
            key = "short_range_" + key
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
    key = key.replace("l_i_o_r", "lior")
    return key


def jsonify_stats(parsed_data):
    """
    Does more fine grained parsing on the qc values, including converting some flat
    strings into structured objects.
    """
    output = {}
    for k, v in parsed_data.items():
        if "%" in v and "-" in v:
            if "bias" in k:
                fields = ("pct_long_range", "pct_short_range")
            elif "pair_type" in k:
                fields = ("pct_left", "pct_inner", "pct_outer", "pct_right")
            try:
                values = parse_to_multiple_percentages(v, fields)
            except ValueError:
                continue
            output[k] = values
        elif k in ("unique_reads", "duplicates"):
            output[k] = parse_to_total_and_two_percentages(
                v,
                first_percentage_name="pct_of_alignable_reads",
                second_percentage_name="pct_of_sequenced_reads",
            )
        elif v.count("%") == 1:
            output[k] = parse_to_total_and_percentage(v)
        elif v.count("%") == 2:
            output[k] = parse_to_total_and_two_percentages(v)
        else:
            try:
                output[k] = parse_to_int_or_float(v)
            except ValueError:
                continue
    return output


if __name__ == "__main__":
    main()
