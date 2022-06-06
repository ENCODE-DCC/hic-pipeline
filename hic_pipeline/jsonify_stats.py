import argparse
import json


def main():
    parser = get_parser()
    args = parser.parse_args()
    data = load_data(args.infile)
    processed = process_data(data)
    write_json_to_file(processed, args.outfile)


def load_data(infile):
    with open(infile, encoding="utf-8") as f:
        return f.read()


def write_json_to_file(data, outfile):
    with open(outfile, "w") as f:
        f.write(json.dumps(data, indent=4, sort_keys=True))


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="path to the stats file to convert to json")
    parser.add_argument("outfile", help="Name of file to output stats to")
    return parser


def process_data(data):
    split = data.splitlines()
    parsed_data = parse_to_dict(split)
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
    value, first_percentage_name="pct_sequenced", second_percentage_name="pct_unique"
):
    """
    For unique and duplicates, they are percentage of alignable / percentage of
    sequenced reads. For the rest, they are percentage of sequenced reads / percentage of
    unique reads.

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


def parse_run_type(value):
    return value.lower().replace(" ", "-") + "ed"


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
    key = key.replace(";", "")
    key = key.replace("*", "")
    key = key.replace("...", "_")
    key = key.replace(".", "_")
    key = key.replace("'", "_prime_")
    key = key.replace("%", "_percent_")
    key = key.replace(">", "_greater_than_")
    key = key.replace("<", "_less_than_")
    key = key.lstrip("_")
    key = key.replace("__", "_")
    key = key.replace("l_i_o_r", "lior")
    key = key.replace("average", "avg")
    if key == "library_complexity_estimate_1_and_2_above":
        key = key.replace("_above", "_alignments")
    if key == "unmapped":
        key += "_reads"
    if key == "read_type":
        key = "run_type"
    return key


def update_with_total_and_percentages(output, parsed, parent_key):
    for key, value in parsed.items():
        if key == "total":
            output[parent_key] = value
        else:
            output["{}_{}".format(key, parent_key)] = value


def jsonify_stats(parsed_data):
    """
    Does more fine grained parsing on the qc values, including converting some flat
    strings into structured objects.
    """
    output = {}
    for k, v in parsed_data.items():
        if "%" in v and "-" in v:
            if "bias" in k:
                fields = ("pct_3_prime_bias", "pct_5_prime_bias")
            elif "pair_type" in k:
                fields = ("pct_left", "pct_inner", "pct_outer", "pct_right")
            try:
                parsed = parse_to_multiple_percentages(v, fields)
            except ValueError:
                continue
            if "bias" in k:
                update_with_total_and_percentages(output, parsed, "long_range")
            elif "pair_type" in k:
                update_with_total_and_percentages(output, parsed, "pair_type")
        elif k in ("unique_reads", "duplicates"):
            parsed = parse_to_total_and_two_percentages(
                v,
                first_percentage_name="pct_alignable",
                second_percentage_name="pct_sequenced",
            )
            update_with_total_and_percentages(output, parsed, k)
        elif k in (
            "1_alignment_unique",
            "1_alignment_duplicates",
            "2_alignment_unique",
            "2_alignment_duplicates",
        ):
            parsed = parse_to_total_and_two_percentages(
                v,
                first_percentage_name="pct",
                second_percentage_name="pct_sequenced",
            )
            update_with_total_and_percentages(output, parsed, k)
        elif v.count("%") == 1:
            parsed = parse_to_total_and_percentage(v)
            update_with_total_and_percentages(output, parsed, k)
        elif v.count("%") == 2:
            parsed = parse_to_total_and_two_percentages(v)
            update_with_total_and_percentages(output, parsed, k)
        elif k == "run_type":
            output[k] = parse_run_type(v)
        else:
            try:
                output[k] = parse_to_int_or_float(v)
            except ValueError:
                continue
    try:
        pct_sequenced_total_duplicates = output["pct_sequenced_total_duplicates"]
        pct_unique_total_duplicates = output["pct_unique_total_duplicates"]
        output["pct_sequenced_total_duplicates"] = pct_unique_total_duplicates
        output["pct_unique_total_duplicates"] = pct_sequenced_total_duplicates
    except KeyError:
        pass

    try:
        pct_sequenced_total_unique = output["pct_sequenced_total_unique"]
        pct_unique_total_unique = output["pct_unique_total_unique"]
        output["pct_sequenced_total_unique"] = pct_unique_total_unique
        output["pct_unique_total_unique"] = pct_sequenced_total_unique
    except KeyError:
        pass

    return output


if __name__ == "__main__":
    main()
