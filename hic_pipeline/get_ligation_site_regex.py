import argparse

# These ligation junctions are consistent with mega.sh
RESTRICTION_ENZYME_TO_SITE = {
    "HindIII": "AAGCTAGCTT",
    "DpnII": "GATCGATC",
    "MboI": "GATCGATC",
    "none": "XXXX",
}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--enzymes",
        required=True,
        nargs="+",
        choices=list(RESTRICTION_ENZYME_TO_SITE.keys()),
        help="The names of retriction enzymes to use to generate the regex",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="The name of the file to output the regex to",
    )
    return parser


def get_ligation_site_regex(enzymes):
    """
    If only one enzyme is inputted, then will return a single ligation site with no
    pipe (`|`), otherwise will return a regex suitable for use with `grep -E` containing
    the sites separated by pipes and surrounded in parentheses, e.g. `(AT|GC)`.

    The sites are already checked to be in the dict by the argument parser. We need to
    sort the list so the regex is the same if the list is shuffled.
    """
    sites = [RESTRICTION_ENZYME_TO_SITE[enzyme] for enzyme in enzymes]
    sites = sorted(list(set(sites)))
    if len(sites) == 1:
        return sites[0]
    return "({})".format("|".join(sites))


def main():  # pragma: nocover
    parser = get_parser()
    args = parser.parse_args()
    regex = get_ligation_site_regex(args.enzymes)
    with open(args.outfile, "w") as f:
        f.write(regex)


if __name__ == "__main__":
    main()
