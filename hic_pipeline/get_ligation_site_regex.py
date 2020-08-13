import argparse

# These ligation junctions are consistent with mega.sh
RESTRICTION_ENZYME_TO_SITE = {
    "HindIII": "AAGCTAGCTT",
    "DpnII": "GATCGATC",
    "MboI": "GATCGATC",
}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--enzymes",
        required=True,
        nargs="+",
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
    sites = []
    for enzyme in enzymes:
        site = RESTRICTION_ENZYME_TO_SITE.get(enzyme)
        if site is None:
            raise ValueError(
                "Restriction enzyme not recognized, valid options are {}".format(
                    list(RESTRICTION_ENZYME_TO_SITE.keys())
                )
            )
        sites.append(site)
    sites = list(set(sites))
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
