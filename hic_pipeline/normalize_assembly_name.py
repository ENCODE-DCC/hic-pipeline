import argparse

VALID_ASSEMBLIES = (
    "hg18",
    "hg19",
    "hg38",
    "dMel",
    "mm9",
    "mm10",
    "anasPlat1",
    "bTaurus3",
    "canFam3",
    "equCab2",
    "galGal4",
    "Pf3D7",
    "sacCer3",
    "sCerS288c",
    "susScr3",
    "TAIR10",
)


def main():
    parser = get_parser()
    args = parser.parse_args()
    normalized_name, assembly_is_supported = normalize_assembly_name(args.assembly_name)
    write_string_to_file(normalized_name, args.normalized_name_outfile)
    write_string_to_file(
        get_wdl_boolean_string(assembly_is_supported),
        args.assembly_is_supported_outfile,
    )


def normalize_assembly_name(assembly_name):
    """
    Returns a tuple of the possibly normalized assembly name and a boolean indicating
    whether or not the assembly is supported.
    """
    assembly_name = assembly_name.lower()
    if "grch" in assembly_name:
        assembly_name = normalize_grch_name(assembly_name)
    for canonical_name in VALID_ASSEMBLIES:
        if canonical_name.lower() == assembly_name:
            return canonical_name, True
    return assembly_name, False


def normalize_grch_name(assembly_name):
    """
    Convert `GRCh` names to `hg` ones, GRCh38 = hg38, GRCh37 = hg19, and GRCh36 = hg18
    """
    assembly_version = int(assembly_name.lower().replace("grch", ""))
    if assembly_version < 38:
        hg_assembly_version = assembly_version - 18
    else:
        hg_assembly_version = assembly_version
    return "hg" + str(hg_assembly_version)


def get_wdl_boolean_string(boolean):
    """
    WDL expects `true` or `false` strings for `read_boolean`, Python `str` doesn't work
    """
    return str(boolean).lower()


def write_string_to_file(data, filename):
    with open(filename, "w") as f:
        f.write(data)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_name", help="Assembly name")
    parser.add_argument(
        "normalized_name_outfile", help="Name of file to write normalized name"
    )
    parser.add_argument(
        "assembly_is_supported_outfile",
        help=(
            "Name for file to write boolean indicating if the assembly is supported by "
            "Juicer"
        ),
    )
    return parser


if __name__ == "__main__":
    main()
