import argparse


def main():
    parser = get_parser()
    args = parser.parse_args()
    with open(args.chrom_sizes) as chrom_sizes, open(
        args.output_file, "w"
    ) as output_file:
        filter_chrom_sizes(chrom_sizes=chrom_sizes, output_file=output_file)


def filter_chrom_sizes(chrom_sizes, output_file):
    for line in chrom_sizes:
        chrom_name = line.split()[0]
        if "_" in chrom_name or chrom_name == "chrEBV":
            continue
        output_file.write(line)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("chrom_sizes")
    parser.add_argument("output_file")
    return parser


if __name__ == "__main__":
    main()
