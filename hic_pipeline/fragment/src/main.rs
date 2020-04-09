// Intended to be a drop-in replacement for
// https://github.com/aidenlab/juicer/blob/encode/CPU/common/fragment.pl

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use structopt::StructOpt;

// Adapted from Rust CLI tutorial
#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    infile: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    outfile: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    site_file: std::path::PathBuf,
}

fn bsearch(target_value: i64, positions: &[i64]) -> i64 {
    // Assumes that the positions are already sorted
    let mut lower_index: i64 = 0;
    let mut upper_index = positions.len() as i64 - 1;
    let mut current_index;
    while lower_index <= upper_index {
        current_index = (lower_index + upper_index) / 2;
        let current_value = positions[current_index as usize];
        match current_value.cmp(&target_value) {
            Ordering::Greater => upper_index = current_index - 1,
            Ordering::Less => lower_index = current_index + 1,
            Ordering::Equal => return current_index + 1,
        }
    }
    lower_index
}

fn process_site_file_line(
    sites_by_chr: &mut HashMap<std::string::String, std::vec::Vec<i64>>,
    line: &io::Result<std::string::String>,
) {
    // https://users.rust-lang.org/t/borrowed-value-does-not-live-long-enough/7225/2
    let mut split_line = line.as_ref().unwrap().split_whitespace();
    let chr = split_line.next().unwrap();
    let values: Vec<i64> = split_line.map(|x| x.parse::<i64>().unwrap()).collect();
    // I have no idea why 14 is special cased, doesn't apply to human chrom sizes files
    if chr == "14" {
        sites_by_chr.insert(
            format!("{}{}", chr, "m"),
            // Can't use values twice, so need to copy
            values.clone(),
        );
        sites_by_chr.insert(format!("{}{}", chr, "p"), values);
    } else {
        sites_by_chr.insert(chr.to_string(), values);
    }
}

fn process_line(
    line: &io::Result<std::string::String>,
    sites_by_chr: &HashMap<std::string::String, std::vec::Vec<i64>>,
) -> std::string::String {
    let split_line = line
        .as_ref()
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>();

    let index1 = bsearch(
        split_line[2].parse::<i64>().unwrap(),
        &sites_by_chr[split_line[1]],
    );
    let index2 = bsearch(
        split_line[5].parse::<i64>().unwrap(),
        &sites_by_chr[split_line[4]],
    );
    format!(
        "{} {} {} {} {} \n",
        &split_line[0..3].join(" "),
        index1,
        &split_line[3..6].join(" "),
        index2,
        &split_line[6..].join(" "),
    )
}

fn main() -> io::Result<()> {
    let args = Cli::from_args();
    let site_file = File::open(&args.site_file)?;
    let site_file_reader = BufReader::new(site_file);
    let mut sites_by_chr = HashMap::new();

    for line in site_file_reader.lines() {
        process_site_file_line(&mut sites_by_chr, &line)
    }

    let infile = File::open(&args.infile)?;
    let infile_reader = BufReader::new(infile);
    let outfile = File::create(&args.outfile)?;
    let mut outfile_writer = BufWriter::new(outfile);

    for line in infile_reader.lines() {
        outfile_writer.write_all(process_line(&line, &sites_by_chr).as_bytes())?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bsearch_target_not_in_positions() {
        let target_value = 3;
        let positions = [1, 4, 5];
        let result = bsearch(target_value, &positions);
        assert_eq!(result, 1)
    }

    #[test]
    fn test_bsearch_target_in_positions() {
        let target_value = 4;
        let positions = [1, 4, 5];
        let result = bsearch(target_value, &positions);
        assert_eq!(result, 2)
    }

    #[test]
    fn test_process_site_file_line() {
        let mut sites_by_chr = HashMap::new();
        let line = Ok("chrIII 968 2310 2541".to_string());
        process_site_file_line(&mut sites_by_chr, &line);
        assert_eq!(sites_by_chr["chrIII"], [968, 2310, 2541])
    }

    #[test]
    fn test_process_site_file_line_chr14_duplicates() {
        let mut sites_by_chr = HashMap::new();
        let line = Ok("14 1 2".to_string());
        process_site_file_line(&mut sites_by_chr, &line);
        assert_eq!(sites_by_chr["14m"], [1, 2]);
        assert_eq!(sites_by_chr["14p"], [1, 2])
    }

    #[test]
    fn test_process_line() {
        let mut sites_by_chr = HashMap::new();
        sites_by_chr.insert("chrIII".to_string(), vec![1, 3, 5]);
        let line = Ok("0       chrIII  2        0       chrIII  3        15      17H30M1D3M      T       8       30M1D3M T       SRR    SRR".to_string());
        let result = process_line(&line, &sites_by_chr);
        assert_eq!(
            result,
            "0 chrIII 2 1 0 chrIII 3 2 15 17H30M1D3M T 8 30M1D3M T SRR SRR \n".to_string()
        );
    }
}
