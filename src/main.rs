extern crate clap;
extern crate rust_htslib;

use clap::{Command,Arg};
use rust_htslib::bam;
use rust_htslib::bam::Read;

// the main logic of the app. perform conversion from SAM/BAM to GTF
fn convert(input_fname: String, output_fname: String) {
    let mut reader = bam::Reader::from_path(input_fname).unwrap();
    
    for r in reader.records() {
        let record = r.unwrap();

        // if unmapped = skip
        if record.is_unmapped() {
            continue;
        }

        let qname = String::from_utf8_lossy(record.qname());
        let pos = record.pos();
        let cigar = record.cigar();

        // determine strand of the read
        let strand = match record.is_reverse() {
            true => "-",
            false => "+"
        };

        // extract contiguous blocks of exons from the CIGAR string
        let mut exons: Vec<(u32, u32)> = Vec::new();
        let mut exon_start: u32 = pos as u32;
        let mut exon_end: u32 = pos as u32;
        let mut started = false; // set to false until the first modifier operation is found (M)
        for c in &cigar {
            match c {
                bam::record::Cigar::Match(len) => {
                    exon_end += len;
                },
                bam::record::Cigar::Equal(len) => {
                    exon_end += len;
                },
                bam::record::Cigar::Diff(len) => {
                    exon_end += len;
                },
                bam::record::Cigar::Ins(len) => {
                    // do nothing
                },
                bam::record::Cigar::Del(len) => {
                    exon_end += len;
                },
                bam::record::Cigar::RefSkip(len) => {
                    exons.push((exon_start, exon_end));
                    exon_start = exon_end + len;
                    exon_end = exon_start;
                },
                bam::record::Cigar::SoftClip(len) => {
                    // if started = true, then this is the end of the read
                    if started {
                        exons.push((exon_start, exon_end));
                    }
                    exon_start += len;
                    started = true;
                },
                bam::record::Cigar::HardClip(len) => {
                    // do nothing
                },
                _ => {
                    println!("Unknown CIGAR operation {}", qname);
                }
            }
        }
        for exon in &exons {
            println!("{}-{}", exon.0, exon.1);
        }
    }
}

fn main() {
    
    let matches = Command::new("SAM 2 GTF")
        .version("0.1.0")
        .author("Ales Varabyou")
        .about("Convert SAM/BAM to GTF directly. useful for investigating long-read alignments such as PacBio ISO-Seq. or ONT.")
        .arg(
            Arg::new("input")
            .short('i')
            .long("input")
        )
        .arg(
            Arg::new("output")
            .short('o')
            .long("output")
        )
        .after_help("--help or -h")
        .get_matches();

    // Extract and use the value of the "input" argument
    let input_fname: &String = matches.get_one("input").unwrap();
    let output_fname: &String = matches.get_one("output").unwrap();

    convert(input_fname.to_string(), output_fname.to_string());
}
