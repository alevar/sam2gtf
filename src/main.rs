extern crate clap;
extern crate rust_htslib;

use std::fs;
use std::io::Write;

use clap::{Command,Arg, ArgAction};
use rust_htslib::bam;
use rust_htslib::bam::Read;

// the main logic of the app. perform conversion from SAM/BAM to GTF
fn convert(input_fname: String, output_fname: String, keep_multi: bool) {
    let mut reader = bam::Reader::from_path(input_fname).unwrap();
    let header =  reader.header().to_owned();

    // prepare the output file
    let mut gtf_writer = fs::File::create(output_fname).expect("Could not create output file");
    
    let mut read_count: u64 = 0;
    for r in reader.records() {
        read_count += 1;

        let record = r.unwrap();

        // if unmapped = skip
        if record.is_unmapped() {
            continue;
        }

        // check if multimapper and if keep_multi = false, skip
        if !keep_multi && record.is_secondary() {
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
                    if (exon_end - exon_start) > 0 {
                        exons.push((exon_start+1, exon_end));
                    }
                    exon_start = exon_end + len;
                    exon_end = exon_start;
                },
                bam::record::Cigar::SoftClip(len) => {
                    // do nothing
                },
                bam::record::Cigar::HardClip(len) => {
                    // do nothing
                },
                _ => {
                    println!("Unknown CIGAR operation {}", qname);
                }
            }
        }
        // push the last exon
        exons.push((exon_start+1, exon_end));
        
        // write transcript line. qname as transcript_id
        let sequence_id = record.tid();
        let seq_name_bytes = header.tid2name(sequence_id as u32);
        let seq_name = std::str::from_utf8(seq_name_bytes).unwrap();
        let tx_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttranscript_id \"{}\"; read_name \"{}\";\n", seq_name, "SAM2GTF", "transcript", exons.first().unwrap().0, exons.last().unwrap().1, ".", strand, ".", read_count, qname);
        gtf_writer.write_all(tx_line.as_bytes()).expect("Could not write to output file");

        // write exon lines. qname as transcript_id
        for exon in exons {
            let exon_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttranscript_id \"{}\"; read_name \"{}\";\n", seq_name, "SAM2GTF", "exon", exon.0, exon.1, ".", strand, ".", read_count, qname);
            gtf_writer.write_all(exon_line.as_bytes()).expect("Could not write to output file");
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
        .arg(
            Arg::new("keep_multi")
            .short('k')
            .long("keep_multi")
            .required(false)
            .action(ArgAction::SetTrue)
            .help("If enable, will convert all multimappers. By default this is disabled and only primary alignments are converted."),

        )
        .after_help("--help or -h")
        .get_matches();

    // Extract and use the value of the "input" argument
    let input_fname: &String = matches.get_one("input").unwrap();
    let output_fname: &String = matches.get_one("output").unwrap();
    let keep_multi: bool = matches.get_flag("keep_multi");

    convert(input_fname.to_string(), output_fname.to_string(), keep_multi);
}
