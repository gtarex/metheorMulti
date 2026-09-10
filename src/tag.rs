use rust_htslib::faidx;
use rust_htslib::{
    bam,
    bam::ext::BamRecordExtensions,
    bam::record::{Aux, Cigar, Record},
    bam::Read,
};
use std::cmp::{max, min};
use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::{mpsc, Arc};
use std::thread;

use crate::bamutil;

fn reverse_complement(seq: &str, rcmapping: &HashMap<char, char>) -> String {
    let mut res: String = String::with_capacity(seq.len());
    for base in seq.chars().rev() {
        res.push(rcmapping[&base]);
    }
    res
}

fn char_at(seq: &str, i: usize) -> char {
    seq.chars().nth(i).unwrap()
}

fn is_chg_context(seq: &str) -> bool {
    matches!(seq, "CAG" | "CTG" | "CCG")
}

fn is_chh_context(seq: &str) -> bool {
    matches!(
        seq,
        "CAA" | "CAT" | "CAC" | "CTA" | "CTT" | "CTC" | "CCA" | "CCT" | "CCC"
    )
}

fn is_unknown_context(seq: &str) -> bool {
    let mut flag = false;
    for base in seq.chars() {
        if base == '-' || base == 'N' {
            flag = true;
            break;
        }
    }
    flag
}

pub fn get_header_template_from_bam(input: &str) -> bam::Header {
    let bam = bam::Reader::from_path(input).unwrap();
    bam::Header::from_template(bam.header())
}

pub fn get_tid2size_from_bam(input: &str) -> HashMap<usize, usize> {
    let header_ = get_header_template_from_bam(input);

    let mut tid2size: HashMap<usize, usize> = HashMap::new();
    for (key, records) in header_.to_hashmap() {
        // SQ = Reference sequence dictionary.
        // The order of @SQ lines defines the alignment sorting order.
        // Each @SQ entry has "SN" (Reference sequence name) and "LN" (Reference sequence length) field.
        // https://samtools.github.io/hts-specs/SAMv1.pdf
        if key != "SQ" {
            continue;
        }

        for (tid, record) in records.iter().enumerate() {
            tid2size.insert(tid, record["LN"].parse().unwrap());
        }
    }
    tid2size
}

pub fn get_rcmapping() -> HashMap<char, char> {
    // reverse-complement mapping
    let mut rcmapping: HashMap<char, char> = HashMap::new();
    rcmapping.insert('A', 'T');
    rcmapping.insert('C', 'G');
    rcmapping.insert('G', 'C');
    rcmapping.insert('T', 'A');
    rcmapping.insert('N', 'N');
    rcmapping.insert('M', 'K');
    rcmapping.insert('R', 'Y');
    rcmapping.insert('W', 'W');
    rcmapping.insert('S', 'S');
    rcmapping.insert('Y', 'R');
    rcmapping.insert('K', 'M');
    rcmapping.insert('V', 'B');
    rcmapping.insert('H', 'D');
    rcmapping.insert('D', 'H');
    rcmapping.insert('B', 'V');
    rcmapping.insert('-', '-');

    rcmapping
}

/// Lightweight payload extracted from a Record, containing only the fields
/// needed for XM tag computation. This is Send + Sync so it can cross thread boundaries.
#[derive(Clone)]
struct ReadPayload {
    tid: i32,
    reference_start: i64,
    reference_end: i64,
    seq_bytes: Vec<u8>,
    cigar_ops: Vec<Cigar>,
    is_reverse: bool,
    is_first_in_template: bool,
    is_last_in_template: bool,
}

impl ReadPayload {
    fn from_record(r: &Record) -> Self {
        Self {
            tid: r.tid(),
            reference_start: r.reference_start(),
            reference_end: r.reference_end(),
            seq_bytes: r.seq().as_bytes().to_vec(),
            cigar_ops: r.cigar().iter().copied().collect(),
            is_reverse: r.is_reverse(),
            is_first_in_template: r.is_first_in_template(),
            is_last_in_template: r.is_last_in_template(),
        }
    }
}

/// Compute XM tag string from a ReadPayload (used by worker threads).
fn compute_xm_tag_from_payload(
    payload: &ReadPayload,
    refgenome: &HashMap<usize, Vec<u8>>,
    tid2size: &HashMap<usize, usize>,
    rcmapping: &HashMap<char, char>,
    is_paired_end: bool,
) -> String {
    let tid = payload.tid;
    let start = payload.reference_start;
    let end = payload.reference_end;

    // Reverse-complement decision (same logic as need_reverse_complement for paired-end,
    // but extracting the needed flags from payload)
    let flag_reverse_complement = match is_paired_end {
        true => {
            !((!payload.is_reverse && payload.is_first_in_template)
                || (payload.is_reverse && payload.is_last_in_template))
        }
        false => payload.is_reverse,
    };

    // Extract read sequence from payload.
    let read_seq = match str::from_utf8(&payload.seq_bytes) {
        Ok(read_seq) => read_seq.to_string().to_uppercase(),
        Err(error) => panic!("Error parsing alignment record: {}", error),
    };

    // For reference sequence,
    // we should additionally consider upstream & downstream 2-bp positions,
    // to determine the cytosine context near the left & right edge of the alignment.
    let chromsize = tid2size[&(tid as usize)] as i64;
    let clipped_start = max(start - 2, 0) as usize;
    let clipped_end = min(end + 2, chromsize) as usize;

    let ref_seq_result =
        str::from_utf8(&refgenome[&(tid as usize)].as_slice()[clipped_start..clipped_end]);
    let ref_seq = match ref_seq_result {
        Ok(ref_seq) => ref_seq.to_string().to_uppercase(),
        Err(error) => panic!("Error extracting reference sequence: {}", error),
    };
    // For reads aligned at the edge of the reference genome,
    // we may not be able to extract flanking 2bp. In that case, just pad with N as much as needed.
    let padding = ["", "N", "NN"];
    let pad_nbases_start = max(2 - start, 0) as usize;
    let pad_nbases_end = max(end - chromsize + 2, 0) as usize;
    let ref_seq = format!(
        "{}{}{}",
        padding[pad_nbases_start], ref_seq, padding[pad_nbases_end]
    );

    // Separate leading and trailing soft-clips from aligned operations.
    // Soft-clipped bases do not align to the reference genome and must be marked with '.' in XM.
    // Keeping soft-clips separate prevents inserting '-' into the reference sequence,
    // which would otherwise sever true adjacent CpG (CG) contexts at alignment boundaries.
    let mut leading_clip: usize = 0;
    let mut trailing_clip: usize = 0;
    let mut aligned_cigar_ops: Vec<Cigar> = Vec::new();

    for &cigar in payload.cigar_ops.iter() {
        match cigar {
            Cigar::SoftClip(length) => {
                if aligned_cigar_ops.is_empty() {
                    leading_clip += length as usize;
                } else {
                    trailing_clip += length as usize;
                }
            }
            Cigar::HardClip(_) => {}
            other => {
                aligned_cigar_ops.push(other);
            }
        }
    }

    if aligned_cigar_ops.is_empty() {
        return std::iter::repeat('.').take(read_seq.len()).collect();
    }

    let aligned_read_seq = &read_seq[leading_clip..read_seq.len() - trailing_clip];

    let mut tmp_read_seq: Vec<char> = Vec::new();
    let mut tmp_ref_seq: Vec<char> = Vec::new();

    tmp_read_seq.push('-');
    tmp_read_seq.push('-');

    tmp_ref_seq.push(ref_seq.chars().next().unwrap());
    tmp_ref_seq.push(ref_seq.chars().nth(1).unwrap());

    let mut used_read_len: usize = 0;
    let mut used_ref_len: usize = 2;

    for cigar in aligned_cigar_ops.iter() {
        match cigar {
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                tmp_read_seq.append(
                    &mut aligned_read_seq
                        .chars()
                        .skip(used_read_len)
                        .take(*length as usize)
                        .collect(),
                );
                tmp_ref_seq.append(
                    &mut ref_seq
                        .chars()
                        .skip(used_ref_len)
                        .take(*length as usize)
                        .collect(),
                );

                used_read_len += *length as usize;
                used_ref_len += *length as usize;
            }
            Cigar::Ins(length) => {
                tmp_read_seq.append(
                    &mut aligned_read_seq
                        .chars()
                        .skip(used_read_len)
                        .take(*length as usize)
                        .collect(),
                );
                tmp_ref_seq.extend(std::iter::repeat('-').take(*length as usize));

                used_read_len += *length as usize;
            }
            Cigar::Del(length) => {
                tmp_read_seq.extend(std::iter::repeat('-').take(*length as usize));
                tmp_ref_seq.append(
                    &mut ref_seq
                        .chars()
                        .skip(used_ref_len)
                        .take(*length as usize)
                        .collect(),
                );

                used_ref_len += *length as usize;
            }
            _ => {}
        }
    }

    tmp_read_seq.push('-');
    tmp_read_seq.push('-');

    tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 2).unwrap());
    tmp_ref_seq.push(ref_seq.chars().nth(ref_seq.len() - 1).unwrap());

    let target_read_seq;
    let target_ref_seq;

    if flag_reverse_complement {
        let read = tmp_read_seq
            .iter()
            .take(tmp_read_seq.len() - 2)
            .collect::<String>();
        let reference = tmp_ref_seq
            .iter()
            .take(tmp_ref_seq.len() - 2)
            .collect::<String>();

        target_read_seq = reverse_complement(&read, rcmapping);
        target_ref_seq = reverse_complement(&reference, rcmapping);
    } else {
        target_read_seq = tmp_read_seq.iter().skip(2).collect::<String>();
        target_ref_seq = tmp_ref_seq.iter().skip(2).collect::<String>();
    }

    let mut xm_tag: Vec<char> = Vec::new();
    for idx in 0..target_read_seq.len() - 2 {
        if char_at(&target_read_seq, idx) == '-' {
            continue;
        } else if char_at(&target_read_seq, idx) == 'N' {
            xm_tag.push('.');
        } else if char_at(&target_ref_seq, idx) == 'C' {
            if (char_at(&target_read_seq, idx + 1) == '-'
                || char_at(&target_read_seq, idx + 2) == '-')
                && ((idx != target_read_seq.len() - 3) && (idx != target_read_seq.len() - 4))
            {
                let mut tmp_target_read_seq: Vec<char> = Vec::new();
                let mut tmp_target_ref_seq: Vec<char> = Vec::new();

                tmp_target_read_seq.push(char_at(&target_read_seq, idx));
                tmp_target_ref_seq.push(char_at(&target_ref_seq, idx));

                let mut flag_tmp = 0;
                let mut tmp_count = 1;

                while flag_tmp != 2 {
                    if idx + tmp_count > target_read_seq.len() - 1 {
                        break;
                    }
                    if char_at(&target_read_seq, idx + tmp_count) != '-' {
                        tmp_target_read_seq.push(char_at(&target_read_seq, idx + tmp_count));
                        tmp_target_ref_seq.push(char_at(&target_ref_seq, idx + tmp_count));
                        flag_tmp += 1;
                    }

                    tmp_count += 1;
                }

                let ref_context = tmp_target_ref_seq.iter().collect::<String>();

                if (tmp_target_ref_seq[0] == 'C') && (tmp_target_ref_seq[1] == 'G') {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('Z');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('z');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chg_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('X');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('x');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chh_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('H');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('h');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_unknown_context(&ref_context) {
                    if tmp_target_read_seq[0] == 'C' {
                        xm_tag.push('U');
                    } else if tmp_target_read_seq[0] == 'T' {
                        xm_tag.push('u');
                    } else {
                        xm_tag.push('.');
                    }
                }
            }
            // No deletion
            else {
                let ref_context = target_ref_seq.chars().skip(idx).take(3).collect::<String>();

                if (char_at(&target_ref_seq, idx) == 'C')
                    && (char_at(&target_ref_seq, idx + 1) == 'G')
                {
                    // Reference context is CG, read 'C' -> 'methylated in CG context (Z)'
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('Z');
                    }
                    // Reference context is CG, read 'T' -> 'unmethylated in CG context (z)'
                    else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('z');
                    }
                    // Reference context is CG, read 'A or G' -> Nothing.
                    else {
                        xm_tag.push('.');
                    }
                } else if is_chg_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('X');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('x');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_chh_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('H');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('h');
                    } else {
                        xm_tag.push('.');
                    }
                } else if is_unknown_context(&ref_context) {
                    if char_at(&target_read_seq, idx) == 'C' {
                        xm_tag.push('U');
                    } else if char_at(&target_read_seq, idx) == 'T' {
                        xm_tag.push('u');
                    } else {
                        xm_tag.push('.');
                    }
                }
            }
        } else {
            xm_tag.push('.');
        }
    }

    let aligned_xm: String = match flag_reverse_complement {
        true => xm_tag.iter().rev().collect::<String>(),
        false => xm_tag.iter().collect::<String>(),
    };

    let mut full_xm = String::with_capacity(read_seq.len());
    for _ in 0..leading_clip {
        full_xm.push('.');
    }
    full_xm.push_str(&aligned_xm);
    for _ in 0..trailing_clip {
        full_xm.push('.');
    }
    full_xm
}

/// Original determine_xm_tag_string (kept for backward compatibility and single-threaded path).
/// Delegates to compute_xm_tag_from_payload after extracting fields from the Record.
pub fn determine_xm_tag_string(
    r: &Record,
    refgenome: &HashMap<usize, Vec<u8>>,
    tid2size: &HashMap<usize, usize>,
    rcmapping: &HashMap<char, char>,
    is_paired_end: bool,
) -> String {
    let payload = ReadPayload::from_record(r);
    compute_xm_tag_from_payload(&payload, refgenome, tid2size, rcmapping, is_paired_end)
}

/// Single-threaded implementation. Reads BAM, processes each record sequentially,
/// writes output BAM with XM tag.
fn run_helper(input: &str, output: &str, genome: &str) {
    let mut reader = bamutil::get_reader(input);
    let is_paired_end = bamutil::is_paired_end(input);
    let header = bamutil::get_header(&reader);
    let tid2size: HashMap<usize, usize> = get_tid2size_from_bam(input);

    let rcmapping = get_rcmapping();

    // Assert if the output directory exists.
    let path = PathBuf::from(&output);
    let dir = path.parent().unwrap();

    if !dir.is_dir() {
        panic!(
            "No such directory for output alignment file: {}",
            dir.to_str().unwrap()
        )
    }
    // Prepare output writer.
    let header_tmpl = get_header_template_from_bam(input);
    let mut writer = match bam::Writer::from_path(output, &header_tmpl, bam::Format::Bam) {
        Ok(writer) => writer,
        Err(error) => panic!("Error opening alignment file to write: {}", error),
    };
    // Prepare reference genome.
    let refgenome_reader = match faidx::Reader::from_path(genome) {
        Ok(refgenome_reader) => refgenome_reader,
        Err(error) => {
            panic!("Error opening reference genome file: {}", error);
        }
    };
    println!("Parsing reference genome...");
    let mut refgenome: HashMap<usize, Vec<u8>> = HashMap::new();
    for (tid, _size) in tid2size.iter() {
        let ref_array = refgenome_reader
            .fetch_seq(bamutil::tid2chrom(*tid as i32, &header), 0, tid2size[tid])
            .expect("Error fetching reference genome sequence.");

        refgenome.insert(*tid, ref_array);
    }
    println!("Done!");

    // Main loop
    // Iterate aligned reads and determine xm tag string.
    for mut r in reader.records().map(|r| r.unwrap()) {
        // Determine XM tag string by comparing read sequence and reference sequence.
        let xm_tag_string =
            determine_xm_tag_string(&r, &refgenome, &tid2size, &rcmapping, is_paired_end);
        // Attach XM tag to the record.
        let add_result = r.push_aux("XM".as_bytes(), Aux::String(&xm_tag_string));
        match add_result {
            Ok(_) => (),
            Err(e) => panic!("Error adding XM tag to alignment record. {}", e),
        }
        // Write record to output.
        writer.write(&r).expect("Error writing to output file.");
    }
}

/// Multi-threaded implementation using channel-based worker pool.
/// Produces exactly the same results as `run_helper`.
fn run_helper_mt(input: &str, output: &str, genome: &str, threads: usize) {
    let mut reader = bamutil::get_reader(input);
    let is_paired_end = bamutil::is_paired_end(input);
    let header = bamutil::get_header(&reader);
    let tid2size: HashMap<usize, usize> = get_tid2size_from_bam(input);

    let rcmapping = get_rcmapping();

    // Assert if the output directory exists.
    let path = PathBuf::from(&output);
    let dir = path.parent().unwrap();

    if !dir.is_dir() {
        panic!(
            "No such directory for output alignment file: {}",
            dir.to_str().unwrap()
        )
    }
    // Prepare output writer.
    let header_tmpl = get_header_template_from_bam(input);
    let mut writer = match bam::Writer::from_path(output, &header_tmpl, bam::Format::Bam) {
        Ok(writer) => writer,
        Err(error) => panic!("Error opening alignment file to write: {}", error),
    };
    // Prepare reference genome.
    let refgenome_reader = match faidx::Reader::from_path(genome) {
        Ok(refgenome_reader) => refgenome_reader,
        Err(error) => {
            panic!("Error opening reference genome file: {}", error);
        }
    };
    println!("Parsing reference genome...");
    let mut refgenome: HashMap<usize, Vec<u8>> = HashMap::new();
    for (tid, _size) in tid2size.iter() {
        let ref_array = refgenome_reader
            .fetch_seq(bamutil::tid2chrom(*tid as i32, &header), 0, tid2size[tid])
            .expect("Error fetching reference genome sequence.");

        refgenome.insert(*tid, ref_array);
    }
    println!("Done!");

    // Wrap read-only data in Arc for sharing across threads.
    let refgenome = Arc::new(refgenome);
    let tid2size = Arc::new(tid2size);
    let rcmapping = Arc::new(rcmapping);

    // Collect all records and payloads first.
    // This is necessary because Record is not Clone/Send, so we extract payloads
    // for workers and keep records for writing in correct order.
    println!("Reading alignment records...");
    let mut records: Vec<Record> = Vec::new();
    let mut payloads: Vec<ReadPayload> = Vec::new();
    for r in reader.records().map(|r| r.unwrap()) {
        payloads.push(ReadPayload::from_record(&r));
        records.push(r);
    }
    println!("Loaded {} records.", records.len());

    let total_records = records.len();

    // Create channels: workers receive Option<(usize, ReadPayload)>, send back (usize, String)
    let mut txs: Vec<mpsc::SyncSender<Option<(usize, ReadPayload)>>> = Vec::with_capacity(threads);

    let workers: Vec<thread::JoinHandle<Vec<(usize, String)>>> = (0..threads)
        .map(|_| {
            let (tx, rx) = mpsc::sync_channel::<Option<(usize, ReadPayload)>>(4);
            txs.push(tx);
            let refgenome = Arc::clone(&refgenome);
            let tid2size = Arc::clone(&tid2size);
            let rcmapping = Arc::clone(&rcmapping);
            thread::spawn(move || {
                let mut local_results: Vec<(usize, String)> = Vec::new();
                while let Some((idx, payload)) = rx.recv().unwrap() {
                    let xm_tag = compute_xm_tag_from_payload(
                        &payload,
                        &refgenome,
                        &tid2size,
                        &rcmapping,
                        is_paired_end,
                    );
                    local_results.push((idx, xm_tag));
                }
                local_results
            })
        })
        .collect();

    // Dispatch work: round-robin distribution of payload indices.
    let mut worker_idx: usize = 0;
    for (i, payload) in payloads.into_iter().enumerate() {
        txs[worker_idx]
            .send(Some((i, payload)))
            .expect("Error sending to worker thread");
        worker_idx = (worker_idx + 1) % threads;
    }

    // Signal all workers to finish.
    for tx in &txs {
        tx.send(None).expect("Error sending termination signal");
    }

    // Collect results from all workers and place into ordered Vec.
    let mut xm_tags: Vec<Option<String>> = vec![None; total_records];
    for worker in workers {
        let local_results = worker.join().expect("Worker thread panicked");
        for (idx, xm_tag) in local_results {
            xm_tags[idx] = Some(xm_tag);
        }
    }

    // Apply XM tags to records in original order and write.
    println!("Writing tagged records...");
    for (i, mut r) in records.into_iter().enumerate() {
        let xm_tag_string = xm_tags[i]
            .as_ref()
            .expect("Missing XM tag for record index");
        let add_result = r.push_aux("XM".as_bytes(), Aux::String(xm_tag_string));
        match add_result {
            Ok(_) => (),
            Err(e) => panic!("Error adding XM tag to alignment record. {}", e),
        }
        writer.write(&r).expect("Error writing to output file.");
    }
    println!("Done writing!");
}

pub fn run(input: &str, output: &str, genome: &str, threads: usize) {
    if threads <= 1 {
        run_helper(input, output, genome);
    } else {
        run_helper_mt(input, output, genome, threads);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let rcmapping = get_rcmapping();
        assert_eq!("ATGC", reverse_complement("GCAT", &rcmapping));
    }
    #[test]
    #[should_panic]
    fn error_when_output_directory_is_not_found() {
        run(
            "tests/test1.bam",
            "tests/no_such_directory/out.bam",
            "tests/tinyref.fa",
            1,
        )
    }
    #[test]
    #[should_panic]
    fn error_when_reference_genome_is_not_found() {
        run(
            "tests/test1.bam",
            "tests/out.tagged.bam",
            "tests/there_is_no_such.fa",
            1,
        )
    }

    #[test]
    fn test_mt_matches_st() {
        // Verify that multi-threaded produces the same tagged BAM as single-threaded.
        use std::process::Command;

        let input = "tests/test1.bam";
        let genome = "tests/tinyref.fa";

        // Run single-threaded
        run(input, "tests/out.tagged.bam", genome, 1);

        // Run multi-threaded
        run(input, "tests/out.tagged.mt.bam", genome, 4);

        // Compare BAM files using samtools view
        let st_output = Command::new("samtools")
            .args(["view", "tests/out.tagged.bam"])
            .output()
            .expect("Failed to run samtools view on ST output");
        let mt_output = Command::new("samtools")
            .args(["view", "tests/out.tagged.mt.bam"])
            .output()
            .expect("Failed to run samtools view on MT output");

        assert_eq!(
            st_output.stdout, mt_output.stdout,
            "Single-threaded and multi-threaded BAM outputs differ"
        );
    }
}