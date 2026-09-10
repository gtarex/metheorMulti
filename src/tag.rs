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

const BATCH_SIZE: usize = 8192;

// Decompression, calculation, compression, main/I/O threads.
fn thread_allocation(threads: usize) -> (usize, usize, usize, usize) {
    assert!((1..=100).contains(&threads), "threads must be from 1 to 100");
    match threads {
        1..=2 => (0, 0, 0, 1),
        3..=7 => (0, threads - 2, 0, 2),
        _ => {
            let remaining = threads - 4;
            let decompress = (remaining / 10).max(1).min(8);
            let compress = (remaining - decompress) / 2;
            (decompress, remaining - decompress - compress, compress, 4)
        }
    }
}

/// One bounded input/output channel per worker; consume in dispatch order.
/// Record ownership moves between stages, with no whole-file accumulation.
fn stream_batches<I, F, W>(
    records: I,
    workers: usize,
    batch_size: usize,
    transform: F,
    mut write: W,
) -> Result<(), String>
where
    I: Iterator<Item = Result<Record, String>> + Send,
    F: Fn(&mut Record) -> Result<(), String> + Sync,
    W: FnMut(&Record) -> Result<(), String>,
{
    assert!(workers > 0 && batch_size > 0);
    thread::scope(|scope| {
        let mut inputs = Vec::with_capacity(workers);
        let mut outputs = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);

        for _ in 0..workers {
            let (input_tx, input_rx) = mpsc::sync_channel::<(usize, Vec<Record>)>(1);
            let (output_tx, output_rx) = mpsc::sync_channel::<(usize, Vec<Record>)>(1);
            inputs.push(input_tx);
            outputs.push(output_rx);
            let transform = &transform;
            handles.push(scope.spawn(move || -> Result<(), String> {
                while let Ok((index, mut batch)) = input_rx.recv() {
                    for record in &mut batch {
                        transform(record)?;
                    }
                    output_tx.send((index, batch))
                        .map_err(|_| "Tag output channel disconnected".to_string())?;
                }
                Ok(())
            }));
        }

        let reader = scope.spawn(move || -> Result<usize, String> {
            let mut index = 0;
            let mut batch = Vec::with_capacity(batch_size);
            for record in records {
                batch.push(record?);
                if batch.len() == batch_size {
                    let ready = std::mem::replace(&mut batch, Vec::with_capacity(batch_size));
                    inputs[index % workers].send((index, ready))
                        .map_err(|_| "Tag input channel disconnected".to_string())?;
                    index += 1;
                }
            }
            if !batch.is_empty() {
                inputs[index % workers].send((index, batch))
                    .map_err(|_| "Tag input channel disconnected".to_string())?;
                index += 1;
            }
            Ok(index)
        });

        let mut expected = 0;
        let written = (|| -> Result<(), String> {
            // The first disconnected channel is EOF because dispatch is round-robin.
            // Reader/worker joins below distinguish EOF from a failed stage.
            while let Ok((index, batch)) = outputs[expected % workers].recv() {
                if index != expected {
                    return Err("Tag batch order mismatch".to_string());
                }
                for record in &batch {
                    write(record)?;
                }
                expected += 1;
            }
            Ok(())
        })();

        // Release blocked sends before joining, including on a write failure.
        drop(outputs);
        let mut failure = written.err();
        for handle in handles {
            let result = handle.join()
                .unwrap_or_else(|_| Err("Tag calculation worker panicked".to_string()));
            if failure.is_none() {
                failure = result.err();
            }
        }
        match reader.join() {
            Ok(Ok(dispatched)) => {
                if failure.is_none() && dispatched != expected {
                    failure = Some("Tag pipeline did not write every batch".to_string());
                }
            }
            Ok(Err(error)) => {
                if failure.is_none() {
                    failure = Some(error);
                }
            }
            Err(_) => {
                if failure.is_none() {
                    failure = Some("Tag reader thread panicked".to_string());
                }
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    })
}

/// Stream ordered batches through XM workers and HTSlib compression workers.
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

    let (decompress, workers, compress, _) = thread_allocation(threads);
    if decompress > 0 {
        reader.set_threads(decompress).expect("Failed to set BAM decompression threads");
    }
    if compress > 0 {
        writer.set_threads(compress).expect("Failed to set BAM compression threads");
    }

    let refgenome = Arc::new(refgenome);
    let tid2size = Arc::new(tid2size);
    let rcmapping = Arc::new(rcmapping);
    let result = stream_batches(
        reader.records().map(|r| r.map_err(|e| format!("Error reading BAM record: {}", e))),
        workers,
        BATCH_SIZE,
        |record| {
            let xm = determine_xm_tag_string(
                record, &refgenome, &tid2size, &rcmapping, is_paired_end,
            );
            record.push_aux(b"XM", Aux::String(&xm))
                .map_err(|e| format!("Error adding XM tag: {}", e))
        },
        |record| writer.write(record).map_err(|e| format!("Error writing BAM record: {}", e)),
    );
    // Finish HTSlib's background work before reporting completion.
    drop(reader);
    drop(writer);
    result.unwrap_or_else(|error| panic!("{}", error));
    println!("Done writing!");
}

pub fn run(input: &str, output: &str, genome: &str, threads: usize) {
    let (decompress, workers, compress, overhead) = thread_allocation(threads);
    eprintln!(
        "tag threads: budget={}, decompression={}, calculation={}, compression={}, main/I/O={}",
        threads, decompress, workers, compress, overhead,
    );
    if threads <= 2 {
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
        let dir = TestDir::new();
        run(
            "tests/test1.bam",
            dir.0.join("missing/out.bam").to_str().unwrap(),
            "tests/tinyref.fa",
            1,
        )
    }
    #[test]
    #[should_panic]
    fn error_when_reference_genome_is_not_found() {
        let dir = TestDir::new();
        run(
            "tests/test1.bam",
            dir.0.join("out.bam").to_str().unwrap(),
            dir.0.join("missing.fa").to_str().unwrap(),
            1,
        )
    }

    struct TestDir(std::path::PathBuf);

    impl TestDir {
        fn new() -> Self {
            static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!("metheor-tag-{}-{}", std::process::id(), id));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn record(index: usize) -> Record {
        let mut r = Record::new();
        r.set(
            index.to_string().as_bytes(),
            Some(&bam::record::CigarString(vec![Cigar::Match(8)])),
            b"CGCGCGCG",
            &[30; 8],
        );
        r.set_tid(0);
        r.set_pos(((index % 80) * 2) as i64);
        r.set_flags([99, 147, 83, 163][index % 4]);
        r.set_mapq(30);
        r.set_mtid(0);
        r.set_mpos(r.pos() + 2);
        r.set_insert_size(16);
        r.push_aux(b"RG", Aux::String("test")).unwrap();
        r
    }

    fn decoded(path: &std::path::Path) -> (Vec<u8>, Vec<Record>) {
        let mut reader = bam::Reader::from_path(path).unwrap();
        let header = reader.header().as_bytes().to_vec();
        let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        (header, records)
    }

    #[test]
    fn test_thread_allocations() {
        for budget in 1..=100 {
            let (d, w, c, overhead) = thread_allocation(budget);
            assert!(d + w + c + overhead <= budget);
            if budget >= 3 {
                assert_eq!(d + w + c + overhead, budget);
                assert!(w > 0);
            }
            if budget >= 8 {
                assert!(d > 0 && c > 0);
                assert_eq!(overhead, 4);
            }
        }
        assert_eq!(thread_allocation(100), (8, 44, 44, 4));
    }

    #[test]
    fn test_mt_matches_st() {
        let dir = TestDir::new();
        let input = dir.0.join("input.bam");
        let genome = dir.0.join("ref.fa");
        let st = dir.0.join("single.bam");
        let mt = dir.0.join("parallel.bam");
        std::fs::write(&genome, format!(">chr1\n{}\n", "CG".repeat(128))).unwrap();
        std::fs::write(dir.0.join("ref.fa.fai"), "chr1\t256\t6\t256\t257\n").unwrap();
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ").push_tag(b"SN", "chr1").push_tag(b"LN", 256),
        );

        // Empty BAM, a partial batch, and multiple full batches plus a partial batch.
        for count in [0, 7, BATCH_SIZE * 2 + 3] {
            let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
            for i in 0..count {
                writer.write(&record(i)).unwrap();
            }
            drop(writer);
            run(input.to_str().unwrap(), st.to_str().unwrap(), genome.to_str().unwrap(), 1);
            let expected = decoded(&st);
            assert_eq!(expected.1.len(), count);
            for budget in [2, 3, 8, 16] {
                run(input.to_str().unwrap(), mt.to_str().unwrap(), genome.to_str().unwrap(), budget);
                assert_eq!(expected, decoded(&mt), "budget={}, records={}", budget, count);
            }
        }
    }

    #[test]
    fn test_slow_first_worker_preserves_order() {
        let mut written = Vec::new();
        stream_batches(
            (0..19).map(|i| Ok(record(i))),
            3,
            2,
            |r| {
                if r.qname() == b"0" {
                    std::thread::sleep(std::time::Duration::from_millis(25));
                }
                Ok(())
            },
            |r| {
                written.push(String::from_utf8(r.qname().to_vec()).unwrap());
                Ok(())
            },
        ).unwrap();
        assert_eq!(written, (0..19).map(|i| i.to_string()).collect::<Vec<_>>());
    }

    #[test]
    fn test_failures_disconnect_and_join() {
        // A timeout makes a channel shutdown regression fail instead of hanging the suite.
        for failure in ["read", "worker", "panic", "write"] {
            let (tx, rx) = mpsc::channel();
            std::thread::spawn(move || {
                let result = stream_batches(
                    (0..50).map(|i| {
                        if failure == "read" && i == 9 {
                            Err("injected read failure".to_string())
                        } else {
                            Ok(record(i))
                        }
                    }),
                    3,
                    2,
                    |r| {
                        if r.qname() == b"3" {
                            if failure == "panic" {
                                panic!("injected worker panic");
                            }
                            if failure == "worker" {
                                return Err("injected worker failure".to_string());
                            }
                        }
                        Ok(())
                    },
                    |_| {
                        if failure == "write" {
                            Err("injected write failure".to_string())
                        } else {
                            Ok(())
                        }
                    },
                );
                tx.send(result.is_err()).unwrap();
            });
            assert!(rx.recv_timeout(std::time::Duration::from_secs(5)).unwrap(), "{}", failure);
        }
    }
}
