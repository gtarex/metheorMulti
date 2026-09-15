use rust_htslib::{faidx, htslib};
use rust_htslib::{
    bam,
    bam::ext::BamRecordExtensions,
    bam::record::{Aux, Cigar, Record},
    bam::Read,
};
use std::cmp::{max, min};
use std::collections::{HashMap, VecDeque};
use std::ffi::CString;
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
    is_top_strand: bool,
    is_unmapped: bool,
    seq_len: usize,
}

impl ReadPayload {
    fn from_record(r: &Record) -> Self {
        if r.is_unmapped() {
            Self {
                tid: -1,
                reference_start: -1,
                reference_end: -1,
                seq_bytes: Vec::new(),
                cigar_ops: Vec::new(),
                is_top_strand: true,
                is_unmapped: true,
                seq_len: r.seq_len(),
            }
        } else {
            Self {
                tid: r.tid(),
                reference_start: r.reference_start(),
                reference_end: r.reference_end(),
                seq_bytes: r.seq().as_bytes().to_vec(),
                cigar_ops: r.cigar().iter().copied().collect(),
                is_top_strand: crate::readutil::is_top_strand(r),
                is_unmapped: false,
                seq_len: r.seq_len(),
            }
        }
    }
}

/// Compute XM tag string from a ReadPayload (used by worker threads).
fn compute_xm_tag_from_payload(
    payload: &ReadPayload,
    refgenome: &HashMap<usize, Vec<u8>>,
    tid2size: &HashMap<usize, usize>,
    rcmapping: &HashMap<char, char>,
) -> String {
    if payload.is_unmapped {
        return ".".repeat(payload.seq_len);
    }
    let tid = payload.tid;
    let start = payload.reference_start;
    let end = payload.reference_end;

    let flag_reverse_complement = !payload.is_top_strand;

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

    // CIGAR maps read bases to genomic positions. Context always comes from
    // adjacent reference bases, never from a reference with read gaps removed.
    let read_bases = read_seq.as_bytes();
    let mut xm = vec![b'.'; read_bases.len()];
    let mut read_pos = 0usize;
    // ref_seq has two upstream bases (or N padding) before the alignment.
    let mut ref_pos = 2usize;

    for cigar in &payload.cigar_ops {
        match *cigar {
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                for offset in 0..length as usize {
                    let query_pos = read_pos + offset;
                    let genome_pos = ref_pos + offset;
                    let base = match (flag_reverse_complement, read_bases[query_pos]) {
                        (false, b'C') | (true, b'G') => b'C',
                        (false, b'T') | (true, b'A') => b'T',
                        _ => continue,
                    };
                    let context = if flag_reverse_complement {
                        reverse_complement(&ref_seq[genome_pos - 2..genome_pos + 1], rcmapping)
                    } else {
                        ref_seq[genome_pos..genome_pos + 3].to_string()
                    };
                    if !context.starts_with('C') {
                        continue;
                    }
                    let code = if context.starts_with("CG") {
                        b'Z'
                    } else if is_chg_context(&context) {
                        b'X'
                    } else if is_chh_context(&context) {
                        b'H'
                    } else if is_unknown_context(&context) {
                        b'U'
                    } else {
                        // Ambiguous reference context: retain a dot at this base.
                        b'.'
                    };
                    xm[query_pos] = if base == b'T' {
                        code.to_ascii_lowercase()
                    } else {
                        code
                    };
                }
                read_pos += length as usize;
                ref_pos += length as usize;
            }
            Cigar::Ins(length) | Cigar::SoftClip(length) => {
                // These read bases have no reference position; leave dots.
                read_pos += length as usize;
            }
            Cigar::Del(length) | Cigar::RefSkip(length) => {
                ref_pos += length as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    assert_eq!(read_pos, read_bases.len(), "CIGAR/read length mismatch");
    String::from_utf8(xm).expect("XM tags contain only ASCII characters")
}

/// Original determine_xm_tag_string (kept for backward compatibility and single-threaded path).
/// Delegates to compute_xm_tag_from_payload after extracting fields from the Record.
pub fn determine_xm_tag_string(
    r: &Record,
    refgenome: &HashMap<usize, Vec<u8>>,
    tid2size: &HashMap<usize, usize>,
    rcmapping: &HashMap<char, char>,
    _is_paired_end: bool,
) -> String {
    // Unmapped records have no reference context, even if a mate position is set.
    if r.is_unmapped() {
        return ".".repeat(r.seq_len());
    }
    let payload = ReadPayload::from_record(r);
    compute_xm_tag_from_payload(&payload, refgenome, tid2size, rcmapping)
}

/// BAM output with a checked final close. rust-htslib's Writer only closes in
/// Drop, which discards errors from the final BGZF flush and EOF write.
struct CheckedBamWriter {
    file: *mut htslib::htsFile,
    header: bam::HeaderView,
}

impl CheckedBamWriter {
    fn from_path(output: &str, header: &bam::Header) -> Result<Self, String> {
        let path = CString::new(output).map_err(|e| format!("Invalid BAM output path: {}", e))?;
        // Use the same header construction and default BAM compression as bam::Writer.
        let header = bam::HeaderView::from_header(header);
        // SAFETY: both strings are NUL-terminated; this writer owns the returned handle.
        let file = unsafe { htslib::hts_open(path.as_ptr(), b"wb\0".as_ptr().cast()) };
        if file.is_null() {
            return Err(format!("Cannot open BAM output {}: {}", output, std::io::Error::last_os_error()));
        }
        let writer = Self { file, header };
        // SAFETY: the file is open and the header remains owned by this writer.
        if unsafe { htslib::sam_hdr_write(writer.file, writer.header.inner_ptr()) } < 0 {
            return Err("Error writing BAM header".to_string());
        }
        Ok(writer)
    }

    fn set_threads(&mut self, threads: usize) -> Result<(), String> {
        if threads == 0 || threads > std::os::raw::c_int::MAX as usize {
            return Err("Invalid BAM compression thread count".to_string());
        }
        // SAFETY: the handle is open and the worker count fits HTSlib's integer type.
        if unsafe { htslib::hts_set_threads(self.file, threads as std::os::raw::c_int) } != 0 {
            return Err("Failed to set BAM compression threads".to_string());
        }
        Ok(())
    }

    fn write(&mut self, record: &Record) -> Result<(), String> {
        // SAFETY: the writer, header and borrowed record all remain valid for this call.
        if unsafe { htslib::sam_write1(self.file, self.header.inner_ptr(), record.inner()) } < 0 {
            return Err("Error writing BAM record".to_string());
        }
        Ok(())
    }

    fn finish(mut self) -> Result<(), String> {
        // hts_close consumes the handle even on failure. Clear it before closing
        // so Drop cannot close it again; the header stays alive through the call.
        let file = std::mem::replace(&mut self.file, std::ptr::null_mut());
        // SAFETY: this is the single close of our owned, open handle.
        if unsafe { htslib::hts_close(file) } != 0 {
            return Err(format!("Error finalizing BAM output: {}", std::io::Error::last_os_error()));
        }
        Ok(())
    }
}

impl Drop for CheckedBamWriter {
    fn drop(&mut self) {
        if !self.file.is_null() {
            // Fallback cleanup on an earlier error/panic. Success paths use finish().
            // SAFETY: a non-null handle here has not yet been closed.
            unsafe { htslib::hts_close(self.file); }
        }
    }
}

/// Single-threaded implementation. Reads BAM, processes each record sequentially,
/// writes output BAM with XM tag.
fn run_helper(input: &str, output: &str, genome: &str) {
    let mut reader = bamutil::get_reader(input);
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
    let mut writer = match CheckedBamWriter::from_path(output, &header_tmpl) {
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
            determine_xm_tag_string(&r, &refgenome, &tid2size, &rcmapping, r.is_paired());
        // Attach XM tag to the record.
        if r.aux(b"XM").is_ok() {
            r.remove_aux(b"XM").expect("Error replacing existing XM tag");
        }
        let add_result = r.push_aux("XM".as_bytes(), Aux::String(&xm_tag_string));
        match add_result {
            Ok(_) => (),
            Err(e) => panic!("Error adding XM tag to alignment record. {}", e),
        }
        // Write record to output.
        writer.write(&r).expect("Error writing to output file.");
    }
    writer.finish().unwrap_or_else(|error| panic!("{}", error));
}

const BATCH_SIZE: usize = 8192;

// Decompression, calculation, compression, main/I/O threads.
fn thread_allocation(threads: usize) -> (usize, usize, usize, usize) {
    assert!(threads >= 1, "threads must be at least 1");
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
#[cfg(test)]
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

/// Stream ordered batches through XM calculation workers without crossing BAM records across threads.
/// Record creation, modification, writing, and dropping remain strictly confined to the main I/O thread.
fn stream_tag_records(
    reader: &mut bam::Reader,
    writer: &mut CheckedBamWriter,
    refgenome: &Arc<HashMap<usize, Vec<u8>>>,
    tid2size: &Arc<HashMap<usize, usize>>,
    rcmapping: &Arc<HashMap<char, char>>,
    workers: usize,
    batch_size: usize,
) -> Result<(), String> {
    assert!(workers > 0 && batch_size > 0);
    thread::scope(|scope| {
        let mut inputs = Vec::with_capacity(workers);
        let mut outputs = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);

        for _ in 0..workers {
            let (input_tx, input_rx) = mpsc::sync_channel::<(usize, Vec<ReadPayload>)>(1);
            let (output_tx, output_rx) = mpsc::sync_channel::<(usize, Vec<String>)>(1);
            inputs.push(input_tx);
            outputs.push(output_rx);
            let refgenome = Arc::clone(refgenome);
            let tid2size = Arc::clone(tid2size);
            let rcmapping = Arc::clone(rcmapping);
            handles.push(scope.spawn(move || -> Result<(), String> {
                while let Ok((index, payloads)) = input_rx.recv() {
                    let mut tags = Vec::with_capacity(payloads.len());
                    for payload in &payloads {
                        tags.push(compute_xm_tag_from_payload(
                            payload,
                            &refgenome,
                            &tid2size,
                            &rcmapping,
                        ));
                    }
                    output_tx
                        .send((index, tags))
                        .map_err(|_| "Tag output channel disconnected".to_string())?;
                }
                Ok(())
            }));
        }

        let max_in_flight = workers.max(1) * 2;
        let mut in_flight: VecDeque<(usize, Vec<Record>)> = VecDeque::new();
        let mut dispatched = 0usize;
        let mut expected = 0usize;
        let mut current_batch = Vec::with_capacity(batch_size);

        let write_next_batch = |expected_idx: usize,
                                in_flight: &mut VecDeque<(usize, Vec<Record>)>,
                                outputs: &[mpsc::Receiver<(usize, Vec<String>)>],
                                writer: &mut CheckedBamWriter|
         -> Result<(), String> {
            let (index, tags) = outputs[expected_idx % workers]
                .recv()
                .map_err(|_| "Tag output channel disconnected unexpectedly".to_string())?;
            if index != expected_idx {
                return Err(format!(
                    "Tag batch order mismatch: expected {}, got {}",
                    expected_idx, index
                ));
            }
            let (rec_index, mut records) = in_flight
                .pop_front()
                .ok_or_else(|| "In-flight record queue underflow".to_string())?;
            if rec_index != expected_idx {
                return Err(format!(
                    "Record batch order mismatch: expected {}, got {}",
                    expected_idx, rec_index
                ));
            }
            if records.len() != tags.len() {
                return Err(format!(
                    "Batch length mismatch: {} records vs {} tags",
                    records.len(),
                    tags.len()
                ));
            }
            for (record, xm) in records.iter_mut().zip(&tags) {
                if record.aux(b"XM").is_ok() {
                    record
                        .remove_aux(b"XM")
                        .map_err(|e| format!("Error replacing XM tag: {}", e))?;
                }
                record
                    .push_aux(b"XM", Aux::String(xm))
                    .map_err(|e| format!("Error adding XM tag: {}", e))?;
                writer
                    .write(record)
                    .map_err(|e| format!("Error writing BAM record: {}", e))?;
            }
            Ok(())
        };

        let run_res = (|| -> Result<(), String> {
            for record_res in reader.records() {
                let record = record_res.map_err(|e| format!("Error reading BAM record: {}", e))?;
                current_batch.push(record);

                if current_batch.len() == batch_size {
                    while in_flight.len() >= max_in_flight {
                        write_next_batch(expected, &mut in_flight, &outputs, writer)?;
                        expected += 1;
                    }
                    let batch_records =
                        std::mem::replace(&mut current_batch, Vec::with_capacity(batch_size));
                    let payloads: Vec<ReadPayload> =
                        batch_records.iter().map(ReadPayload::from_record).collect();
                    inputs[dispatched % workers]
                        .send((dispatched, payloads))
                        .map_err(|_| "Tag input channel disconnected".to_string())?;
                    in_flight.push_back((dispatched, batch_records));
                    dispatched += 1;
                }
            }

            if !current_batch.is_empty() {
                while in_flight.len() >= max_in_flight {
                    write_next_batch(expected, &mut in_flight, &outputs, writer)?;
                    expected += 1;
                }
                let batch_records = current_batch;
                let payloads: Vec<ReadPayload> =
                    batch_records.iter().map(ReadPayload::from_record).collect();
                inputs[dispatched % workers]
                    .send((dispatched, payloads))
                    .map_err(|_| "Tag input channel disconnected".to_string())?;
                in_flight.push_back((dispatched, batch_records));
                dispatched += 1;
            }

            while expected < dispatched {
                write_next_batch(expected, &mut in_flight, &outputs, writer)?;
                expected += 1;
            }

            Ok(())
        })();

        // Drop inputs to signal workers to terminate gracefully.
        drop(inputs);

        // Join workers and collect any error
        let mut failure = run_res.err();
        for handle in handles {
            let res = handle
                .join()
                .unwrap_or_else(|_| Err("Tag calculation worker panicked".to_string()));
            if failure.is_none() {
                failure = res.err();
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
    let mut writer = match CheckedBamWriter::from_path(output, &header_tmpl) {
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
    let result = stream_tag_records(
        &mut reader,
        &mut writer,
        &refgenome,
        &tid2size,
        &rcmapping,
        workers,
        BATCH_SIZE,
    );
    // Finish HTSlib's background work before reporting completion.
    drop(reader);
    let finished = writer.finish();
    result.unwrap_or_else(|error| panic!("{}", error));
    finished.unwrap_or_else(|error| panic!("{}", error));
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

    include!("tag_context_tests.rs");

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
    fn test_checked_writer_matches_bam_writer() {
        let dir = TestDir::new();
        let original = dir.0.join("original.bam");
        let checked = dir.0.join("checked.bam");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"HD").push_tag(b"VN", "1.6").push_tag(b"SO", "unsorted"),
        );
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ").push_tag(b"SN", "chr1").push_tag(b"LN", 256),
        );
        header.push_record(
            bam::header::HeaderRecord::new(b"RG").push_tag(b"ID", "test").push_tag(b"SM", "sample"),
        );
        header.push_record(
            bam::header::HeaderRecord::new(b"PG").push_tag(b"ID", "test").push_tag(b"PN", "metheor"),
        );
        // Empty output and enough records to span several BGZF blocks.
        for count in [0, 2049] {
            for threads in [0, 2] {
                let mut expected = bam::Writer::from_path(&original, &header, bam::Format::Bam).unwrap();
                let mut actual = CheckedBamWriter::from_path(checked.to_str().unwrap(), &header).unwrap();
                if threads > 0 {
                    expected.set_threads(threads).unwrap();
                    actual.set_threads(threads).unwrap();
                }
                for i in 0..count {
                    let mut r = record(i);
                    r.push_aux(b"XM", Aux::String("ZzZzZzZz")).unwrap();
                    expected.write(&r).unwrap();
                    actual.write(&r).unwrap();
                }
                drop(expected);
                actual.finish().unwrap();
                let expected = decoded(&original);
                assert_eq!(expected.1.len(), count);
                assert_eq!(expected, decoded(&checked), "threads={}, records={}", threads, count);
            }
        }
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_checked_writer_reports_final_flush_failure() {
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ").push_tag(b"SN", "chr1").push_tag(b"LN", 256),
        );
        // The small header/record fit in the buffers; /dev/full fails at final flush.
        let mut writer = CheckedBamWriter::from_path("/dev/full", &header).unwrap();
        writer.write(&record(0)).unwrap();
        let error = writer.finish().unwrap_err();
        assert!(error.contains("Error finalizing BAM output"), "{}", error);
    }

    #[test]
    fn test_thread_allocations() {
        for budget in 1..=256 {
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
        assert_eq!(thread_allocation(128), (8, 58, 58, 4));
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
