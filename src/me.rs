use rust_htslib::{bam, bam::Read};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::str;
use std::sync::mpsc;
use std::thread;

use crate::{bamutil, progressbar, readutil};

pub struct QuartetStat {
    pos1: readutil::CpGPosition,
    pos2: readutil::CpGPosition,
    pos3: readutil::CpGPosition,
    pos4: readutil::CpGPosition,
    quartet_pattern_counts: [u32; 16],
}

impl QuartetStat {
    fn new(q: readutil::Quartet) -> Self {
        let pos1 = q.pos1;
        let pos2 = q.pos2;
        let pos3 = q.pos3;
        let pos4 = q.pos4;

        let quartet_pattern_counts = [0; 16];
        Self {
            pos1,
            pos2,
            pos3,
            pos4,
            quartet_pattern_counts,
        }
    }

    fn get_read_depth(&self) -> u32 {
        self.quartet_pattern_counts.iter().sum()
    }

    fn add_quartet_pattern(&mut self, p: readutil::QuartetPattern) {
        self.quartet_pattern_counts[p] += 1;
    }

    fn merge(&mut self, other: &QuartetStat) {
        for i in 0..16 {
            self.quartet_pattern_counts[i] += other.quartet_pattern_counts[i];
        }
    }

    fn compute_me(&self) -> f32 {
        let mut me: f32 = 0.0;

        let total: u32 = self.quartet_pattern_counts.iter().sum();
        for count in self.quartet_pattern_counts.iter() {
            let p: f32 = (*count as f32) / (total as f32);
            if *count > 0 {
                me += p * p.log2();
            }
        }
        me *= -0.25;

        me
    }

    fn to_bedgraph_field(&self, header: &bam::HeaderView) -> String {
        let chrom = bamutil::tid2chrom(self.pos1.tid, header);
        let me = self.compute_me();
        let c = self.quartet_pattern_counts;

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            chrom, self.pos1.pos, self.pos2.pos, self.pos3.pos, self.pos4.pos, me,
            c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
            c[8], c[9], c[10], c[11], c[12], c[13], c[14], c[15],
        )
    }
}

const BATCH_SIZE: usize = 8192;
const OUTPUT_BUFFER_SIZE: usize = 64 * 1024 * 1024;

// Decompression, calculation, main/I/O threads.
fn thread_allocation(threads: usize) -> (usize, usize, usize) {
    assert!((1..=100).contains(&threads), "threads must be from 1 to 100");
    match threads {
        1 => (0, 0, 1),
        2..=3 => (0, threads - 1, 1),
        _ => {
            let decompress = ((threads - 2) / 10).max(1).min(8);
            (decompress, threads - 2 - decompress, 2)
        }
    }
}

fn merge_maps(
    mut a: HashMap<readutil::Quartet, QuartetStat>,
    b: HashMap<readutil::Quartet, QuartetStat>,
) -> HashMap<readutil::Quartet, QuartetStat> {
    for (q, stat) in b {
        a.entry(q)
            .or_insert_with(|| QuartetStat::new(q))
            .merge(&stat);
    }
    a
}

pub fn compute(
    input: &str,
    output: &str,
    min_depth: u32,
    min_qual: u8,
    cpg_set: &Option<String>,
    threads: usize,
) {
    let (decompress, workers, overhead) = thread_allocation(threads);
    eprintln!(
        "me threads: budget={}, decompression={}, calculation={}, main/I/O={}",
        threads, decompress, workers, overhead,
    );
    let header = {
        let reader = bamutil::get_reader(input);
        bamutil::get_header(&reader)
    };

    let result = if threads <= 1 {
        compute_helper(input, min_qual, cpg_set)
    } else {
        compute_helper_mt(input, min_qual, cpg_set, threads)
    };

    let file = fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();
    let mut out = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

    let mut to_write: Vec<&QuartetStat> = result
        .values()
        .filter(|stat| stat.get_read_depth() >= min_depth)
        .collect();
    to_write.sort_by_key(|s| (s.pos1.tid, s.pos1.pos, s.pos2.pos, s.pos3.pos, s.pos4.pos));

    for stat in to_write {
        writeln!(out, "{}", stat.to_bedgraph_field(&header))
            .expect("Error writing to output file.");
    }
    out.flush().expect("Error flushing entropy output.");
}

/// Original single-threaded implementation (unchanged, used when threads <= 1).
pub fn compute_helper(
    input: &str,
    min_qual: u8,
    cpg_set: &Option<String>,
) -> HashMap<readutil::Quartet, QuartetStat> {
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = &readutil::get_target_cpgs(cpg_set, &header);
    let mut quartet2stat: HashMap<readutil::Quartet, QuartetStat> = HashMap::new();

    let mut readcount = 0;
    let mut valid_readcount = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        if let Some(target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs);
        }

        readcount += 1;

        if r.mapq() < min_qual {
            continue;
        }
        valid_readcount += 1;

        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
        for (q, p) in quartets.iter().zip(patterns.iter()) {
            let stat = quartet2stat.entry(*q).or_insert(QuartetStat::new(*q));

            stat.add_quartet_pattern(*p);
        }

        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }
    quartet2stat
}

/// Parse XM and accumulate counts inside workers; never filter batch-local depth.
fn count_batches<I>(
    records: I,
    min_qual: u8,
    target_cpgs: &Option<HashSet<readutil::CpGPosition>>,
    workers: usize,
    batch_size: usize,
) -> Result<HashMap<readutil::Quartet, QuartetStat>, String>
where
    I: Iterator<Item = Result<bam::Record, String>>,
{
    assert!(workers > 0 && batch_size > 0);
    thread::scope(|scope| {
        let mut inputs = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<Vec<bam::Record>>(1);
            inputs.push(tx);
            handles.push(scope.spawn(move || {
                let mut local_map: HashMap<readutil::Quartet, QuartetStat> = HashMap::new();
                while let Ok(batch) = rx.recv() {
                    for record in batch {
                        // Preserve the single-thread path's parsing/filter order.
                        let mut br = readutil::BismarkRead::new(&record);
                        if let Some(target) = target_cpgs {
                            br.filter_isin(target);
                        }
                        if record.mapq() < min_qual {
                            continue;
                        }
                        let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
                        for (q, p) in quartets.iter().zip(patterns.iter()) {
                            local_map.entry(*q).or_insert_with(|| QuartetStat::new(*q))
                                .add_quartet_pattern(*p);
                        }
                    }
                }
                local_map
            }));
        }

        let dispatched = (|| -> Result<(), String> {
            let mut worker = 0;
            let mut batch = Vec::with_capacity(batch_size);
            for record in records {
                batch.push(record?);
                if batch.len() == batch_size {
                    let ready = std::mem::replace(&mut batch, Vec::with_capacity(batch_size));
                    inputs[worker].send(ready)
                        .map_err(|_| "ME input channel disconnected".to_string())?;
                    worker = (worker + 1) % workers;
                }
            }
            if !batch.is_empty() {
                inputs[worker].send(batch)
                    .map_err(|_| "ME input channel disconnected".to_string())?;
            }
            Ok(())
        })();

        // Closing inputs also releases waiting workers when BAM reading fails.
        drop(inputs);
        let mut failure = dispatched.err();
        let mut result = HashMap::new();
        for handle in handles {
            match handle.join() {
                Ok(local_map) => {
                    if failure.is_none() {
                        result = merge_maps(result, local_map);
                    }
                }
                Err(_) => {
                    if failure.is_none() {
                        failure = Some("ME calculation worker panicked".to_string());
                    }
                }
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(result),
        }
    })
}

/// Overlap BAM input/decompression with batched parsing and quartet counting.
pub fn compute_helper_mt(
    input: &str,
    min_qual: u8,
    cpg_set: &Option<String>,
    threads: usize,
) -> HashMap<readutil::Quartet, QuartetStat> {
    let (decompress, workers, _) = thread_allocation(threads);
    if workers == 0 {
        return compute_helper(input, min_qual, cpg_set);
    }
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);
    let target_cpgs = readutil::get_target_cpgs(cpg_set, &header);
    if decompress > 0 {
        reader.set_threads(decompress).expect("Failed to set BAM decompression threads");
    }

    let mut readcount = 0;
    let mut valid_readcount = 0;
    let bar = progressbar::ProgressBar::new();
    let records = reader.records().map(|r| {
        let record = r.map_err(|e| format!("Error reading BAM record: {}", e))?;
        readcount += 1;
        if record.mapq() >= min_qual {
            valid_readcount += 1;
        }
        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount);
        }
        Ok(record)
    });
    let result = count_batches(records, min_qual, &target_cpgs, workers, BATCH_SIZE);
    drop(reader);
    result.unwrap_or_else(|error| panic!("{}", error))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let input = "tests/test1.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.get_read_depth(), 16);
            assert_eq!(reads.compute_me(), 1.0);
        }
    }

    #[test]
    fn test2() {
        let input = "tests/test2.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 0.25);
        }
    }
    #[test]
    fn test3() {
        let input = "tests/test3.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 1);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 0.25);
        }
    }
    #[test]
    fn test4() {
        let input = "tests/test4.bam";
        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 2);

        for (_, reads) in quartet2stat.iter() {
            assert_eq!(reads.compute_me(), 1.0);
        }
    }
    #[test]
    fn test5() {
        // No reads pass quality cutoff.
        let input = "tests/test5.bam";

        let min_qual = 10;
        let cpg_set = None;

        let quartet2stat = compute_helper(input, min_qual, &cpg_set);

        assert_eq!(quartet2stat.len(), 0);
    }

    fn assert_counts_equal(
        expected: &HashMap<readutil::Quartet, QuartetStat>,
        actual: &HashMap<readutil::Quartet, QuartetStat>,
    ) {
        assert_eq!(expected.len(), actual.len());
        for (q, stat) in expected {
            let other = actual.get(q).expect("Quartet missing from parallel result");
            assert_eq!(stat.quartet_pattern_counts, other.quartet_pattern_counts);
            assert_eq!(stat.compute_me().to_bits(), other.compute_me().to_bits());
        }
    }

    #[test]
    fn test_mt_matches_st() {
        for input in ["tests/test1.bam", "tests/test2.bam", "tests/test3.bam", "tests/test4.bam", "tests/test5.bam"] {
            let expected = compute_helper(input, 10, &None);
            for threads in [2, 3, 4, 8] {
                let actual = compute_helper_mt(input, 10, &None, threads);
                assert_counts_equal(&expected, &actual);
            }
        }
    }

    #[test]
    fn test_thread_allocations() {
        for budget in 1..=100 {
            let (d, w, overhead) = thread_allocation(budget);
            assert_eq!(d + w + overhead, budget);
            if budget > 1 {
                assert!(w > 0);
            }
            if budget >= 4 {
                assert!(d > 0);
                assert_eq!(overhead, 2);
            }
        }
        assert_eq!(thread_allocation(100), (8, 90, 2));
    }

    struct TestDir(std::path::PathBuf);

    impl TestDir {
        fn new() -> Self {
            static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!("metheor-me-{}-{}", std::process::id(), id));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn record(pos: i64, mapq: u8, index: usize) -> bam::Record {
        use bam::record::{Aux, Cigar, CigarString};
        let mut r = bam::Record::new();
        r.set(
            format!("pair{}", index / 2).as_bytes(),
            Some(&CigarString(vec![Cigar::Match(8)])),
            b"CGCGCGCG",
            &[30; 8],
        );
        r.set_tid(0);
        r.set_pos(pos);
        r.set_mapq(mapq);
        r.set_flags(if index % 2 == 0 { 99 } else { 147 });
        r.set_mtid(0);
        r.set_mpos(pos);
        let mut xm = String::new();
        for bit in (0..4).rev() {
            xm.push(if index & (1 << bit) != 0 { 'Z' } else { 'z' });
            xm.push('.');
        }
        r.push_aux(b"XM", Aux::String(&xm)).unwrap();
        r.push_aux(b"XG", Aux::String("CT")).unwrap();
        r
    }

    fn depth_records() -> Vec<bam::Record> {
        let mut records = Vec::new();
        // Interleaved windows; some mate alignments overlap and must count separately.
        for index in 0..11 {
            for (pos, depth) in [(0, 9), (20, 10), (40, 11)] {
                if index < depth {
                    records.push(record(pos, 10, index));
                }
            }
        }
        for pos in [0, 20, 40] {
            records.push(record(pos, 9, 15));
        }
        records
    }

    fn write_fixture(path: &std::path::Path, records: &[bam::Record]) {
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ").push_tag(b"SN", "chr1").push_tag(b"LN", 100),
        );
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        for record in records {
            writer.write(record).unwrap();
        }
    }

    #[test]
    fn test_depth_and_targets_across_batches() {
        let dir = TestDir::new();
        let input = dir.0.join("input.bam");
        write_fixture(&input, &depth_records());
        let expected = compute_helper(input.to_str().unwrap(), 10, &None);
        let actual = count_batches(depth_records().into_iter().map(Ok), 10, &None, 3, 2).unwrap();
        assert_counts_equal(&expected, &actual);
        let mut depths: Vec<_> = actual.values().map(|s| s.get_read_depth()).collect();
        depths.sort_unstable();
        assert_eq!(depths, vec![9, 10, 11]);

        let bed = dir.0.join("target.bed");
        std::fs::write(&bed, "chr1\t20\t21\nchr1\t22\t23\nchr1\t24\t25\nchr1\t26\t27\n").unwrap();
        let target_path = Some(bed.to_str().unwrap().to_string());
        let reader = bamutil::get_reader(input.to_str().unwrap());
        let target = readutil::get_target_cpgs(&target_path, reader.header());
        let expected = compute_helper(input.to_str().unwrap(), 10, &target_path);
        let actual = count_batches(depth_records().into_iter().map(Ok), 10, &target, 3, 2).unwrap();
        assert_counts_equal(&expected, &actual);
        assert_eq!(actual.len(), 1);
    }

    #[test]
    fn test_complete_text_and_batch_boundaries() {
        let dir = TestDir::new();
        let input = dir.0.join("input.bam");
        let single = dir.0.join("single.tsv");
        let parallel = dir.0.join("parallel.tsv");
        let base = depth_records();
        let repeated: Vec<_> = base.iter().cycle().take(BATCH_SIZE * 2 + 3).cloned().collect();
        for records in [&[][..], base.as_slice(), repeated.as_slice()] {
            write_fixture(&input, records);
            compute(input.to_str().unwrap(), single.to_str().unwrap(), 10, 10, &None, 1);
            let expected = std::fs::read_to_string(&single).unwrap();
            for budget in [2, 4, 8] {
                compute(input.to_str().unwrap(), parallel.to_str().unwrap(), 10, 10, &None, budget);
                assert_eq!(expected, std::fs::read_to_string(&parallel).unwrap());
            }
            if records.len() == base.len() {
                let rows: Vec<_> = expected.lines().collect();
                assert_eq!(rows.len(), 2);
                for (line, depth) in rows.iter().zip([10u32, 11]) {
                    let fields: Vec<_> = line.split('\t').collect();
                    assert_eq!(fields.len(), 22);
                    assert_eq!(fields[6..].iter().map(|x| x.parse::<u32>().unwrap()).sum::<u32>(), depth);
                }
            }
        }
    }

    #[test]
    fn test_read_and_worker_failure_shutdown() {
        for fail_read in [true, false] {
            let (tx, rx) = mpsc::channel();
            std::thread::spawn(move || {
                let records = (0..50).map(|index| {
                    if fail_read && index == 9 {
                        return Err("injected read failure".to_string());
                    }
                    let mut r = record(0, 9, index);
                    if !fail_read && index == 3 {
                        // Missing XM must still fail even below the MAPQ cutoff.
                        r.remove_aux(b"XM").unwrap();
                    }
                    Ok(r)
                });
                let result = count_batches(records, 10, &None, 3, 2);
                tx.send(result.is_err()).unwrap();
            });
            assert!(rx.recv_timeout(std::time::Duration::from_secs(5)).unwrap());
        }
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_final_text_flush_failure_is_reported() {
        let dir = TestDir::new();
        let input = dir.0.join("input.bam");
        write_fixture(&input, &depth_records());
        let result = std::panic::catch_unwind(|| {
            compute(input.to_str().unwrap(), "/dev/full", 10, 10, &None, 4);
        });
        assert!(result.is_err());
    }
}
