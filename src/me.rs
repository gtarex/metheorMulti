use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
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

/// A light-weight record sent from the main thread to worker threads.
/// Contains everything a worker needs to process a single read.
struct ReadPayload {
    cpgs: Vec<readutil::CpG>,
}

impl ReadPayload {
    fn get_cpg_quartets_and_patterns(&self) -> (Vec<readutil::Quartet>, Vec<readutil::QuartetPattern>) {
        let mut quartets: Vec<readutil::Quartet> = Vec::new();
        let mut patterns: Vec<readutil::QuartetPattern> = Vec::new();

        if self.cpgs.len() < 4 {
            return (quartets, patterns);
        }

        for i in 0..self.cpgs.len() - 3 {
            let q = readutil::Quartet {
                pos1: self.cpgs[i].abspos,
                pos2: self.cpgs[i + 1].abspos,
                pos3: self.cpgs[i + 2].abspos,
                pos4: self.cpgs[i + 3].abspos,
            };
            let mut p = 0;

            if self.cpgs[i].methylated {
                p += 8;
            }
            if self.cpgs[i + 1].methylated {
                p += 4;
            }
            if self.cpgs[i + 2].methylated {
                p += 2;
            }
            if self.cpgs[i + 3].methylated {
                p += 1;
            }

            quartets.push(q);
            patterns.push(p);
        }

        (quartets, patterns)
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
    let reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let result = if threads <= 1 {
        compute_helper(input, min_qual, cpg_set)
    } else {
        compute_helper_mt(input, min_qual, cpg_set, threads)
    };

    let mut out = fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(output)
        .unwrap();

    let mut to_write: Vec<&QuartetStat> = result
        .values()
        .filter(|stat| stat.get_read_depth() >= min_depth)
        .collect();
    to_write.sort_by_key(|s| (s.pos1.tid, s.pos1.pos, s.pos2.pos, s.pos3.pos, s.pos4.pos));

    for stat in to_write {
        writeln!(out, "{}", stat.to_bedgraph_field(&header))
            .expect("Error writing to output file.");
    }
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

/// Multi-threaded implementation using channel-based worker pool.
/// Produces exactly the same results as `compute_helper`.
pub fn compute_helper_mt(
    input: &str,
    min_qual: u8,
    cpg_set: &Option<String>,
    threads: usize,
) -> HashMap<readutil::Quartet, QuartetStat> {
    let mut reader = bamutil::get_reader(input);
    let header = bamutil::get_header(&reader);

    let target_cpgs = readutil::get_target_cpgs(cpg_set, &header);

    // Create one SPSC channel per worker (Receiver doesn't implement Clone).
    let mut txs: Vec<mpsc::SyncSender<Option<ReadPayload>>> = Vec::with_capacity(threads);

    // Spawn worker threads, each owning its own local HashMap and dedicated channel.
    let workers: Vec<thread::JoinHandle<HashMap<readutil::Quartet, QuartetStat>>> = (0..threads)
        .map(|_| {
            let (tx, rx) = mpsc::sync_channel::<Option<ReadPayload>>(4);
            txs.push(tx);
            thread::spawn(move || {
                let mut local_map: HashMap<readutil::Quartet, QuartetStat> = HashMap::new();
                while let Some(payload) = rx.recv().unwrap() {
                    let (quartets, patterns) = payload.get_cpg_quartets_and_patterns();
                    for (q, p) in quartets.iter().zip(patterns.iter()) {
                        let stat = local_map.entry(*q).or_insert(QuartetStat::new(*q));
                        stat.add_quartet_pattern(*p);
                    }
                }
                local_map
            })
        })
        .collect();

    // Main thread: read BAM, filter, round-robin dispatch to workers.
    let mut readcount = 0;
    let mut valid_readcount = 0;
    let mut worker_idx: usize = 0;

    let bar = progressbar::ProgressBar::new();

    for r in reader.records().map(|r| r.unwrap()) {
        let mut br = readutil::BismarkRead::new(&r);

        if let Some(ref target_cpgs) = target_cpgs {
            br.filter_isin(target_cpgs);
        }

        readcount += 1;

        if r.mapq() < min_qual {
            continue;
        }
        valid_readcount += 1;

        let payload = ReadPayload {
            cpgs: br.get_cpgs().clone(),
        };
        txs[worker_idx].send(Some(payload)).expect("Error sending to worker thread");
        worker_idx = (worker_idx + 1) % threads;

        if readcount % 10000 == 0 {
            bar.update(readcount, valid_readcount)
        };
    }

    // Signal all workers to finish.
    for tx in &txs {
        tx.send(None).expect("Error sending termination signal");
    }

    // Merge all local maps.
    let mut final_map: HashMap<readutil::Quartet, QuartetStat> = HashMap::new();
    for worker in workers {
        let local_map = worker.join().expect("Worker thread panicked");
        final_map = merge_maps(final_map, local_map);
    }

    final_map
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

    #[test]
    fn test_mt_matches_st() {
        // Verify that multi-threaded produces the same results as single-threaded.
        for &input in &["tests/test1.bam", "tests/test2.bam", "tests/test3.bam", "tests/test4.bam", "tests/test5.bam"] {
            let min_qual = 10u8;
            let cpg_set = None;

            let st = compute_helper(input, min_qual, &cpg_set);
            let mt = compute_helper_mt(input, min_qual, &cpg_set, 4);

            assert_eq!(st.len(), mt.len());
            for (q, st_stat) in st.iter() {
                let mt_stat = mt.get(q).expect("Quartet missing in MT result");
                assert_eq!(st_stat.get_read_depth(), mt_stat.get_read_depth(),
                    "Read depth mismatch for quartet {:?}", q);
                assert!((st_stat.compute_me() - mt_stat.compute_me()).abs() < 1e-6,
                    "ME mismatch for quartet {:?}: ST={} MT={}", q, st_stat.compute_me(), mt_stat.compute_me());
                for i in 0..16 {
                    assert_eq!(st_stat.quartet_pattern_counts[i], mt_stat.quartet_pattern_counts[i],
                        "Pattern count mismatch at index {} for quartet {:?}", i, q);
                }
            }
        }
    }
}