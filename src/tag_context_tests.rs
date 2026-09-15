// Included in tag.rs's test module to share its temporary-directory helpers.

fn context_record(reference: &str, read: &str, cigar: Vec<Cigar>, start: i64, flags: u16) -> Record {
    let mut r = Record::new();
    r.set(
        b"context",
        Some(&bam::record::CigarString(cigar)),
        read.as_bytes(),
        &vec![30; read.len()],
    );
    r.set_tid(0);
    r.set_pos(start);
    r.set_flags(flags);
    r.set_mapq(60);
    let refgenome = HashMap::from([(0usize, reference.as_bytes().to_vec())]);
    let sizes = HashMap::from([(0usize, reference.len())]);
    let xm = determine_xm_tag_string(&r, &refgenome, &sizes, &get_rcmapping(), r.is_paired());
    assert_eq!(xm.len(), read.len(), "XM must have one character per read base");
    r.push_aux(b"XM", Aux::String(&xm)).unwrap();
    r
}

#[test]
fn test_reference_context_across_indels_and_strands() {
    use Cigar::*;
    // Explicit expected tags describe genomic context, independent of CIGAR.
    let cases = vec![
        ("AACAGTAA", "CGT", vec![Match(1), Del(1), Match(2)], 2, "X.."),
        ("AACAGTAA", "TGT", vec![Match(1), Del(1), Match(2)], 2, "x.."),
        ("AACGATAA", "CAT", vec![Match(1), Del(1), Match(2)], 2, "Z.."),
        ("AACGATAA", "TAT", vec![Match(1), Del(1), Match(2)], 2, "z.."),
        ("AACGATAA", "CAGAT", vec![Match(1), Ins(1), Match(3)], 2, "Z...."),
        ("AACGATAA", "TAGAT", vec![Match(1), Ins(1), Match(3)], 2, "z...."),
        // Both inserted Cs, including an apparent inserted CpG, stay uncalled.
        ("AACGATAA", "CCCGAT", vec![Match(1), Ins(2), Match(3)], 2, "Z....."),
        ("AACAGTAA", "CGT", vec![Match(1), RefSkip(1), Match(2)], 2, "X.."),
        ("AACGAA", "C", vec![Match(1), Del(1)], 2, "Z"),
        ("AACAGAA", "C", vec![Match(1), Del(2)], 2, "X"),
        ("AACGATAA", "CGAT", vec![Match(4)], 2, "Z..."),
        ("AACGATAA", "TCGATT", vec![SoftClip(1), Match(4), SoftClip(1)], 2, ".Z...."),
        ("AACGATAA", "CGAT", vec![Ins(1), Del(1), Match(3)], 2, "...."),
        ("AACGATAA", "CTAT", vec![HardClip(2), Equal(1), Diff(1), Pad(1), Match(2), HardClip(1)], 2, "Z..."),
        ("AACGACAGACATAA", "CGACAGACAT", vec![Match(10)], 2, "Z..X...H.."),
        ("AACGACAGACATAA", "TGATAGATAT", vec![Match(10)], 2, "z..x...h.."),
        ("AACGAA", "AG", vec![Match(2)], 2, ".."),
        ("AACGAA", "NG", vec![Match(2)], 2, ".."),
        ("AACRAA", "CR", vec![Match(2)], 2, ".."),
        // Context outside the read still comes from the genome.
        ("AACGAA", "C", vec![Match(1)], 2, "Z"),
        ("CG", "CG", vec![Match(2)], 0, "Z."),
        ("C", "C", vec![Match(1)], 0, "U"),
        ("C", "T", vec![Match(1)], 0, "u"),
    ];
    let rc = get_rcmapping();
    for (reference, read, cigar, start, expected) in cases {
        let reference_span: i64 = cigar.iter().map(|op| match *op {
            Match(n) | Equal(n) | Diff(n) | Del(n) | RefSkip(n) => n as i64,
            _ => 0,
        }).sum();
        // SE top, PE top read 1, PE top read 2.
        for flags in [0, 99, 147] {
            let r = context_record(reference, read, cigar.clone(), start, flags);
            assert_eq!(r.aux(b"XM").unwrap(), Aux::String(expected),
                "reference={reference}, read={read}, flags={flags}");
        }
        // Mirror the complete fixture for SE bottom and both PE bottom mates.
        let bottom_reference = reverse_complement(reference, &rc);
        let bottom_read = reverse_complement(read, &rc);
        let bottom_cigar: Vec<Cigar> = cigar.iter().rev().copied().collect();
        let bottom_expected: String = expected.chars().rev().collect();
        let bottom_start = reference.len() as i64 - start - reference_span;
        for flags in [16, 83, 163] {
            let r = context_record(
                &bottom_reference, &bottom_read, bottom_cigar.clone(), bottom_start, flags,
            );
            assert_eq!(r.aux(b"XM").unwrap(), Aux::String(&bottom_expected),
                "mirrored reference={reference}, read={read}, flags={flags}");
        }
    }
}

#[test]
fn test_reference_context_preserves_cpg_coordinates_on_both_strands() {
    for flags in [0, 16, 99, 147, 83, 163] {
        let r = context_record("AACGCGCGCGAA", "CGCGCGCG", vec![Cigar::Match(8)], 2, flags);
        let br = crate::readutil::BismarkRead::new(&r);
        let positions: Vec<i32> = br.get_cpg_positions().iter().map(|p| p.pos).collect();
        assert_eq!(positions, vec![2, 4, 6, 8], "flags={flags}");
        let (_, patterns) = br.get_cpg_quartets_and_patterns();
        assert_eq!(patterns, vec![15]);
    }
}

#[test]
fn test_reference_context_deletion_does_not_create_a_quartet() {
    let r = context_record(
        "AACAGTCGCGCGAA", "CGTCGCGCG",
        vec![Cigar::Match(1), Cigar::Del(1), Cigar::Match(8)], 2, 99,
    );
    let br = crate::readutil::BismarkRead::new(&r);
    assert_eq!(br.get_num_cpgs(), 3);
    let (quartets, patterns) = br.get_cpg_quartets_and_patterns();
    assert!(quartets.is_empty());
    assert!(patterns.is_empty());
}

#[test]
fn test_reference_context_tag_to_entropy_single_and_multi_thread() {
    let dir = TestDir::new();
    let reference = "AACAGTCGCGCGCGAA";
    let genome = dir.0.join("ref.fa");
    std::fs::write(&genome, format!(">chr1\n{reference}\n")).unwrap();
    std::fs::write(dir.0.join("ref.fa.fai"), format!(
        "chr1\t{}\t6\t{}\t{}\n", reference.len(), reference.len(), reference.len() + 1,
    )).unwrap();

    let input = dir.0.join("input.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1").push_tag(b"LN", reference.len()),
    );
    let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
    for i in 0..10 {
        let read = if i < 5 { "CGTCGCGCGCG" } else { "TGTTGTGTGTG" };
        let mut r = context_record(
            reference, read, vec![Cigar::Match(1), Cigar::Del(1), Cigar::Match(10)], 2, 99,
        );
        // Exercise the real tag command below, starting with untagged records.
        r.remove_aux(b"XM").unwrap();
        writer.write(&r).unwrap();
    }
    drop(writer);

    let mut baseline = None;
    let mut baseline_bam = None;
    for tag_threads in [1, 8] {
        let tagged = dir.0.join(format!("tagged-{tag_threads}.bam"));
        run(input.to_str().unwrap(), tagged.to_str().unwrap(), genome.to_str().unwrap(), tag_threads);
        let records = decoded(&tagged);
        assert_eq!(records.1.len(), 10);
        if let Some(expected) = &baseline_bam {
            assert_eq!(&records, expected);
        } else {
            baseline_bam = Some(records);
        }
        for me_threads in [1, 8] {
            let output = dir.0.join(format!("entropy-{tag_threads}-{me_threads}.txt"));
            crate::me::compute(tagged.to_str().unwrap(), output.to_str().unwrap(), 10, 10, &None, me_threads);
            let text = std::fs::read_to_string(output).unwrap();
            let rows: Vec<&str> = text.lines().collect();
            assert_eq!(rows.len(), 1, "Only the four genuine reference CpGs form a quartet");
            let fields: Vec<&str> = rows[0].split('\t').collect();
            assert_eq!(fields.len(), 22);
            assert_eq!(fields[0], "chr1");
            assert_eq!(&fields[1..5], &["6", "8", "10", "12"]);
            assert_eq!(fields[5].parse::<f64>().unwrap(), 0.25);
            assert_eq!(fields[6], "5");  // c0: all unmethylated
            assert!(fields[7..21].iter().all(|count| *count == "0"));
            assert_eq!(fields[21], "5"); // c15: all methylated
            if let Some(expected) = &baseline {
                assert_eq!(&text, expected);
            } else {
                baseline = Some(text);
            }
        }
    }
}


// These expectations are fixed by the reference sequence and Shannon entropy,
// not by V1/V2 or by comparing two paths through the same algorithm.
fn truth_record(start: i64, pattern: usize, top: bool, flags: u16, xg: Option<&str>) -> Record {
    let mut bases = Vec::new();
    for bit in [8, 4, 2, 1] {
        let methylated = pattern & bit != 0;
        if top {
            bases.extend_from_slice(&[if methylated { b'C' } else { b'T' }, b'G']);
        } else {
            bases.extend_from_slice(&[b'C', if methylated { b'G' } else { b'A' }]);
        }
    }
    let mut r = Record::new();
    r.set(b"truth", Some(&bam::record::CigarString(vec![Cigar::Match(8)])), &bases, &[30; 8]);
    r.set_tid(0);
    r.set_pos(start);
    r.set_flags(flags);
    r.set_mapq(60);
    if let Some(xg) = xg {
        r.push_aux(b"XG", Aux::String(xg)).unwrap();
    }
    r
}

#[test]
fn test_truth_all_patterns_strands_flags_and_xg() {
    let reference = HashMap::from([(0usize, b"AACGCGCGCGAA".to_vec())]);
    let sizes = HashMap::from([(0usize, 12usize)]);
    // Include extra flag bits and XG values that override the flag fallback.
    let orientations = [
        (true, 0, None), (false, 16, None),
        (true, 99, None), (true, 147, None),
        (false, 83, None), (false, 163, None),
        (true, 1123, None), (false, 1187, None),
        (true, 83, Some("CT")), (false, 99, Some("GA")),
    ];
    for (top, flags, xg) in orientations {
        for pattern in 0..16 {
            let mut r = truth_record(2, pattern, top, flags, xg);
            // A wrong file-wide pairing hint must not affect per-record handling.
            let xm = determine_xm_tag_string(&r, &reference, &sizes, &get_rcmapping(), !r.is_paired());
            let mut expected_xm = String::new();
            for bit in [8, 4, 2, 1] {
                let call = if pattern & bit != 0 { 'Z' } else { 'z' };
                if top { expected_xm.push(call); expected_xm.push('.'); }
                else { expected_xm.push('.'); expected_xm.push(call); }
            }
            assert_eq!(xm, expected_xm, "pattern={pattern}, flags={flags}, XG={xg:?}");
            r.push_aux(b"XM", Aux::String(&xm)).unwrap();
            let br = crate::readutil::BismarkRead::new(&r);
            assert_eq!(br.get_cpg_positions().iter().map(|p| p.pos).collect::<Vec<_>>(), vec![2, 4, 6, 8]);
            assert_eq!(br.get_cpg_quartets_and_patterns().1, vec![pattern]);
        }
    }
}

#[test]
fn test_truth_rejects_malformed_xm_and_xg() {
    for xm in ["Z.Z.Z.", "Z.Z.Z.Z...", "Z.Z.Z.é"] {
        let mut r = truth_record(2, 15, true, 0, None);
        r.push_aux(b"XM", Aux::String(xm)).unwrap();
        assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            crate::readutil::BismarkRead::new(&r);
        })).is_err(), "Malformed XM was silently accepted: {xm}");
    }
    let mut r = truth_record(2, 15, true, 0, Some("invalid"));
    r.push_aux(b"XM", Aux::String("Z.Z.Z.Z.")).unwrap();
    assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        crate::readutil::BismarkRead::new(&r);
    })).is_err());
}

fn assert_truth_entropy(text: &str, expected: &[(i64, [u32; 16], f64)]) {
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows.len(), expected.len());
    for (row, (start, counts, entropy)) in rows.iter().zip(expected) {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields.len(), 22);
        assert_eq!(fields[0], "chr1");
        let positions: Vec<i64> = fields[1..5].iter().map(|x| x.parse().unwrap()).collect();
        assert_eq!(positions, vec![*start, start + 2, start + 4, start + 6]);
        let actual_counts: Vec<u32> = fields[6..].iter().map(|x| x.parse().unwrap()).collect();
        assert_eq!(actual_counts.as_slice(), counts.as_slice());
        let actual_entropy: f64 = fields[5].parse().unwrap();
        assert!(actual_entropy.is_finite());
        assert!((actual_entropy - entropy).abs() < 1e-6,
            "start={start}: expected {entropy}, got {actual_entropy}");
    }
}

#[test]
fn test_truth_bam_tag_me_known_entropies_and_batch_boundary() {
    let dir = TestDir::new();
    let genome = dir.0.join("truth.fa");
    let reference = format!("AACGCGCGCG{}", "A".repeat(10)).repeat(5);
    std::fs::write(&genome, format!(">chr1\n{reference}\n")).unwrap();
    std::fs::write(dir.0.join("truth.fa.fai"), format!(
        "chr1\t{}\t6\t{}\t{}\n", reference.len(), reference.len(), reference.len() + 1,
    )).unwrap();
    let input = dir.0.join("truth.bam");
    let mut header = bam::Header::new();
    header.push_record(bam::header::HeaderRecord::new(b"SQ")
        .push_tag(b"SN", "chr1").push_tag(b"LN", reference.len()));
    let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();

    let mut one = [0u32; 16]; one[15] = 10;
    let mut two = [0u32; 16]; two[0] = 5; two[15] = 5;
    let mut four = [0u32; 16]; for p in [0, 5, 10, 15] { four[p] = 3; }
    let mut biased = [0u32; 16]; biased[0] = 9; biased[15] = 3;
    // 16 * 513 = 8208 reads crosses the real 8192-record dispatch boundary.
    let all = [513u32; 16];
    let expected = [
        (2i64, one, 0.0), (22, two, 0.25), (42, four, 0.5),
        (62, biased, 0.2028195311147832), (82, all, 1.0),
    ];

    // An unplaced unmapped first record must neither panic nor determine pairing.
    let mut unmapped = Record::new();
    unmapped.set(b"unmapped", None, b"ACGT", &[30; 4]);
    unmapped.set_flags(4);
    unmapped.set_tid(-1);
    unmapped.set_pos(-1);
    writer.write(&unmapped).unwrap();
    let orientations = [
        (true, 0, None), (false, 16, None),
        (true, 99, None), (true, 147, None),
        (false, 83, None), (false, 163, None),
        (true, 83, Some("CT")), (false, 99, Some("GA")),
    ];
    let mut index = 0usize;
    for (start, counts, _) in &expected {
        for (pattern, count) in counts.iter().enumerate() {
            for _ in 0..*count {
                let (top, flags, xg) = orientations[index % orientations.len()];
                let mut r = truth_record(*start, pattern, top, flags, xg);
                // Re-tagging must replace stale calls, not fail or retain them.
                if index % 3 == 0 { r.push_aux(b"XM", Aux::String("........")).unwrap(); }
                writer.write(&r).unwrap();
                index += 1;
            }
        }
        let mut low_mapq = truth_record(*start, 7, true, 99, None);
        low_mapq.set_mapq(9);
        writer.write(&low_mapq).unwrap();
    }
    // Mapped coordinates on an unmapped record must not create methylation calls.
    let mut placed_unmapped = truth_record(2, 0, true, 4, None);
    placed_unmapped.push_aux(b"XM", Aux::String("z.z.z.z.")).unwrap();
    writer.write(&placed_unmapped).unwrap();
    drop(writer);

    for tag_threads in [1, 8] {
        let tagged = dir.0.join(format!("truth-tag-{tag_threads}.bam"));
        run(input.to_str().unwrap(), tagged.to_str().unwrap(), genome.to_str().unwrap(), tag_threads);
        let (_, records) = decoded(&tagged);
        assert_eq!(records.len(), index + expected.len() + 2);
        assert_eq!(records[0].aux(b"XM").unwrap(), Aux::String("...."));
        assert_eq!(records.last().unwrap().aux(b"XM").unwrap(), Aux::String("........"));
        for me_threads in [1, 8] {
            let output = dir.0.join(format!("truth-me-{tag_threads}-{me_threads}.txt"));
            crate::me::compute(tagged.to_str().unwrap(), output.to_str().unwrap(), 10, 10, &None, me_threads);
            assert_truth_entropy(&std::fs::read_to_string(&output).unwrap(), &expected);
        }
    }

    // me also accepts unmapped input records without requiring an XM tag.
    let mapped_input = dir.0.join("truth-me-input.bam");
    let mut writer = bam::Writer::from_path(&mapped_input, &header, bam::Format::Bam).unwrap();
    writer.write(&unmapped).unwrap();
    for _ in 0..10 {
        let mut r = truth_record(2, 15, true, 99, None);
        r.push_aux(b"XM", Aux::String("Z.Z.Z.Z.")).unwrap();
        writer.write(&r).unwrap();
    }
    writer.write(&placed_unmapped).unwrap();
    drop(writer);
    for threads in [1, 8] {
        let output = dir.0.join(format!("unmapped-me-{threads}.txt"));
        crate::me::compute(mapped_input.to_str().unwrap(), output.to_str().unwrap(), 10, 10, &None, threads);
        assert_truth_entropy(&std::fs::read_to_string(output).unwrap(), &[(2, one, 0.0)]);
    }
}
