
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio, exit};
use std::io::Cursor;
use std::collections::HashMap;
use polars::lazy::dsl::GetOutput;
use polars::prelude::*;
use clap::Parser;

#[derive(Parser)]
#[command(name = "readloc")]
#[command(about = "Remove region ambiguity for reads", long_about = None)]
struct Cli {
    #[arg(short, value_name = "ALIGNED_BAM", help = "Aligned BAM files")]
    a: PathBuf,
    #[arg(short, value_name = "REGION_BED_DIR", help = "BED files generated from gtf annotations")]
    r: PathBuf,
    #[arg(short, value_name = "OUTPUT", help = "Output directory")]
    o: Option<PathBuf>,
}


fn merge_range(vec: &Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    let mut merged = vec![];
    if vec.is_empty() {
        return merged;
    }
    
    let mut start = vec[0].0;
    let mut end = vec[0].1;
    
    for &(s, e) in vec.iter().skip(1) {
        if s <= end {
            end = end.max(e);
        } else {
            merged.push((start, end));
            start = s;
            end = e;
        }
    }
    
    merged.push((start, end));
    
    merged
}
// the largest chromosome chr1 size 248,956,422 is smaller than i64 max 4294,967,295
fn cigar_parser(cigar: &str, offset: i64) -> Vec<(i64, i64)> {
    let mut valid_ranges = vec![];
    let mut start: i64 = 0;
    let mut end: i64 = 0;
    let mut i: usize = 0;
    while i < cigar.len() {
        let mut j = i;
        while j < cigar.len() && cigar[j..=j].chars().next().unwrap().is_digit(10) {
            j += 1;
        }
        let num = cigar[i..j].parse::<i64>().unwrap();
        
        if cigar[j..=j] == "M".to_string() || cigar[j..=j] == "D".to_string() {
            end += num;
            valid_ranges.push((start + offset, end + offset));
            start = end;
        } else if cigar[j..=j] == "N".to_string() {
            end += num;
            start = end;
        } 
        i = j + 1;
    }
    valid_ranges = merge_range(&valid_ranges);
    valid_ranges
}

fn calc_coverage(a: i64, b: i64, c: i64, d: i64, cigar: &str) -> i64{
    let ranges: Vec<(i64, i64)> = cigar_parser(cigar, a);
    assert_eq!(b, ranges[ranges.len()-1].1);
    let mapped_vec: Vec<i64> = ranges.iter().map(|range| {
        let btm = range.0.max(c);
        let top = range.1.min(d);
        let cov: i64;
        if top > btm {
            cov = top - btm;
        } else {
            cov = 0;
        }
        cov
    }).collect();
    let coverage = mapped_vec.iter().sum();
    coverage
}

fn check_command(cmd: &str) {
    let cmd_check = Command::new(cmd)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status().unwrap();
    if !cmd_check.success() {
        println!("{} is not installed.", cmd);
        exit(1);
    }
}

fn main() -> PolarsResult<()> {
    let cli = Cli::parse();
    check_command("bedtools");

    let abs_align_file = cli.a.canonicalize()?;
    let abs_region_dir = cli.r.canonicalize()?;

    let merged: LazyFrame;
    if let Ok(entries) = fs::read_dir(abs_region_dir.clone()) {
        let mut dfs = vec![];
        for entry in entries {
            if let Ok(entry) = entry {
                let mut abs_region_path: PathBuf = abs_region_dir.clone();
                let file_name = entry.file_name().to_string_lossy().to_string();
                let extension = Path::new(&file_name).extension();
                if let Some(ext) = extension {
                    if ext == "bed" {
                        println!("\x1b[44mStart screening overlap to {}\x1b[m", file_name);
                        abs_region_path.push(file_name.clone());
                        let bamraw = Command::new("bedtools")
                            .args(&["intersect", "-s", "-a", abs_align_file.to_str().unwrap(), 
                                "-b", abs_region_path.to_str().unwrap(), "-wa", "-split", "-ubam"])
                            .stdout(Stdio::piped())
                            .stderr(Stdio::null())
                            .spawn()
                            .expect("Error: bedtools intersect failed");
                        let bamraw_out = bamraw.stdout.expect("Error: failed to open bedtools intersect stdout");

                        let bamview = Command::new("samtools")
                            .args(&["view", "-"])
                            .stdin(Stdio::from(bamraw_out))
                            .stdout(Stdio::piped())
                            .spawn()
                            .expect("Error: failed to samtools view bam result");
                    
                        let bamout = bamview.wait_with_output().expect("Error: failed to open samtools view stdout");
                        let bamreader = Cursor::new(&bamout.stdout);
                    
                        let bedraw = Command::new("bedtools")
                            .args(&["intersect", "-s", "-a", abs_align_file.to_str().unwrap(), 
                                        "-b", abs_region_path.to_str().unwrap(), "-wo", "-split", "-bed"])
                            .output()
                            .expect("Error: bedtools intersect for BED failed");
                        let bedreader = Cursor::new(&bedraw.stdout);
                    
                        let bamdf = CsvReader::new(bamreader)
                            .with_delimiter(b'\t')
                            .has_header(false)
                            .with_projection(Some(vec![0, 2, 3, 5]))
                            .finish().expect("failed to create dataframe from bam");
                        let beddf = CsvReader::new(bedreader)
                            .with_delimiter(b'\t')
                            .has_header(false)
                            .with_projection(Some(vec![0, 1, 2, 3, 13, 14, 15]))
                            .finish().expect("failed to create dataframe from bed");
                        println!("\x1b[42mFinished overlapping\x1b[m");
                        let joined_df = beddf.lazy()
                            .join(bamdf.lazy(), 
                                [col("column_4"), col("column_1")], 
                                [col("column_1"), col("column_3")], 
                                JoinType::Inner)
                            .filter(
                                (col("column_2") + lit(1)).eq(col("column_4_right")),
                            )
                            .rename(&[
                                "column_1", "column_2", 
                                "column_3", "column_4", 
                                "column_14", "column_15", 
                                "column_16", "column_6", ],
                            &[
                                "chr", "align_0", "align_1", "read", 
                                "region_0", "region_1", "region", "cigar"
                            ])
                            .select(&[
                                col("chr"), col("align_0"), col("align_1"), col("read"), 
                                col("region_0"), col("region_1"), col("region"), col("cigar"), ]
                            );
                        dfs.push(joined_df);
                    }
                }
            }
        }

        merged = concat(&dfs, false, true)?;
        println!("\nshape: {:?} (unrefined)", merged.clone().collect()?.shape());
    } else {
        println!("Failed to read directory.");
        exit(1);
    }

    let mut abs_output_file: PathBuf;
    if let Some(output) = cli.o.as_deref() {
        if ! output.is_dir() {
            println!("Output directory do not exist.");
            exit(1);
        } 
        abs_output_file = std::fs::canonicalize(output)?;
    } else {
        abs_output_file = std::fs::canonicalize(".")?;
    }
    let task_name = abs_align_file
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap();
    abs_output_file.push(format!("{}.reloc.bed", task_name));

    let uniq: LazyFrame = merged.clone().unique(Some(vec![String::from("read")]), UniqueKeepStrategy::None);
    let duplicated: LazyFrame = merged.filter(
        col("read").is_in(lit(uniq.clone().collect()?["read"].clone())).not()
    );
    let dupcov: LazyFrame = duplicated.clone()
        .groupby_stable([col("read")])
        .agg([
            as_struct(&[col("align_0"), col("align_1"), col("region_0"), col("region_1"), col("cigar")])
            .apply(|s| {
                    let ca = s.struct_()?;
                    let s_a = &ca.fields()[0];
                    let s_b = &ca.fields()[1];
                    let s_c = &ca.fields()[2];
                    let s_d = &ca.fields()[3];
                    let s_cigar = &ca.fields()[4];

                    let ca_a = s_a.i64()?;
                    let ca_b = s_b.i64()?;
                    let ca_c = s_c.i64()?;
                    let ca_d = s_d.i64()?;
                    let ca_cigar = s_cigar.utf8()?;

                    // iterate both `ChunkedArrays`
                    let out: Int64Chunked = ca_a
                        .into_iter()
                        .zip(ca_b)
                        .zip(ca_c)
                        .zip(ca_d)
                        .zip(ca_cigar)
                        .map(|
                            ((((opt_a, 
                            opt_b),
                            opt_c), 
                            opt_d),
                            opt_cigar
                        )| match ((((
                            opt_a, 
                            opt_b),
                            opt_c), 
                            opt_d),
                            opt_cigar
                        ) {
                            ((((Some(a), Some(b)), Some(c)), Some(d)), Some(cigar)) => Some(calc_coverage(a, b ,c, d, cigar)),
                            _ => None
                        })
                        .collect();
                    Ok(Some(out.into_series()))
            },
            GetOutput::from_type(DataType::Int64),
        ).arg_max().alias("cov_idx"),
        ]);

    let dedup: LazyFrame = duplicated.join(dupcov, [col("read")], [col("read")], 
        JoinType::Inner).groupby_stable([col("read")])
        .agg([
            col("*").exclude(&["cov_idx"]).take(col("cov_idx")).first()
        ]).select([
            col("chr"),
            col("align_0"),
            col("align_1"),
            col("read"),
            col("region_0"),
            col("region_1"),
            col("region"),
            col("cigar")
        ]);

    let mut chr_map: HashMap<String, u32> = HashMap::new();
    for i in 1..=22 {
        let chr = format!("chr{}", i);
        chr_map.entry(chr).or_insert(i);
    }
    chr_map.insert("chrX".to_string(), 97);
    chr_map.insert("chrY".to_string(), 98);
    chr_map.insert("chrM".to_string(), 99);

    let mut result = concat(&[
        uniq,
        dedup
    ], false, false)?
        .with_columns([
            col("chr").map(move |x: Series|{
                let y:Series = x.utf8()?.into_iter().map(|c| {
                    let v = chr_map.get(c.unwrap()).unwrap();
                    v
                }).collect();
                Ok(Some(y))
            }, GetOutput::from_type(DataType::UInt32)).alias("chr_n")
        ])
        .sort_by_exprs(
            vec![
                col("chr_n"), 
                col("align_0"), 
                col("align_1"), 
                col("region_0"), 
                col("region_1")],
            vec![false, false, false, false, false],
            false
        )
        .select([
            col("*").exclude(&["chr_n"])
        ]).collect()?;
    println!("{:?}", result);

    
    let mut outfile = std::fs::File::create(abs_output_file.clone()).unwrap();
    CsvWriter::new(&mut outfile).has_header(false).with_delimiter(b'\t').finish(&mut result)?;
    println!("Results in \x1b[33m{}\x1b[m", abs_output_file.to_string_lossy().to_string());
    
    Ok(())
}