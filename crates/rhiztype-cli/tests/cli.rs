use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rhiztype"))
}

fn demos() -> String {
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/examples/demos.fasta").to_string()
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("rhiztype_cli_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

#[test]
fn no_args_exit_2() {
    let o = bin().output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(!o.stderr.is_empty());
}

#[test]
fn empty_file_exit_1() {
    let p = tmp("empty.fasta");
    std::fs::write(&p, "").unwrap();
    let o = bin().arg(&p).output().unwrap();
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("no FASTA records found"));
}

#[test]
fn nan_p_thresh_exit_2() {
    let o = bin().args([&demos(), "--p-thresh", "nan"]).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
}

#[test]
fn unreadable_input_exit_2() {
    let o = bin().arg("/nonexistent/in.fasta").output().unwrap();
    assert_eq!(o.status.code(), Some(2));
}

#[test]
fn unwritable_out_exit_2() {
    let o = bin().args([&demos(), "--out", "/nonexistent/x.csv"]).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(!String::from_utf8_lossy(&o.stdout).contains("RHIZTYPE"));
}

#[test]
fn csv_quoting() {
    let p = tmp("q.fasta");
    std::fs::write(&p, ">=\"x\"\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let out = tmp("q.csv");
    let o = bin().arg(&p).arg("--out").arg(&out).arg("-q").output().unwrap();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    let csv = std::fs::read_to_string(&out).unwrap();
    let line = csv.lines().nth(1).unwrap();
    assert!(line.starts_with("\"'=\"\"x\"\"\","), "{}", line);
}

#[test]
fn broken_pipe_no_panic() {
    let out = tmp("bp.json");
    let mut c = bin()
        .arg(demos())
        .arg("--json")
        .arg(&out)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut rd = BufReader::new(c.stdout.take().unwrap());
    let mut line = String::new();
    rd.read_line(&mut line).unwrap();
    drop(rd);
    let st = c.wait().unwrap();
    assert_ne!(st.code(), Some(101));
    assert_eq!(st.code(), Some(0));
    assert!(out.exists());
}

#[test]
fn demos_five_results() {
    let out = tmp("demos.json");
    let o = bin().arg(demos()).arg("--json").arg(&out).arg("-q").output().unwrap();
    assert_eq!(o.status.code(), Some(0));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert_eq!(v.as_array().unwrap().len(), 5);
}
