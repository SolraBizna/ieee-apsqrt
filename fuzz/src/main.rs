use std::{
    ops::AddAssign,
    sync::Arc,
    time::Instant,
};

use clap::Parser;
use liso::{liso, println};
use rand::prelude::*;
use rustc_apfloat::{Round, Status};
use parking_lot::Mutex;

const ATTEMPTS_PER_REPORT: usize = 10000;

#[derive(Default, Clone)]
struct TrialResults {
    trials: u64,
    exacts: u64,
    inexact_disagreements: u64,
    inexact_error_sum: u64,
    inexact_error_max: u32,
    trials_f32: u64,
    trials_f64: u64,
    total_iterations_f32: u64,
    total_iterations_f64: u64,
    fail: u64,
}

impl TrialResults {
    fn report(&self, start_time: &Instant) -> String {
        let now = Instant::now();
        let elapsed = now - *start_time;
        format!(
            " {} trials, {} exact, {}={:.8}% disagreement\n inexact error mean={:.2} max={} (ulp); {} total fails\n {:.2}/{:.2} mean iterations per trial, {:.2}/{:.2} answers per second ",
            self.trials,
            self.exacts,
            self.inexact_disagreements,
            self.inexact_disagreements as f64 * 100.0 / self.trials as f64,
            self.inexact_error_sum as f64 / self.inexact_disagreements as f64,
            self.inexact_error_max,
            self.fail,
            self.total_iterations_f32 as f64 / self.trials_f32 as f64,
            self.total_iterations_f64 as f64 / self.trials_f64 as f64,
            self.trials_f32 as f64 / elapsed.as_secs_f64(),
            self.trials_f64 as f64 / elapsed.as_secs_f64(),
        )
    }
}

impl AddAssign<&TrialResults> for TrialResults {
    fn add_assign(&mut self, rhs: &Self) {
        self.trials += rhs.trials;
        self.exacts += rhs.exacts;
        self.inexact_disagreements += rhs.inexact_disagreements;
        self.inexact_error_sum += rhs.inexact_error_sum;
        self.inexact_error_max = self.inexact_error_max.max(rhs.inexact_error_max);
        self.trials_f32 += rhs.trials_f32;
        self.trials_f64 += rhs.trials_f64;
        self.total_iterations_f32 += rhs.total_iterations_f32;
        self.total_iterations_f64 += rhs.total_iterations_f64;
        self.fail += rhs.fail;
    }
}

macro_rules! make_attempt {
    ($name:ident, $inner:ident, $error_threshold:expr, $u:ident, $f:ident, $fstr:expr, $trials:ident, $total_iterations:ident) => {
fn $name(results: &mut TrialResults) {
    let mut rng = thread_rng();
    let bits = rng.gen::<$u>() & ($u::MAX >> 1); // don't generate negative floats
    if !$f::from_bits(bits).is_finite() {
        // we're not fuzzing NaNs / infinities
        return
    }
    let native_sqrt = (($f::from_bits(bits) as f64).sqrt() as $f).to_bits();
    //println!(format!("{:X} -> {:X}", bits, native_sqrt));
    let (evil_sqrt, iterations) = ieee_apsqrt::$inner::<$u>(bits, Round::NearestTiesToEven);
    results.trials += 1;
    results.$trials += 1;
    if evil_sqrt.status == Status::OK {
        // exact result!
        let evil_sqrt = evil_sqrt.value;
        if evil_sqrt != native_sqrt {
            println!(format!("exact result but our sqrt was wrong! sqrt({bits:08X}) should be {native_sqrt:08X}, is {evil_sqrt:08X}"));
            results.fail += 1;
        }
        results.exacts += 1;
    } else {
        if evil_sqrt.status != Status::INEXACT {
            println!(format!("WARNING: result or sqrt({bits:08X}) not merely inexact but {:?}", evil_sqrt.status));
        }
        let evil_sqrt = evil_sqrt.value;
        if evil_sqrt != native_sqrt {
            results.inexact_disagreements += 1;
            let error = (if evil_sqrt > native_sqrt { evil_sqrt - native_sqrt } else { native_sqrt - evil_sqrt }) as u32;
            results.inexact_error_sum += error as u64;
            results.inexact_error_max = results.inexact_error_max.max(error);
            if error > $error_threshold {
                println!(
                    format!(concat!($fstr, " square {:08X}/{:e}, true root = {:08X}/{:e}, bad root = {:08X}/{:e}"),
                    bits, $f::from_bits(bits),
                    native_sqrt, $f::from_bits(native_sqrt),
                    evil_sqrt, $f::from_bits(evil_sqrt))
                );
            }
        }
    }
    results.$total_iterations += iterations as u64;
}}}

make_attempt!(make_attempt_32_slower, sqrt_accurate, 0, u32, f32, "f32", trials_f32, total_iterations_f32);
make_attempt!(make_attempt_64_slower, sqrt_accurate, 0, u64, f64, "f64", trials_f64, total_iterations_f64);
make_attempt!(make_attempt_32_slow, sqrt_fast, 2, u32, f32, "f32", trials_f32, total_iterations_f32);
make_attempt!(make_attempt_64_slow, sqrt_fast, 2, u64, f64, "f64", trials_f64, total_iterations_f64);

#[derive(Debug, Parser)]
#[clap(about, author)]
struct Invocation {
    /// Number of threads that will be fuzzing f32s. Default is number of CPUs
    /// divided by two (rounding up).
    #[clap(short, long)]
    float_threads: Option<usize>,
    /// Number of threads that will be fuzzing f64s. Default is number of CPUs
    /// divided by two (rounding up).
    #[clap(short, long)]
    double_threads: Option<usize>,
    /// Use the faster, less accurate sqrt.
    #[clap(short, long)]
    inexact: bool,
}

fn worker(f: impl Fn(&mut TrialResults), results_out: Arc<Mutex<(TrialResults,bool)>>) {
    let mut results = TrialResults::default();
    loop {
        for _ in 0 .. ATTEMPTS_PER_REPORT {
            f(&mut results);
        }
        let mut lock = results_out.lock();
        lock.0 = results.clone();
        if lock.1 { break }
    }
}

fn set_liso_panic_hook() {
    std::panic::set_hook(Box::new(|info| {
        let reason = if let Some(s) = info.payload().downcast_ref::<String>() {
            s.clone()
        } else if let Some(s) = info.payload().downcast_ref::<&str>() {
            s.to_string()
        } else { "(unknown payload)".to_string() };
        let location = match info.location() {
            None => "(unknown location)".to_string(),
            Some(loc) => format!("{}:{}", loc.file(), loc.line()),
        };
        println!(fg=red, bold, "PANIC! ", location, "\n", -bold, reason);
    }));
}

fn main() {
    let invocation = Invocation::parse();
    let mut io = liso::InputOutput::new();
    io.prompt("", false, true);
    set_liso_panic_hook();
    let mut result_pools = vec![];
    let mut float_handles = vec![];
    let inexact = invocation.inexact;
    for n in 0 .. invocation.float_threads.unwrap_or_else(|| (num_cpus::get() + 1) / 2) {
        let response_ptr = Arc::new(Mutex::new((TrialResults::default(),false)));
        result_pools.push(response_ptr.clone());
        float_handles.push(std::thread::Builder::new()
            .name(format!("f32[{n}]"))
            .spawn(move || {
                set_liso_panic_hook();
                if inexact {
                    worker(make_attempt_32_slow, response_ptr)
                } else {
                    worker(make_attempt_32_slower, response_ptr)
                }
            }).unwrap());
    }
    let mut double_handles = vec![];
    for n in 0 .. invocation.double_threads.unwrap_or_else(|| (num_cpus::get() + 1) / 2) {
        let response_ptr = Arc::new(Mutex::new((TrialResults::default(),false)));
        result_pools.push(response_ptr.clone());
        double_handles.push(std::thread::Builder::new()
            .name(format!("f64[{n}]"))
            .spawn(move || {
                set_liso_panic_hook();
                if inexact {
                    worker(make_attempt_64_slow, response_ptr)
                } else {
                    worker(make_attempt_64_slower, response_ptr)
                }
            }).unwrap());
    }
    if float_handles.len() == 0 && double_handles.len() == 0 {
        println!("Please actually ask for at least one trial thread.");
        drop(io);
        std::process::exit(1);
    }
    let start_time = Instant::now();
    'outer: loop {
        let results = result_pools.iter().fold(TrialResults::default(),
        |mut a,l| {
            a += &l.lock().0;
            a
        });
        io.status(Some(liso!(
            bg = blue, fg = white, +bold,
            results.report(&start_time),
        )));
        while let Some(response) = io.read_timeout(std::time::Duration::from_secs(1)) {
            use liso::Response;
            match response {
                Response::Dead => return,
                Response::Quit => {
                    io.println("^C");
                    break 'outer;
                },
                _ => (),
            }
        }
    }
    for mutex in result_pools.iter() {
        mutex.lock().1 = true;
    }
    for handle in float_handles.into_iter() { let _ = handle.join(); }
    for handle in double_handles.into_iter() { let _ = handle.join(); }
    let results = result_pools.iter().fold(TrialResults::default(),
    |mut a,l| {
        a += &l.lock().0;
        a
    });
    println!(format!("Final report:\n{}", results.report(&start_time)));
}
