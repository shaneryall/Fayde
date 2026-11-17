//! The `fayde` binary. All of the command lives in the library, [`fayde_cli::run`].

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    ExitCode::from(fayde_cli::run(&args))
}
