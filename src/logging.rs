//! Structured logging for falak via `FALAK_LOG` env var.

/// Initialize tracing with the default `info` level.
///
/// Uses `FALAK_LOG` env var to override the filter.
pub fn init() {
    init_with_level("info");
}

/// Initialize tracing with a custom default level (e.g., `"debug"`, `"trace"`).
///
/// Uses `FALAK_LOG` env var to override the filter.
pub fn init_with_level(default_level: &str) {
    use tracing_subscriber::EnvFilter;
    use tracing_subscriber::fmt;
    use tracing_subscriber::prelude::*;

    let filter =
        EnvFilter::try_from_env("FALAK_LOG").unwrap_or_else(|_| EnvFilter::new(default_level));

    let _ = tracing_subscriber::registry()
        .with(fmt::layer().with_target(true).with_thread_ids(true))
        .with(filter)
        .try_init();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_does_not_panic() {
        init();
        init();
    }
}
