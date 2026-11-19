use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

pub fn create_temp_beside(path: &Path) -> io::Result<(File, PathBuf)> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.subsec_nanos());
    let mut last = None;
    for attempt in 0..16u32 {
        let tmp = path.with_file_name(format!(
            ".{name}.tmp-{}-{nanos:09}-{attempt}",
            std::process::id()
        ));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&tmp)
        {
            Ok(f) => return Ok((f, tmp)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => last = Some(e),
            Err(e) => return Err(e),
        }
    }
    Err(last.unwrap_or_else(|| std::io::Error::other("no free temporary name")))
}

pub fn write_atomic(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let (mut file, tmp) = create_temp_beside(path)?;
    let write = || -> io::Result<()> {
        file.write_all(bytes)?;
        file.sync_all()?;
        drop(file);
        std::fs::rename(&tmp, path)
    };
    let result = write();
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::create_temp_beside;

    #[test]
    fn never_reuses_an_existing_name() {
        let dir = std::env::temp_dir().join(format!("fayde-fsutil-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let target = dir.join("out.wcs");
        let (_f, a) = create_temp_beside(&target).unwrap();
        let (_g, b) = create_temp_beside(&target).unwrap();
        assert_ne!(a, b);
        assert!(a.exists() && b.exists());
        assert_eq!(a.parent(), target.parent());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
