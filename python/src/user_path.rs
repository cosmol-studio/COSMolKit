//! Source-defined path transport at the Python boundary only.
use std::ffi::OsString;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MissingHome;
impl std::fmt::Display for MissingHome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("cannot expand '~': HOME is not set")
    }
}
impl std::error::Error for MissingHome {}

pub(crate) fn expand_user_path_with_home(
    path: &str,
    home: impl FnOnce() -> Option<OsString>,
) -> Result<PathBuf, MissingHome> {
    // COSMolKit❗✔️: project pin d892ec3507c5b568c5ed5d86ae44e466f7d03855,
    // python/src/lib.rs:1075-1087, SHA
    // bdfea4d9e69db241fec7bd0615e89d64ace47d426d51d2f4023eaf533c307652.
    // fn expand_user_path(path: &str) -> PyResult<PathBuf> {
    //     if path == "~" || path.starts_with("~/") {
    //         let home = std::env::var_os("HOME")
    //             .ok_or_else(|| PyValueError::new_err("cannot expand '~': HOME is not set"))?;
    //         let mut expanded = PathBuf::from(home);
    //         if let Some(rest) = path.strip_prefix("~/") {
    //             expanded.push(rest);
    //         }
    //         Ok(expanded)
    //     } else {
    //         Ok(PathBuf::from(path))
    //     }
    // }
    // Preserve source PathBuf::push, including absolute suffix replacement,
    // empty HOME and non-Unicode HOME bytes. No account lookup or fallback.
    // Local cost remains O(path bytes), with one owned PathBuf as in the source.
    if path == "~" || path.starts_with("~/") {
        let home = home().ok_or(MissingHome)?;
        let mut expanded = PathBuf::from(home);
        if let Some(rest) = path.strip_prefix("~/") {
            expanded.push(rest);
        }
        Ok(expanded)
    } else {
        Ok(PathBuf::from(path))
    }
}

pub(crate) fn expand_user_path(path: &str) -> Result<PathBuf, MissingHome> {
    expand_user_path_with_home(path, || std::env::var_os("HOME"))
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum ImagePathError<E> {
    Directory(MissingHome),
    Export(E),
    Report(MissingHome),
    ReportWrite(E),
}

pub(crate) fn with_image_user_paths<T, E>(
    directory: &str,
    report_path: Option<&str>,
    mut home: impl FnMut() -> Option<OsString>,
    export: impl FnOnce(&Path) -> Result<T, E>,
    write_report: impl FnOnce(&Path, &T) -> Result<(), E>,
) -> Result<T, ImagePathError<E>> {
    // COSMolKit❗✔️: same pinned Python source:4852-4885.
    //         let out_dir = expand_user_path(out_dir)?;
    //         let filenames = complete_batch_filenames(filenames, self.inner.len(), &image_format)?;
    //         let report = self
    //             .inner
    //             .write_images_with_options(
    //                 out_dir.as_path(),
    //                 &image_format,
    //                 width,
    //                 height,
    //                 mode,
    //                 filenames.as_deref(),
    //                 validate_n_jobs(n_jobs)?,
    //                 progress_bar,
    //             )
    //             .map_err(batch_validation_pyerr)?;
    //         if let Some(path) = report_path {
    //             write_batch_report(path, &report)?;
    //         }
    // Directory projection precedes export. Only a returned report permits
    // report-path HOME lookup and report writing. Keep both effects sequential;
    // an export error never reads report HOME or opens the report.
    // Two path allocations as in the pinned Python source, no report cloning,
    // buffering, environment mutation, fallback or chemistry implementation.
    let directory =
        expand_user_path_with_home(directory, &mut home).map_err(ImagePathError::Directory)?;
    let result = export(&directory).map_err(ImagePathError::Export)?;
    if let Some(path) = report_path {
        let path = expand_user_path_with_home(path, &mut home).map_err(ImagePathError::Report)?;
        write_report(&path, &result).map_err(ImagePathError::ReportWrite)?;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn exact_tilde_prefix_expands_with_explicit_home() {
        for (input, expected) in [
            ("~", "/controlled/home"),
            ("~/", "/controlled/home/"),
            ("~/output", "/controlled/home/output"),
            ("~/report.json", "/controlled/home/report.json"),
            ("~/dir/../report.CSV", "/controlled/home/dir/../report.CSV"),
            ("~//absolute/report", "/absolute/report"),
        ] {
            let actual =
                expand_user_path_with_home(input, || Some("/controlled/home".into())).unwrap();
            assert_eq!(
                actual.as_os_str(),
                Path::new(expected).as_os_str(),
                "{input}"
            );
        }
    }
    #[test]
    fn all_other_source_paths_preserve_original_bytes_without_home_lookup() {
        for input in [
            "",
            ".",
            "output",
            "a//b/",
            "/absolute/report",
            "~other/output",
            "~\\output",
            "a/~/b",
            "$HOME/output",
        ] {
            let actual =
                expand_user_path_with_home(input, || panic!("literal path must not query HOME"))
                    .unwrap();
            assert_eq!(actual.as_os_str(), std::ffi::OsStr::new(input));
        }
    }
    #[test]
    fn missing_home_retains_source_error_and_empty_home_is_valid() {
        for input in ["~", "~/", "~/output", "~/report.json"] {
            assert_eq!(expand_user_path_with_home(input, || None), Err(MissingHome));
        }
        assert_eq!(
            MissingHome.to_string(),
            "cannot expand '~': HOME is not set"
        );
        assert_eq!(
            expand_user_path_with_home("~", || Some(OsString::new())).unwrap(),
            PathBuf::new()
        );
        assert_eq!(
            expand_user_path_with_home("~/output", || Some(OsString::new())).unwrap(),
            PathBuf::from("output")
        );
    }
    #[test]
    fn directory_missing_home_precedes_every_export_error_and_output_open() {
        for directory in ["~", "~/output"] {
            let result = with_image_user_paths(
                directory,
                Some("~/report.json"),
                || None,
                |_| -> Result<(), &'static str> { panic!("must not export") },
                |_, _| panic!("must not write report"),
            );
            assert_eq!(result, Err(ImagePathError::Directory(MissingHome)));
        }
    }
    #[test]
    fn resolved_directory_and_report_are_passed_once_without_parameter_mutation() {
        let directory = "~/output".to_string();
        let report = "~/report.CSV".to_string();
        let calls = Cell::new(0);
        let result = with_image_user_paths(
            &directory,
            Some(&report),
            || Some("/controlled/home".into()),
            |dir| -> Result<usize, &'static str> {
                calls.set(calls.get() + 1);
                assert_eq!(dir, Path::new("/controlled/home/output"));
                Ok(3)
            },
            |path, report| {
                calls.set(calls.get() + 1);
                assert_eq!(path, Path::new("/controlled/home/report.CSV"));
                assert_eq!(*report, 3);
                Ok(())
            },
        );
        assert_eq!(result, Ok(3));
        assert_eq!(calls.get(), 2);
        assert_eq!(directory, "~/output");
        assert_eq!(report, "~/report.CSV");
    }
    #[test]
    fn missing_report_home_does_not_hide_export_errors_or_suppress_image_export() {
        for error in [
            "invalid filenames",
            "invalid jobs",
            "drawing failure",
            "image write failure",
        ] {
            let called = Cell::new(false);
            let result = with_image_user_paths(
                "output",
                Some("~/report.json"),
                || panic!("export error must not read report HOME"),
                |dir| {
                    called.set(true);
                    assert_eq!(dir, Path::new("output"));
                    Err::<(), _>(error)
                },
                |_, _| panic!("export error must not write report"),
            );
            assert!(called.get());
            assert_eq!(result, Err(ImagePathError::Export(error)));
        }
    }
    #[test]
    fn missing_report_home_is_returned_after_successful_image_export_without_report_open() {
        let images_written = Cell::new(false);
        let result = with_image_user_paths(
            "output",
            Some("~/report.json"),
            || {
                assert!(images_written.get());
                None
            },
            |_| -> Result<(), &'static str> {
                images_written.set(true);
                Ok(())
            },
            |_, _| panic!("missing HOME must not open report"),
        );
        assert!(images_written.get());
        assert_eq!(result, Err(ImagePathError::Report(MissingHome)));
    }
    #[test]
    fn no_report_and_literal_report_keep_source_inputs_without_home_access() {
        for report in [
            None,
            Some("report.json"),
            Some("~other/report.csv"),
            Some("a//b/report.json"),
        ] {
            let writes = Cell::new(0);
            let result = with_image_user_paths(
                "output",
                report,
                || panic!("no tilde prefix"),
                |dir| -> Result<(), &'static str> {
                    assert_eq!(dir, Path::new("output"));
                    Ok(())
                },
                |path, _| {
                    writes.set(writes.get() + 1);
                    assert_eq!(Some(path), report.map(Path::new));
                    Ok(())
                },
            );
            assert_eq!(result, Ok(()));
            assert_eq!(writes.get(), usize::from(report.is_some()));
        }
    }
    #[test]
    fn report_writer_failure_preserves_typed_cause_after_image_export() {
        #[derive(Debug, PartialEq, Eq)]
        struct IoCause(u32);
        let exported = Cell::new(false);
        let result = with_image_user_paths(
            "output",
            Some("~/report.json"),
            || Some("/controlled/home".into()),
            |_| -> Result<(), IoCause> {
                exported.set(true);
                Ok(())
            },
            |path, _| {
                assert!(exported.get());
                assert_eq!(path, Path::new("/controlled/home/report.json"));
                Err(IoCause(13))
            },
        );
        assert_eq!(result, Err(ImagePathError::ReportWrite(IoCause(13))));
    }
    #[test]
    fn report_home_lookup_occurs_only_after_image_export_returns() {
        let image_returned = Cell::new(false);
        let result = with_image_user_paths(
            "output",
            Some("~/report.json"),
            || {
                assert!(
                    image_returned.get(),
                    "report HOME must not be read before image export returns"
                );
                Some("/controlled/home".into())
            },
            |_| -> Result<(), &'static str> {
                image_returned.set(true);
                Ok(())
            },
            |_, _| Ok(()),
        );
        assert_eq!(result, Ok(()));
    }
    #[test]
    fn directory_and_report_read_home_at_their_own_source_phases() {
        let phase = Cell::new(0);
        let result = with_image_user_paths(
            "~/output",
            Some("~/report.json"),
            || match phase.get() {
                0 => {
                    phase.set(1);
                    Some("/directory-home".into())
                }
                2 => {
                    phase.set(3);
                    Some("/report-home".into())
                }
                _ => panic!("wrong HOME read phase"),
            },
            |path| -> Result<(), &'static str> {
                assert_eq!(phase.get(), 1);
                assert_eq!(path, Path::new("/directory-home/output"));
                phase.set(2);
                Ok(())
            },
            |path, _| {
                assert_eq!(phase.get(), 3);
                assert_eq!(path, Path::new("/report-home/report.json"));
                phase.set(4);
                Ok(())
            },
        );
        assert_eq!(result, Ok(()));
        assert_eq!(phase.get(), 4);
    }
    #[cfg(unix)]
    #[test]
    fn home_os_bytes_are_not_lossily_converted() {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let home = OsString::from_vec(b"/controlled/nonunicode-\xff".to_vec());
        let result = with_image_user_paths(
            "~/output",
            Some("~/report.json"),
            || Some(home.clone()),
            |dir| -> Result<(), &'static str> {
                assert_eq!(
                    dir.as_os_str().as_bytes(),
                    b"/controlled/nonunicode-\xff/output"
                );
                Ok(())
            },
            |path, _| {
                assert_eq!(
                    path.as_os_str().as_bytes(),
                    b"/controlled/nonunicode-\xff/report.json"
                );
                Ok(())
            },
        );
        assert_eq!(result, Ok(()));
    }
}
