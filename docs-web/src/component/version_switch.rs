use dioxus::prelude::*;

#[component]
pub(crate) fn VersionSwitch() -> Element {
    // The deployment publishes the same catalog. These fallback links remain
    // usable if JavaScript or the latest site's catalog is unavailable.
    let catalog: serde_json::Value =
        serde_json::from_str(include_str!("../../versions.json")).expect("docs version catalog");
    let versions = catalog["versions"].as_array().expect("docs versions");
    let loader = include_str!("../../assets/version-switch.js");
    rsx! {
        details { id: "docs-version-switch", class: "docs-version-switch",
            summary { aria_label: "Documentation version", span { id: "docs-version-current", "latest" } }
            nav { id: "docs-version-options", aria_label: "Documentation versions",
                for entry in versions {
                    a {
                        href: entry["url"].as_str().expect("version URL"),
                        "data-docs-version": entry["version"].as_str().expect("version name"),
                        {entry["version"].as_str().expect("version name")}
                    }
                }
            }
        }
        document::Script {
            id: "docs-version-loader", r#type: "module",
            "{loader}"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{cell::RefCell, rc::Rc};

    #[derive(Default)]
    struct RecordingDocument(RefCell<Vec<(Vec<(&'static str, String)>, String)>>);

    impl document::Document for RecordingDocument {
        fn eval(&self, js: String) -> document::Eval {
            document::NoOpDocument.eval(js)
        }

        fn create_script(&self, props: document::ScriptProps) {
            self.0.borrow_mut().push((
                props.attributes(),
                props.script_contents().ok().expect("one script text node"),
            ));
        }
    }

    #[test]
    fn version_switch_emits_one_complete_static_module() {
        let document = Rc::new(RecordingDocument::default());
        let mut dom = VirtualDom::new(VersionSwitch);
        dom.provide_root_context(document.clone() as Rc<dyn document::Document>);
        dom.rebuild_in_place();
        let scripts = document.0.borrow();
        assert_eq!(scripts.len(), 1);
        assert_eq!(scripts[0].1, include_str!("../../assets/version-switch.js"));
        assert!(scripts[0].0.contains(&("id", "docs-version-loader".into())));
        assert!(scripts[0].0.contains(&("type", "module".into())));
        assert!(!scripts[0].0.iter().any(|(name, _)| *name == "src"));
    }
}
