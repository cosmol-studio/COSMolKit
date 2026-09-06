use dioxus::prelude::*;

#[component]
pub fn Seo(title: String, description: String, canonical: String) -> Element {
    rsx! {
        document::Title { "{title}" }
        document::Meta { name: "description", content: "{description}" }
        document::Link { rel: "canonical", href: "{canonical}" }
        document::Meta { property: "og:title", content: "{title}" }
        document::Meta { property: "og:description", content: "{description}" }
        document::Meta { property: "og:url", content: "{canonical}" }
    }
}
