use dioxus::prelude::*;

#[component]
pub fn Seo(title: String, description: String) -> Element {
    let route = use_route::<crate::route::Route>();
    let canonical = route.canonical();
    let robots = if !route.metadata().indexable {
        "noindex, follow"
    } else {
        "index, follow"
    };
    let image = crate::route::SOCIAL_IMAGE_URL;
    let image_alt = "COSMolKit documentation — Python guides and API reference";
    rsx! {
        document::Title { "{title}" }
        document::Meta { name: "description", content: "{description}" }
        document::Link { rel: "canonical", href: "{canonical}" }
        document::Meta { name: "robots", content: robots }
        document::Meta { property: "og:type", content: "website" }
        document::Meta { property: "og:site_name", content: "COSMolKit Documentation" }
        document::Meta { property: "og:title", content: "{title}" }
        document::Meta { property: "og:description", content: "{description}" }
        document::Meta { property: "og:url", content: "{canonical}" }
        document::Meta { property: "og:image", content: image }
        document::Meta { property: "og:image:type", content: "image/png" }
        document::Meta { property: "og:image:width", content: "1200" }
        document::Meta { property: "og:image:height", content: "630" }
        document::Meta { property: "og:image:alt", content: image_alt }
        document::Meta { name: "twitter:card", content: "summary_large_image" }
        document::Meta { name: "twitter:title", content: "{title}" }
        document::Meta { name: "twitter:description", content: "{description}" }
        document::Meta { name: "twitter:image", content: image }
        document::Meta { name: "twitter:image:alt", content: image_alt }
    }
}
