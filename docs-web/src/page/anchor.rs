//! Client-only fragment restoration. Static HTML uses browser-native anchors.
use std::{cell::RefCell, rc::Rc};

use dioxus::prelude::*;
use wasm_bindgen::{JsCast, closure::Closure};
use web_sys::{Event, HtmlLinkElement};

use crate::route::Route;

struct PendingScroll {
    links: Vec<HtmlLinkElement>,
    loaded: Closure<dyn FnMut(Event)>,
}

impl Drop for PendingScroll {
    fn drop(&mut self) {
        for link in &self.links {
            let _ = link
                .remove_event_listener_with_callback("load", self.loaded.as_ref().unchecked_ref());
        }
    }
}

pub(super) fn use_fragment_scroll() {
    let route = use_route::<Route>();
    let pending = use_hook(|| Rc::new(RefCell::new(None::<PendingScroll>)));
    use_effect(use_reactive((&route,), move |(_route,)| {
        pending.borrow_mut().take();
        let Some(window) = web_sys::window() else {
            return;
        };
        let Some(document) = window.document() else {
            return;
        };
        let Ok(hash) = window.location().hash() else {
            return;
        };
        let Some(encoded) = hash.strip_prefix('#').filter(|id| !id.is_empty()) else {
            return;
        };
        let Ok(id) = js_sys::decode_uri_component(encoded) else {
            return;
        };
        let id = String::from(id);
        let Ok(nodes) = document.query_selector_all("link[rel='stylesheet']") else {
            return;
        };
        let links: Vec<_> = (0..nodes.length())
            .filter_map(|index| nodes.item(index)?.dyn_into::<HtmlLinkElement>().ok())
            .filter(|link| link.sheet().is_none())
            .collect();
        let waiting = links.clone();
        let scroll = move || {
            // Styles loaded by the client change article height. Wait for them
            // before positioning, and never scroll a page the reader has left.
            if window.location().hash().ok().as_deref() != Some(hash.as_str()) {
                return true;
            }
            if waiting.iter().any(|link| link.sheet().is_none()) {
                return false;
            }
            if let Some(target) = document.get_element_by_id(&id) {
                target.scroll_into_view();
            }
            true
        };
        if !scroll() {
            // A weak reference lets completion release the listeners and callback
            // without a cycle retaining the component after it unmounts.
            let pending_weak = Rc::downgrade(&pending);
            let loaded = Closure::new(move |_: Event| {
                if scroll()
                    && let Some(pending) = pending_weak.upgrade()
                {
                    pending.borrow_mut().take();
                }
            });
            for link in &links {
                let _ =
                    link.add_event_listener_with_callback("load", loaded.as_ref().unchecked_ref());
            }
            *pending.borrow_mut() = Some(PendingScroll { links, loaded });
        }
    }));
}
