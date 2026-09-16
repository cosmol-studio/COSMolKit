use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};
use wasm_bindgen::{JsCast, prelude::*};
use web_sys::{
    Document, Element, Event, EventTarget, HtmlInputElement, MutationObserver,
    MutationObserverInit, Url, Window,
};

use super::engine::SearchIndex;

const PAGE_SIZE: usize = 30;
const DEBOUNCE_MS: i32 = 120;
const INDEX: &str = include_str!(concat!(env!("OUT_DIR"), "/search_index.json"));

thread_local! {
    static INDEX_CACHE: RefCell<Option<Rc<SearchIndex>>> = const { RefCell::new(None) };
    static ACTIVE: RefCell<Option<Rc<SearchUi>>> = const { RefCell::new(None) };
}

struct Listener {
    target: EventTarget,
    name: &'static str,
    callback: Closure<dyn FnMut(Event)>,
}

impl Drop for Listener {
    fn drop(&mut self) {
        let _ = self
            .target
            .remove_event_listener_with_callback(self.name, self.callback.as_ref().unchecked_ref());
    }
}

struct SearchUi {
    window: Window,
    document: Document,
    form: Element,
    input: HtmlInputElement,
    status: Element,
    output: Element,
    index: Rc<SearchIndex>,
    results: RefCell<Vec<usize>>,
    shown: Cell<usize>,
    revision: Cell<u64>,
    listeners: RefCell<Vec<Listener>>,
    removal_observer: MutationObserver,
    _on_removal: Closure<dyn FnMut()>,
}

impl Drop for SearchUi {
    fn drop(&mut self) {
        self.removal_observer.disconnect();
    }
}

impl SearchUi {
    fn listen(
        self: &Rc<Self>,
        target: EventTarget,
        name: &'static str,
        action: impl Fn(&Rc<Self>, Event) -> Result<(), JsValue> + 'static,
    ) -> Result<(), JsValue> {
        let weak = Rc::downgrade(self);
        let callback = Closure::wrap(Box::new(move |event: Event| {
            if let Some(ui) = weak.upgrade() {
                if let Err(error) = action(&ui, event) {
                    ui.fail(error);
                }
            }
        }) as Box<dyn FnMut(Event)>);
        target.add_event_listener_with_callback(name, callback.as_ref().unchecked_ref())?;
        self.listeners.borrow_mut().push(Listener {
            target,
            name,
            callback,
        });
        Ok(())
    }

    fn fail(&self, error: JsValue) {
        self.status
            .set_text_content(Some("Search could not complete. Reload the page to retry."));
        web_sys::console::error_1(&error);
    }

    fn write_query(&self, push: bool) -> Result<(), JsValue> {
        let url = Url::new(&self.window.location().href()?)?;
        let query = self.input.value();
        let query = query.trim();
        if query.is_empty() {
            url.search_params().delete("q");
        } else {
            url.search_params().set("q", query);
        }
        url.set_hash("");
        let history = self.window.history()?;
        let state = history.state()?;
        if push {
            history.push_state_with_url(&state, "", Some(&url.href()))?;
        } else {
            history.replace_state_with_url(&state, "", Some(&url.href()))?;
        }
        self.render(query)
    }

    fn sync_query(&self) -> Result<(), JsValue> {
        self.revision.set(self.revision.get() + 1);
        let url = Url::new(&self.window.location().href()?)?;
        let query = url.search_params().get("q").unwrap_or_default();
        self.input.set_value(&query);
        self.render(query.trim())
    }

    fn render(&self, query: &str) -> Result<(), JsValue> {
        self.output.set_inner_html("");
        self.output.set_attribute("data-query", query)?;
        self.shown.set(0);
        *self.results.borrow_mut() = self.index.search(query);
        let count = self.results.borrow().len();
        let message = if query.is_empty() {
            "Enter a topic, class, or method name to search the Python documentation.".to_string()
        } else if count == 0 {
            format!("No results for “{query}”. Try Molecule, fingerprint, or from_smiles.")
        } else {
            format!("{count} results for “{query}”.")
        };
        self.status.set_text_content(Some(&message));
        let list = self.document.create_element("ul")?;
        list.set_class_name("search");
        self.output.append_child(&list)?;
        let more = self.document.create_element("button")?;
        more.set_class_name("docs-search-more");
        more.set_attribute("type", "button")?;
        self.output.append_child(&more)?;
        self.show_more()
    }

    fn show_more(&self) -> Result<(), JsValue> {
        let list = required(&self.output, "ul.search")?;
        let more = required(&self.output, ".docs-search-more")?;
        let results = self.results.borrow();
        let end = (self.shown.get() + PAGE_SIZE).min(results.len());
        for &id in &results[self.shown.get()..end] {
            let record = &self.index.records[id];
            let row = self.document.create_element("li")?;
            let link = self.document.create_element("a")?;
            link.set_attribute("href", &record.url)?;
            link.set_text_content(Some(&record.title));
            row.append_child(&link)?;
            let location = self.document.create_element("small")?;
            location.set_text_content(Some(&record.url));
            row.append_child(&location)?;
            if !record.summary.is_empty() {
                let summary = self.document.create_element("p")?;
                summary.set_text_content(Some(&record.summary));
                row.append_child(&summary)?;
            }
            list.append_child(&row)?;
        }
        self.shown.set(end);
        more.set_text_content(Some(&format!(
            "Show more ({} remaining)",
            results.len() - end
        )));
        if end == results.len() {
            more.set_attribute("hidden", "")?;
        } else {
            more.remove_attribute("hidden")?;
        }
        Ok(())
    }
}

fn required(parent: &Element, selector: &str) -> Result<Element, JsValue> {
    parent
        .query_selector(selector)?
        .ok_or_else(|| JsValue::from_str(&format!("missing search element: {selector}")))
}

/// Mount the search UI after its module loads. Repeated route visits reuse WASM
/// and the parsed index, replacing listeners belonging to the previous form.
#[wasm_bindgen]
pub fn mount_search() -> Result<(), JsValue> {
    let window = web_sys::window().ok_or_else(|| JsValue::from_str("missing window"))?;
    let document = window
        .document()
        .ok_or_else(|| JsValue::from_str("missing document"))?;
    let Some(form) = document.get_element_by_id("docs-search-form") else {
        return Ok(());
    };
    if ACTIVE.with(|active| active.borrow().as_ref().is_some_and(|ui| ui.form == form)) {
        return Ok(());
    }
    let initial = js_sys::Reflect::get(window.as_ref(), &"cosmolkitInitialSearch".into())?;
    if !initial.is_undefined() && !initial.is_null() {
        let path = js_sys::Reflect::get(&initial, &"path".into())?.as_string();
        let query = js_sys::Reflect::get(&initial, &"query".into())?.as_string();
        if path.as_deref() == Some(&window.location().pathname()?) {
            if let Some(query) = query {
                let url = Url::new(&window.location().href()?)?;
                url.search_params().set("q", &query);
                let history = window.history()?;
                history.replace_state_with_url(&history.state()?, "", Some(&url.href()))?;
            }
        }
        js_sys::Reflect::delete_property(
            &js_sys::Object::from(window.clone()),
            &"cosmolkitInitialSearch".into(),
        )?;
    }
    let index = INDEX_CACHE.with(|cache| {
        if let Some(index) = cache.borrow().as_ref() {
            return Ok(index.clone());
        }
        let index =
            Rc::new(SearchIndex::load(INDEX).map_err(|e| JsValue::from_str(&e.to_string()))?);
        *cache.borrow_mut() = Some(index.clone());
        Ok::<_, JsValue>(index)
    })?;
    let input = required(&form, "input[name=q]")?.dyn_into::<HtmlInputElement>()?;
    // Observe removal from the document, including client-router navigation.
    // The callback holds no UI reference; dropping ACTIVE releases the DOM,
    // listeners and observer while INDEX_CACHE remains available for re-entry.
    let on_removal = Closure::new(|| {
        ACTIVE.with(|active| {
            let mut active = active.borrow_mut();
            if active.as_ref().is_some_and(|ui| !ui.form.is_connected()) {
                active.take();
            }
        });
    });
    let removal_observer = MutationObserver::new(on_removal.as_ref().unchecked_ref())?;
    let ui = Rc::new(SearchUi {
        status: document
            .get_element_by_id("docs-search-status")
            .ok_or_else(|| JsValue::from_str("missing status"))?,
        output: document
            .get_element_by_id("search-results")
            .ok_or_else(|| JsValue::from_str("missing results"))?,
        window,
        document,
        form,
        input,
        index,
        results: RefCell::new(Vec::new()),
        shown: Cell::new(0),
        revision: Cell::new(0),
        listeners: RefCell::new(Vec::new()),
        removal_observer,
        _on_removal: on_removal,
    });
    ui.listen(ui.input.clone().into(), "input", |ui, _| {
        let revision = ui.revision.get() + 1;
        ui.revision.set(revision);
        let weak = Rc::downgrade(ui);
        let callback = Closure::once_into_js(move || {
            if let Some(ui) = weak.upgrade() {
                if ui.revision.get() == revision && ui.form.is_connected() {
                    if let Err(error) = ui.write_query(false) {
                        ui.fail(error);
                    }
                }
            }
        });
        ui.window
            .set_timeout_with_callback_and_timeout_and_arguments_0(
                callback.unchecked_ref(),
                DEBOUNCE_MS,
            )?;
        Ok(())
    })?;
    ui.listen(ui.form.clone().into(), "submit", |ui, event| {
        event.prevent_default();
        ui.revision.set(ui.revision.get() + 1);
        ui.write_query(true)
    })?;
    ui.listen(ui.form.clone().into(), "reset", |ui, event| {
        event.prevent_default();
        ui.revision.set(ui.revision.get() + 1);
        ui.input.set_value("");
        ui.write_query(false)?;
        ui.input.focus()
    })?;
    ui.listen(ui.output.clone().into(), "click", |ui, event| {
        if let Some(target) = event.target().and_then(|t| t.dyn_into::<Element>().ok()) {
            if target.class_name() == "docs-search-more" {
                ui.show_more()?;
            }
        }
        Ok(())
    })?;
    ui.listen(ui.window.clone().into(), "popstate", |ui, _| {
        if ui.form.is_connected() {
            ui.sync_query()?;
        }
        Ok(())
    })?;
    if ui.input.value().is_empty() {
        ui.sync_query()?;
    } else {
        ui.write_query(false)?;
    }
    ui.form.set_attribute("data-ready", "true")?;
    let options = MutationObserverInit::new();
    options.set_child_list(true);
    options.set_subtree(true);
    ui.removal_observer
        .observe_with_options(&ui.document, &options)?;
    ACTIVE.with(|active| *active.borrow_mut() = Some(ui));
    Ok(())
}
