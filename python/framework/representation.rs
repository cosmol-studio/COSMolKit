//! Bounded native formatting: no chemistry, Python callbacks, or object copies.
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyString, PyTuple};
use std::fmt::{self, Write};

/// Explicit adapter fields may contain Python callbacks; never call their repr.
pub(crate) fn python_value(value: &Bound<'_, PyAny>) -> PyResult<String> {
    if value.is_none() {
        return Ok("None".into());
    }
    if value.is_callable() {
        return Ok("<callback>".into());
    }
    if let Ok(text) = value.cast::<PyString>() {
        return Ok(preview(&text.to_cow()?));
    }
    if value.is_instance_of::<PyBool>()
        || value.is_instance_of::<PyInt>()
        || value.is_instance_of::<PyFloat>()
    {
        let text = value.repr()?.to_string_lossy().into_owned();
        return Ok(if text.len() > 160 {
            preview(&text)
        } else {
            text
        });
    }
    if let Ok(values) = value.cast::<PyTuple>() {
        let mut items = values
            .iter()
            .take(5)
            .map(|value| python_value(&value))
            .collect::<PyResult<Vec<_>>>()?;
        if values.len() > 5 {
            items.push(format!("… (len={})", values.len()));
        }
        return Ok(format!(
            "({}{})",
            items.join(", "),
            if values.len() == 1 { "," } else { "" }
        ));
    }
    if let Ok(values) = value.cast::<PyDict>() {
        return Ok(format!("dict(len={})", values.len()));
    }
    if matches!(
        value
            .get_type()
            .getattr("__module__")?
            .extract::<String>()?
            .split('.')
            .next(),
        Some("cosmolkit" | "pathlib")
    ) {
        return Ok(value.repr()?.to_string_lossy().into_owned());
    }
    Ok(format!("<{}>", value.get_type().name()?))
}

pub(crate) fn configuration(
    value: &Bound<'_, PyAny>,
    name: &str,
    fields: &[&str],
) -> PyResult<String> {
    // Stateful adapters may rename native callback fields. The existing registry
    // already owns those public spellings; do not maintain another alias list.
    let registered = cosmolkit::BINDING_CONTRACT
        .iter()
        .find(|entry| entry.python_name == name && entry.python_configuration.is_some())
        .and_then(|entry| entry.python_configuration);
    let public_fields: Vec<&str> = registered.map_or_else(
        || fields.to_vec(),
        |fields| fields.iter().map(|field| field.name).collect(),
    );
    let fields = public_fields
        .iter()
        .map(|field| {
            Ok(format!(
                "{field}={}",
                python_value(&value.getattr(*field)?)?
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(format!("{name}({})", fields.join(", ")))
}

struct Preview {
    text: String,
    limit: usize,
    truncated: bool,
}
impl Write for Preview {
    fn write_str(&mut self, value: &str) -> fmt::Result {
        let available = self.limit.saturating_sub(self.text.len());
        if value.len() <= available {
            self.text.push_str(value);
            return Ok(());
        }
        let mut end = available.min(value.len());
        while !value.is_char_boundary(end) {
            end -= 1;
        }
        self.text.push_str(&value[..end]);
        self.truncated = true;
        Err(fmt::Error)
    }
}

pub(crate) fn preview(value: &impl fmt::Debug) -> String {
    let mut output = Preview {
        text: String::new(),
        limit: 160,
        truncated: false,
    };
    let _ = write!(output, "{value:?}");
    if output.truncated {
        output.text.push('…');
    }
    output.text
}

pub(crate) fn native(name: &str, value: &impl fmt::Debug) -> String {
    format!("{name}({})", preview(value))
}

/// Opt-in on the binding declaration, not a central business-type inventory.
macro_rules! native_repr {
    ($name:ident, $field:ident) => {
        #[pyo3::pymethods]
        impl $name {
            /// Compact preview of the retained native value; does not run chemistry.
            fn __repr__(&self) -> String {
                $crate::representation::native(stringify!($name), &self.$field)
            }
        }
    };
}
pub(crate) use native_repr;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounded_writer_stops_native_traversal_and_preserves_unicode() {
        struct Many(std::cell::Cell<usize>);
        impl fmt::Debug for Many {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                for _ in 0..100_000 {
                    self.0.set(self.0.get() + 1);
                    f.write_str("原子")?;
                }
                Ok(())
            }
        }
        let many = Many(std::cell::Cell::new(0));
        let text = preview(&many);
        assert!(text.ends_with('…'));
        assert!(text.len() <= 163);
        assert!(many.0.get() < 30);
    }
}
