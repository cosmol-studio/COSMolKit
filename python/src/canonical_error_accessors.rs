//! Descriptor publication for native exception payload accessors.
use pyo3::prelude::*;

pub(crate) fn attach(class: &Bound<'_, PyAny>, methods: &[(&str, &str)]) -> PyResult<()> {
    let definitions = PyModule::from_code(
        class.py(),
        c"def accessor(name, field):\n    def read(self):\n        return getattr(self, field)\n    read.__doc__ = 'Return the {} context recorded by this error.'.format(name.replace('_', ' '))\n    read.__name__ = name\n    read.__module__ = 'cosmolkit'\n    return read\n",
        c"_cosmolkit_error_accessors",
        c"_cosmolkit_error_accessors",
    )?;
    for (name, field) in methods {
        class.setattr(
            *name,
            definitions.getattr("accessor")?.call1((*name, *field))?,
        )?;
    }
    Ok(())
}

pub(crate) fn residue_error(class: &Bound<'_, PyAny>) -> PyResult<()> {
    let definitions = PyModule::from_code(
        class.py(),
        c"def initialize(self, input):\n    if not isinstance(input, str):\n        raise TypeError('input must be str')\n    self._input = input\n    ValueError.__init__(self, f\"unknown residue code name '{input}'\")\n",
        c"_cosmolkit_residue_error",
        c"_cosmolkit_residue_error",
    )?;
    class.setattr("__init__", definitions.getattr("initialize")?)?;
    attach(class, &[("input", "_input")])
}
