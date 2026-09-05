# COSMolKit WASM

`@cosmol-studio/cosmolkit` is the WebAssembly projection of the Rust `cosmolkit` facade. The
package uses the `wasm-bindgen` `web` target and exposes an opaque
`Molecule` value. Parsing and chemistry errors are thrown as JavaScript
exceptions; successful methods return strings, numbers, typed arrays, or new
`Molecule` values.

## Install

```bash
npm install @cosmol-studio/cosmolkit
```

The package is browser/bundler-oriented. It can be used from Vite, webpack,
Rollup, and other tools that support the web `WebAssembly` target.

## Native JavaScript

Initialize the module once before using `Molecule`:

```js
import init, { Molecule } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCO");
console.log(molecule.toSmiles());
console.log(molecule.numAtoms());
molecule.free();
```

The published npm package is scoped as `@cosmol-studio/cosmolkit`. The Rust
binding crate remains `cosmolkit-wasm` internally.

## Native JavaScript From A CDN

For a browser page that does not use a bundler, import the generated web
module from jsDelivr. Pin the package version in production and keep the
background `.wasm` file on the same CDN host:

```html
<div id="molecule"></div>
<script type="module">
  import init, { Molecule } from
    "https://cdn.jsdelivr.net/npm/@cosmol-studio/cosmolkit@0.3.0/wasm_wasm.js";

  await init(
    "https://cdn.jsdelivr.net/npm/@cosmol-studio/cosmolkit@0.3.0/wasm_wasm_bg.wasm",
  );
  const molecule = Molecule.fromSmiles("c1ccccc1O");
  document.querySelector("#molecule").textContent = molecule.toSmiles();
  molecule.free();
</script>
```

The generated wasm-bindgen module is currently named `wasm_wasm.js` because
Alef derives the generated Cargo package from the `wasm/` source directory.
The package export remains `Molecule`; the generated filename is an artifact
of the binding build and should be updated here if that output layout changes.

More focused examples are in [`examples/vanilla`](./examples/vanilla/):

- [`basic.js`](./examples/vanilla/basic.js): parse, write, and render SVG.
- [`cdn.html`](./examples/vanilla/cdn.html): run the web module directly from
  jsDelivr without a bundler.
- [`descriptors.js`](./examples/vanilla/descriptors.js): descriptors and typed
  coordinate arrays.
- [`fingerprints.js`](./examples/vanilla/fingerprints.js): typed fingerprint
  bit vectors.
- [`binary.js`](./examples/vanilla/binary.js): binary serialization roundtrip.
- [`substructure.js`](./examples/vanilla/substructure.js): SMARTS matching and
  disconnected fragments.
- [`stereo.js`](./examples/vanilla/stereo.js): stereochemistry and conformers.

## React

The React example keeps initialization and molecule ownership in component
state. It does not call the Rust API during render:

```jsx
import { useEffect, useState } from "react";
import init, { Molecule } from "@cosmol-studio/cosmolkit";

export function MoleculePanel({ smiles }) {
  const [ready, setReady] = useState(false);
  const [molecule, setMolecule] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    let active = true;
    init()
      .then(() => {
        if (active) setReady(true);
      })
      .catch((value) => {
        if (active) setError(String(value));
      });
    return () => {
      active = false;
    };
  }, []);

  useEffect(() => {
    if (!ready) return;
    try {
      const next = Molecule.fromSmiles(smiles);
      setMolecule(next);
      setError(null);
    } catch (value) {
      setError(String(value));
    }
  }, [ready, smiles]);

  useEffect(() => () => molecule?.free(), [molecule]);

  if (error) return <pre>{error}</pre>;
  if (!molecule) return <p>Loading COSMolKit...</p>;
  return <output>{molecule.toSmiles()} ({molecule.numAtoms()} atoms)</output>;
}
```

The complete copyable component is
[`examples/react/MoleculePanel.jsx`](./examples/react/MoleculePanel.jsx).

## Vue

Vue applications should initialize in `onMounted` and replace the previous
opaque value when the input changes:

```vue
<script setup>
import { onBeforeUnmount, onMounted, ref, watch } from "vue";
import init, { Molecule } from "@cosmol-studio/cosmolkit";

const props = defineProps({ smiles: { type: String, default: "CCO" } });
const molecule = ref(null);
const error = ref("");

async function update() {
  try {
    const next = Molecule.fromSmiles(props.smiles);
    molecule.value?.free();
    molecule.value = next;
    error.value = "";
  } catch (value) {
    error.value = String(value);
  }
}

onMounted(async () => {
  try {
    await init();
    await update();
  } catch (value) {
    error.value = String(value);
  }
});
watch(() => props.smiles, update);
onBeforeUnmount(() => molecule.value?.free());
</script>

<template>
  <pre v-if="error">{{ error }}</pre>
  <output v-else-if="molecule">{{ molecule.toSmiles() }}</output>
  <p v-else>Loading COSMolKit...</p>
</template>
```

The complete component is [`examples/vue/MoleculePanel.vue`](./examples/vue/MoleculePanel.vue).

## Ownership And Errors

`Molecule` values own Rust allocations. Keep a value while it is in use and
call `.free()` when a long-lived UI state or worker no longer needs it. Methods
such as `withHydrogens()`, `with2dCoordinates()`, and
`with3dConformer()` return new values and leave the source value unchanged.

Fingerprint methods return `Uint32Array` set-bit lists. Coordinate methods
return flattened `Float64Array` values. Binary serialization returns a
`Uint8Array`, which can be stored or transferred through a worker.

## Examples

The examples are intentionally small and use only the published package:

- [`examples/vanilla/basic.js`](./examples/vanilla/basic.js)
- [`examples/vanilla/cdn.html`](./examples/vanilla/cdn.html)
- [`examples/vanilla/descriptors.js`](./examples/vanilla/descriptors.js)
- [`examples/vanilla/fingerprints.js`](./examples/vanilla/fingerprints.js)
- [`examples/vanilla/binary.js`](./examples/vanilla/binary.js)
- [`examples/vanilla/substructure.js`](./examples/vanilla/substructure.js)
- [`examples/vanilla/stereo.js`](./examples/vanilla/stereo.js)
- [`examples/react/MoleculePanel.jsx`](./examples/react/MoleculePanel.jsx)
- [`examples/vue/MoleculePanel.vue`](./examples/vue/MoleculePanel.vue)
