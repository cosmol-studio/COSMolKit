# COSMolKit for JavaScript and TypeScript

WebAssembly bindings to the public Rust API, with generated TypeScript declarations.

```bash
npm install @cosmol-studio/cosmolkit
```

```typescript
import init, { molFromSmiles } from "@cosmol-studio/cosmolkit";

await init();
const mol = molFromSmiles("CCO");
console.log(mol.toSmiles());
console.log(mol.numAtoms);
mol.free();
```

Use the full package above or select a smaller
[feature package](https://kit.cosmol.org/javascript/).
Native objects retain Rust ownership and copy-on-write behavior.
Params accept typed options or parameter instances; invalid inputs throw typed errors.

[API and examples](https://kit.cosmol.org/javascript/api) ·
[Build and development](https://github.com/cosmol-studio/COSMolKit/blob/main/wasm/dev.md) · [Binding contract](https://github.com/cosmol-studio/COSMolKit/blob/main/wasm/BINDING_CONTRACT.md)
