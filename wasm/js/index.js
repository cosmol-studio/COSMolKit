import init from './cosmolkit_wasm.js';

const wasmUrl = new URL('./cosmolkit_wasm_bg.wasm', import.meta.url);

if (wasmUrl.protocol === 'file:') {
    // Leave the Node-only import out of browser bundles.
    const fsModule = 'node:fs/promises';
    const { readFile } = await import(/* webpackIgnore: true */ /* @vite-ignore */ fsModule);
    await init({ module_or_path: await readFile(wasmUrl) });
} else {
    await init({ module_or_path: wasmUrl });
}

export * from './cosmolkit_wasm.js';
