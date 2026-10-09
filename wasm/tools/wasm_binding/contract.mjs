// Registry-driven final export/descriptor checks; no chemistry API inventory.
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const camel = name => name.replace(/_([a-z])/g, (_, letter) => letter.toUpperCase());
export function checkContract(module, document) {
    const types = new Map(document.entries.filter(row => row.item === 'type').map(row => [row.semantic_id.split('.').at(-1), row.javascript_name]));
    const errors = [];
    for (const row of document.entries.filter(row => row.platform !== 'native')) {
        const ownerName = row.owner === 'molecule' ? 'Molecule' : types.get(row.semantic_id.split('.').slice(0, -1).join('.'));
        const owner = row.item === 'type' || row.owner === 'module' ? module : row.receiver ? module[ownerName]?.prototype : module[ownerName];
        if (row.javascript_name === 'new' && row.item === 'callable') {
            if (typeof module[ownerName] !== 'function') errors.push(`${row.semantic_id}: actual constructor missing`);
        } else if (!owner || !(row.javascript_name in owner)) {
            errors.push(`${row.semantic_id}: actual JavaScript export missing`);
        }
        if (row.item === 'type') {
            for (const field of row.properties || []) {
                if (!module[row.javascript_name]?.prototype || !(camel(field.name) in module[row.javascript_name].prototype)) errors.push(`${row.javascript_name}.${camel(field.name)}: actual property missing`);
            }
        }
        if (row.role !== 'parameter') continue;
        if (row.fields === null) {
            errors.push(`${row.semantic_id}: canonical configuration constructor missing`);
            continue;
        }
        const cls = module[row.javascript_name];
        if (!cls) continue;
        for (const field of row.fields) {
            const descriptor = Object.getOwnPropertyDescriptor(cls.prototype, camel(field.name));
            if (!descriptor?.get || !descriptor?.set) errors.push(`${row.javascript_name}.${camel(field.name)}: actual getter/setter missing`);
        }
    }
    return errors;
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
    const module = await import(pathToFileURL(process.argv[2]).href);
    module.initSync({ module: readFileSync(process.argv[3]) });
    const errors = checkContract(module, JSON.parse(readFileSync(process.argv[4], 'utf8')));
    if (errors.length) {
        console.error(`WASM binding contract failed: ${errors.length} violations\n${errors.join('\n')}`);
        process.exitCode = 1;
    }
}
