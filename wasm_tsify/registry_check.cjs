// Structural registry inspection, not a second business-behavior test suite.
const fs = require('node:fs');
const path = require('node:path');
const [packagePath, reportPath] = process.argv.slice(2);
const tsc = process.env.PATH.split(path.delimiter)
    .map(dir => path.join(dir, process.platform === 'win32' ? 'tsc.cmd' : 'tsc'))
    .find(file => fs.existsSync(file));
if (!tsc) throw new Error('TypeScript parser unavailable');
const parserPath = [
    path.join(path.dirname(fs.realpathSync(tsc)), '../lib/typescript.js'),
    path.join(path.dirname(tsc), '../typescript/lib/typescript.js'),
].find(file => fs.existsSync(file));
if (!parserPath) throw new Error('Cannot locate TypeScript beside its tsc executable');
const ts = require(parserPath);
const api = require(packagePath);
if (typeof api.__bindingRegistrySnapshot !== 'function') {
    throw new Error('Build inspection metadata first: python wasm_tsify/build.py --registry');
}
const snapshot = JSON.parse(api.__bindingRegistrySnapshot());
if (snapshot.version !== 2 || !snapshot.records.length) throw new Error('Rebuild registry metadata with the current build.py --registry');
const declarationPath = packagePath.replace(/\.js$/, '.d.ts');
const source = ts.createSourceFile(declarationPath, fs.readFileSync(declarationPath, 'utf8'),
    ts.ScriptTarget.Latest, true, ts.ScriptKind.TS);
if (source.parseDiagnostics.length) throw new Error('Cannot parse generated declarations');
const declarations = new Map();
for (const node of source.statements) {
    if (!node.modifiers?.some(m => m.kind === ts.SyntaxKind.ExportKeyword)) continue;
    if (node.name) declarations.set(node.name.text, node);
    if (ts.isVariableStatement(node)) {
        for (const item of node.declarationList.declarations) declarations.set(item.name.getText(source), item);
    }
}
function descriptor(object, name) {
    while (object) {
        const found = Object.getOwnPropertyDescriptor(object, name);
        if (found) return found;
        object = Object.getPrototypeOf(object);
    }
}
const errors = [];
let matchedRuntime = 0, matchedDeclaration = 0;
const excluded = [];
const nativeProjections = [];
// Public API naming rule, not an independently maintained shortcut inventory.
const conveniences = snapshot.conveniences.map(row => ({...row,
    name: row.name.replace(/_([a-z])/g, (_, letter) => letter.toUpperCase())}));
const rows = [...snapshot.records, ...conveniences];
for (const row of rows) {
    const label = row.owner ? `${row.owner}.${row.projection === 'error_details' ? 'details.' : ''}${row.name}` : row.name;
    const fail = (category, message) => errors.push({id: row.id, path: label, category, message});
    if (row.projection === 'host_filesystem') {
        excluded.push({id: row.id, path: label, reason: 'host_filesystem'});
        continue;
    }
    if (row.projection === 'native') {
        nativeProjections.push({id: row.id, name: row.name, type: row.native_type});
        continue;
    }
    if (row.unresolved_owner || !row.name) {
        fail('registry_metadata', 'No unambiguous JavaScript path');
        continue;
    }
    let object = row.owner ? (row.instance ? api[row.owner]?.prototype : api[row.owner]) : api;
    let parent = declarations.get(row.owner || row.name);
    if (row.projection === 'error_details') {
        const details = parent?.members?.find(member => member.name?.getText(source) === 'details');
        const detailsType = details?.type && ts.isTypeReferenceNode(details.type)
            ? details.type.typeName.getText(source) : undefined;
        if (!descriptor(object, 'details')?.get || !detailsType) {
            fail('shape', 'Error requires a native details getter and a concrete details type');
            continue;
        }
        object = api[detailsType]?.prototype;
        parent = declarations.get(detailsType);
    }
    const constructor = row.owner && !row.instance && row.name === 'new';
    const runtime = constructor ? descriptor(api, row.owner) : descriptor(object, row.name);
    const members = row.owner ? (parent?.members || []) : [];
    const declaration = constructor ? members.find(ts.isConstructorDeclaration) : row.owner ? members.find(member =>
        member.name?.getText(source) === row.name &&
        Boolean(member.modifiers?.some(m => m.kind === ts.SyntaxKind.StaticKeyword)) === !row.instance) : parent;
    if (row.projection !== 'type_only') {
        if (runtime) matchedRuntime++;
        else fail('missing_runtime', 'No export/member descriptor');
    }
    if (declaration) matchedDeclaration++;
    else fail('missing_declaration', 'No exported declaration/member with matching static/instance placement');
    if (row.projection === 'type_only' && declaration &&
        !ts.isTypeAliasDeclaration(declaration) && !ts.isInterfaceDeclaration(declaration)) {
        fail('shape', 'type_only requires a type alias or interface, not a missing runtime class');
    }
    if (!runtime || !declaration) continue;
    if (row.kind === 'type') {
        if (ts.isClassDeclaration(declaration) && typeof runtime.value !== 'function') {
            fail('shape', 'Declared class is not a runtime constructor');
        }
        if (ts.isEnumDeclaration(declaration) && typeof runtime.value !== 'object') {
            fail('shape', 'Declared enum has no runtime enum object');
        }
        continue;
    }
    const method = ts.isMethodDeclaration(declaration) || ts.isMethodSignature(declaration) || ts.isFunctionDeclaration(declaration);
    if (method && typeof runtime.value !== 'function') fail('shape', 'Declaration is callable but runtime descriptor is not');
    if (ts.isGetAccessorDeclaration(declaration) && !runtime.get) fail('shape', 'Declared getter missing at runtime');
    if (ts.isSetAccessorDeclaration(declaration) && !runtime.set) fail('shape', 'Declared setter missing at runtime');
    const setter = members.find(member => ts.isSetAccessorDeclaration(member) &&
        member.name?.getText(source) === row.name &&
        Boolean(member.modifiers?.some(m => m.kind === ts.SyntaxKind.StaticKeyword)) === !row.instance);
    if (setter && !runtime.set && !ts.isSetAccessorDeclaration(declaration)) {
        fail('shape', 'Declared paired setter missing at runtime');
    }
}
// Validate explicit native type expressions in the declaration's own scope.
// No extra exported class is required for a number, union or readonly array.
const original = fs.readFileSync(declarationPath, 'utf8');
const nativeSource = original + nativeProjections.map((p, i) =>
    `\ntype __RegistryNative${i} = ${p.type};`).join('');
const options = {noEmit: true, strict: true, target: ts.ScriptTarget.ESNext};
const host = ts.createCompilerHost(options);
const readFile = host.readFile.bind(host);
host.readFile = file => path.resolve(file) === path.resolve(declarationPath) ? nativeSource : readFile(file);
const program = ts.createProgram([declarationPath], options, host);
const nativeFile = program.getSourceFile(declarationPath);
for (const statement of nativeFile.statements) {
    if (statement.getStart(nativeFile) < original.length) continue;
    const check = node => {
        if (node.kind === ts.SyntaxKind.AnyKeyword || node.kind === ts.SyntaxKind.UnknownKeyword) {
            errors.push({category: 'registry_metadata', path: 'native projections',
                message: 'Native projections must retain concrete types, not any/unknown'});
        }
        ts.forEachChild(node, check);
    };
    check(statement);
}
for (const diagnostic of ts.getPreEmitDiagnostics(program)) {
    if (diagnostic.file?.fileName === declarationPath && diagnostic.start >= original.length) {
        errors.push({category: 'registry_metadata', path: 'native projections',
            message: ts.flattenDiagnosticMessageText(diagnostic.messageText, '\n')});
    }
}
const counts = {};
for (const error of errors) counts[error.category] = (counts[error.category] || 0) + 1;
const report = {package: packagePath, declarations: declarationPath,
    canonical_records: snapshot.records.length, conveniences: snapshot.conveniences.length,
    matched_runtime: matchedRuntime, matched_declaration: matchedDeclaration, counts, errors,
    excluded_host_filesystem: excluded, native_projections: nativeProjections,
    scope: 'Registered JavaScript names, static/instance placement and runtime/declaration shape. Not behavior, defaults, signature equivalence or full Python coverage.'};
fs.mkdirSync(path.dirname(reportPath), {recursive: true});
fs.writeFileSync(reportPath, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({...report, errors: undefined, excluded_host_filesystem: excluded.length,
    native_projections: nativeProjections.length}, null, 2));
for (const error of errors.slice(0, 10)) console.log(`${error.category}: ${error.path}: ${error.message}`);
console.log(`Full report: ${reportPath}`);
process.exitCode = errors.length ? 1 : 0;
