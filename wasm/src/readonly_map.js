// A ReadonlyMap receiver, not a writable Map with hidden mutator methods.
// Binding declarations remain generated from Rust; this is shared JS ownership
// adaptation only, with no API- or chemistry-specific behavior.
export function readonlyMap(map) {
    const reject = () => { throw new TypeError("readonly snapshot; replace the whole field"); };
    return new Proxy(map, {
        get(target, key, receiver) {
            if (key === "size") return target.size;
            if (key === "set" || key === "delete" || key === "clear") return reject;
            if (key === "get" || key === "has" || key === "keys" ||
                key === "values" || key === "entries" || key === Symbol.iterator) {
                return target[key].bind(target);
            }
            if (key === "forEach") {
                return (callback, thisArg) => target.forEach(
                    (value, key) => callback.call(thisArg, value, key, receiver)
                );
            }
            return Reflect.get(target, key, receiver);
        },
        set: reject,
        defineProperty: reject,
        deleteProperty: reject,
        setPrototypeOf: reject,
    });
}
