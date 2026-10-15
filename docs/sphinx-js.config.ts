import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { ReflectionKind } from "typedoc";

// sphinx-js 5 skips Enum/EnumMember (convertTopLevel.ts). Keep TypeDoc as
// the parser and add only the missing rendering, without rewriting declarations.
export const config = {
  postConvert(_app, project, typedocToIRMap) {
    // sphinx-js treats `const alias: typeof Class.factory` as an attribute.
    // Display the resolved callable's actual signature, keeping the alias's
    // own name, anchor and location. No API names or signatures are duplicated.
    for (const [declaration, entry] of typedocToIRMap) {
      if (!declaration.kindOf(ReflectionKind.Variable) || declaration.type?.type !== "query") continue;
      const target = declaration.type.queryType.reflection;
      const original = target && typedocToIRMap.get(target);
      if (original?.kind !== "function") continue;
      Object.assign(entry, {
        kind: "function", params: original.params, returns: original.returns,
        type_params: original.type_params, exceptions: original.exceptions,
        is_async: original.is_async,
      });
      delete entry.type;
      delete entry.readonly;
    }
    // sphinx-js hides explicit void results; keep every declared return type.
    function retainVoid(declaration, entry) {
      const signature = declaration.signatures?.[0]
        ?? (declaration.type?.type === "reflection" ? declaration.type.declaration.signatures?.[0] : undefined);
      const result = signature?.type;
      if (entry.kind === "function" && !declaration.kindOf(ReflectionKind.Constructor)
          && result?.type === "intrinsic" && result.name === "void") {
        entry.returns = [{ type: ["void"], description: entry.block_tags.returns?.[0] ?? [] }];
      }
      const children = new Map((declaration.children ?? []).map(child => [child.name, child]));
      for (const member of entry.members ?? []) {
        const child = children.get(member.name);
        if (child) retainVoid(child, member);
      }
    }
    for (const [declaration, entry] of typedocToIRMap) {
      retainVoid(declaration, entry);
    }
    const lines = ["Enumerations", "------------", ""];
    for (const enumeration of project.getReflectionsByKind(ReflectionKind.Enum)) {
      lines.push(`.. js:attribute:: ${enumeration.name}`, "", "   TypeScript enum.", "");
      for (const member of enumeration.children ?? []) {
        if (member.type?.type !== "literal") {
          throw new Error(`Expected TypeDoc literal value for ${enumeration.name}.${member.name}`);
        }
        lines.push(`.. js:attribute:: ${enumeration.name}.${member.name}`, "",
          `   Value: \`\`${JSON.stringify(member.type.value)}\`\`.`, "");
      }
    }
    const path = fileURLToPath(new URL("../target/docs-js-package/enums.rst", import.meta.url));
    const content = lines.join("\n");
    // Keep incremental guide edits from invalidating the complete API page.
    if (!existsSync(path) || readFileSync(path, "utf8") !== content) {
      writeFileSync(path, content);
    }
  },
};
