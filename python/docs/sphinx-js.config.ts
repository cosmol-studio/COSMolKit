import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { ReflectionKind } from "typedoc";

// sphinx-js 5 skips Enum/EnumMember (convertTopLevel.ts). Keep TypeDoc as
// the parser and add only the missing rendering, without rewriting declarations.
export const config = {
  postConvert(_app, project) {
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
    const path = fileURLToPath(new URL("../../target/docs-js-package/enums.rst", import.meta.url));
    const content = lines.join("\n");
    // Keep incremental guide edits from invalidating the complete API page.
    if (!existsSync(path) || readFileSync(path, "utf8") !== content) {
      writeFileSync(path, content);
    }
  },
};
