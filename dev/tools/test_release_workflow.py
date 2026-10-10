"""Manual Release workflow regressions; no Git, network or publication."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest


WORKFLOW = Path(__file__).resolve().parents[2] / ".github/workflows/release.yml"


def literal_block(step, key):
    lines = WORKFLOW.read_text().split(f"- name: {step}\n", 1)[1].splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip() == f"{key}: |")
    indent = len(lines[start]) - len(lines[start].lstrip())
    block = []
    for line in lines[start + 1:]:
        if line.strip() and len(line) - len(line.lstrip()) <= indent:
            break
        block.append(line)
    return textwrap.dedent("\n".join(block))


class ReleaseWorkflowTests(unittest.TestCase):
    def select(self, tags=(), releases=(), ref="refs/heads/main", api_error=False):
        script = literal_block("Find an unreleased tag on this commit", "script")
        fixture = json.dumps(dict(tags=tags, releases=releases, ref=ref, error=api_error))
        harness = """
const fixture = JSON.parse(process.argv[1]);
const result = {outputs: {}, notices: [], failures: []};
const core = {
  setOutput: (key, value) => {result.outputs[key] = value;},
  notice: message => result.notices.push(message),
  setFailed: message => result.failures.push(message),
  info: () => {},
};
const context = {repo: {owner: 'test', repo: 'test'}, sha: 'current', ref: fixture.ref};
const github = {
  rest: {repos: {listTags: 'tags', listReleases: 'releases'}},
  paginate: async endpoint => {
    if (fixture.error) throw new Error('API unavailable');
    return fixture[endpoint];
  },
};
"""
        harness += "(async () => {\n" + script + "\n})().then(() => console.log(JSON.stringify(result))).catch(error => {console.error(error.message); process.exitCode = 1;});"
        process = subprocess.run(
            ["node", "-e", harness, fixture], capture_output=True, text=True
        )
        if api_error:
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("API unavailable", process.stderr)
            return None
        self.assertEqual(process.returncode, 0, process.stderr)
        return json.loads(process.stdout)

    @staticmethod
    def tag(name="v0.5.0", sha="current"):
        return dict(name=name, commit=dict(sha=sha))

    def test_manual_only_without_tag_input(self):
        workflow = WORKFLOW.read_text()
        self.assertIn("on:\n  workflow_dispatch:\n", workflow)
        self.assertNotIn("  push:", workflow)
        self.assertNotIn("inputs:", workflow)
        self.assertIn("ref: ${{ github.sha }}", workflow)
        self.assertIn('gh release create "$RELEASE_TAG" --verify-tag', workflow)
        self.assertEqual(workflow.count("if: steps.release_tag.outputs.tag != ''"), 3)

    def test_current_commit_unreleased_tag(self):
        result = self.select(tags=[self.tag()])
        self.assertEqual(result["outputs"], {"tag": "v0.5.0"})
        self.assertFalse(result["failures"])

    def test_untagged_or_different_commit_is_skipped(self):
        for tags in ([], [self.tag(sha="other")], [self.tag(name="non-version")]):
            with self.subTest(tags=tags):
                result = self.select(tags=tags)
                self.assertFalse(result["outputs"])
                self.assertTrue(result["notices"])

    def test_existing_releases_are_skipped_including_drafts(self):
        for flags in ({}, {"draft": True}, {"prerelease": True}):
            with self.subTest(flags=flags):
                result = self.select(
                    tags=[self.tag()], releases=[dict(tag_name="v0.5.0", **flags)]
                )
                self.assertFalse(result["outputs"])
                self.assertTrue(result["notices"])

    def test_multiple_unreleased_tags_fail_without_selecting_arbitrarily(self):
        result = self.select(tags=[self.tag(), self.tag(name="v0.5.1")])
        self.assertFalse(result["outputs"])
        self.assertTrue(result["failures"])

    def test_selected_tag_ref_is_respected(self):
        result = self.select(
            tags=[self.tag(), self.tag(name="v0.5.1")], ref="refs/tags/v0.5.1"
        )
        self.assertEqual(result["outputs"], {"tag": "v0.5.1"})

    def test_already_released_tag_does_not_block_another_pending_tag(self):
        result = self.select(
            tags=[self.tag(), self.tag(name="v0.5.1")],
            releases=[dict(tag_name="v0.5.0")],
        )
        self.assertEqual(result["outputs"], {"tag": "v0.5.1"})

    def test_api_failure_does_not_allow_publication(self):
        self.select(api_error=True)

    def test_changelog_uses_detected_tag_not_branch_name(self):
        script = literal_block("Read release notes from changelog", "run")
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "CHANGELOG.md").write_text(
                "<!-- release-header:start -->\nHeader\n<!-- release-header:end -->\n"
                "## [0.5.1] - 2026-10-11\nNewer notes\n"
                "## [0.5.0] - 2026-10-10\nSelected notes\n"
                "## [0.4.0] - 2026-10-01\nOlder notes\n"
                "<!-- release-footer:start -->\nFooter\n<!-- release-footer:end -->\n"
            )
            process = subprocess.run(
                ["bash", "-c", script], cwd=root,
                env={**os.environ, "RUNNER_TEMP": folder,
                     "RELEASE_TAG": "v0.5.0", "GITHUB_REF_NAME": "main"},
                capture_output=True, text=True,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            notes = (root / "release-notes.md").read_text()
            self.assertIn("Selected notes", notes)
            self.assertIn("Header", notes)
            self.assertIn("Footer", notes)
            self.assertNotIn("Newer notes", notes)
            self.assertNotIn("Older notes", notes)


if __name__ == "__main__":
    unittest.main()
