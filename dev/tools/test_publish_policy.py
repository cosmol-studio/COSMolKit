"""Release-policy regressions; no publication or remote calls."""

import unittest

from publish_policy import check_versions


class PublishPolicyTests(unittest.TestCase):
    def test_manual_rc_publish_and_dry_run(self):
        for requested in (False, True):
            self.assertTrue(check_versions("0.5.0-rc.13", "0.5.0rc13", "0.5.0-rc.13", "workflow_dispatch", "refs/heads/main", requested))

    def test_manual_stable_build_but_not_publish(self):
        self.assertFalse(check_versions("0.5.0", "0.5.0", "0.5.0", "workflow_dispatch", "refs/heads/main", False))
        with self.assertRaisesRegex(ValueError, "RC-only"):
            check_versions("0.5.0", "0.5.0", "0.5.0", "workflow_dispatch", "refs/heads/main", True)

    def test_matching_tags_allow_stable_and_rc(self):
        self.assertFalse(check_versions("0.5.0", "0.5.0", "0.5.0", "push", "refs/tags/v0.5.0", False))
        self.assertTrue(check_versions("0.5.0-rc.13", "0.5.0rc13", "0.5.0-rc.13", "push", "refs/tags/v0.5.0-rc.13", False))

    def test_mismatched_versions_and_tags_fail(self):
        for python, wasm in (("0.5.0rc12", "0.5.0-rc.13"), ("0.5.0rc13", "0.5.0-rc.12")):
            with self.subTest(python=python, wasm=wasm), self.assertRaisesRegex(ValueError, "must match"):
                check_versions("0.5.0-rc.13", python, wasm, "workflow_dispatch", "", True)
        with self.assertRaisesRegex(ValueError, "Release tag"):
            check_versions("0.5.0", "0.5.0", "0.5.0", "push", "refs/tags/v0.4.0", False)

    def test_other_prerelease_versions_are_not_rc(self):
        with self.assertRaisesRegex(ValueError, "stable or RC"):
            check_versions("0.5.0-beta.1", "0.5.0b1", "0.5.0-beta.1", "workflow_dispatch", "", True)


if __name__ == "__main__":
    unittest.main()
