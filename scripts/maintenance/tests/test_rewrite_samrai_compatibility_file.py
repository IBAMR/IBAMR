#!/usr/bin/env python3

import importlib.util
import pathlib
import sys
import tempfile
import unittest


SCRIPT_PATH = pathlib.Path(__file__).resolve().parents[1] / "rewrite_samrai_compatibility_file.py"
SPEC = importlib.util.spec_from_file_location("rewrite_samrai_compatibility_file", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

AliasEntry = MODULE.AliasEntry


class RewriteSamraiCompatibilityFileTests(unittest.TestCase):
    def test_transform_preserves_string_literals_while_rewriting_code_tokens(self):
        entries = [
            AliasEntry(
                alias_name="SAMRAIHyperbolicLevelIntegrator",
                include_name="SAMRAIHyperbolicLevelIntegrator.h",
                include_keys={"HyperbolicLevelIntegrator.h"},
                base_types={"HyperbolicLevelIntegrator"},
                is_template_alias=False,
            ),
            AliasEntry(
                alias_name="SAMRAIPointer",
                include_name="SAMRAIPointer.h",
                include_keys={"Pointer.h"},
                base_types={"Pointer"},
                is_template_alias=True,
            ),
        ]

        old_text = """#include \"HyperbolicLevelIntegrator.h\"\n#include \"Pointer.h\"\n\nvoid f()\n{\n    Pointer<HyperbolicLevelIntegrator<NDIM> > p;\n    const char* key = \"HyperbolicLevelIntegrator\";\n    std::string obj_name = \"::HyperbolicLevelIntegrator\";\n}\n"""

        with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as tmp:
            tmp.write(old_text)
            tmp_path = pathlib.Path(tmp.name)

        try:
            _, new_text = MODULE.transform_file(tmp_path, pathlib.Path("/tmp"), entries, aggressive_unqualified=False)
        finally:
            tmp_path.unlink(missing_ok=True)

        self.assertIn("SAMRAIPointer<SAMRAIHyperbolicLevelIntegrator", new_text)
        self.assertIn('"HyperbolicLevelIntegrator"', new_text)
        self.assertIn('"::HyperbolicLevelIntegrator"', new_text)
        self.assertNotIn('"SAMRAIHyperbolicLevelIntegrator"', new_text)
        self.assertNotIn('"::SAMRAIHyperbolicLevelIntegrator"', new_text)

    def test_remove_ibtk_alias_qualification_skips_non_code_regions(self):
        text = (
            'std::string s = "IBTK::SAMRAIPointer";\n'
            "// IBTK::SAMRAIPointer should stay in comment\n"
            "IBTK::SAMRAIPointer<int> p;\n"
        )

        new_text, changed = MODULE._remove_ibtk_alias_qualification(text, {"SAMRAIPointer"})

        self.assertTrue(changed)
        self.assertIn('std::string s = "IBTK::SAMRAIPointer";', new_text)
        self.assertIn('// IBTK::SAMRAIPointer should stay in comment', new_text)
        self.assertIn('SAMRAIPointer<int> p;', new_text)
        self.assertNotIn('IBTK::SAMRAIPointer<int> p;', new_text)

    def test_collect_alias_includes_ignores_strings_and_comments(self):
        text = (
            'std::string s = "SAMRAIPointer";\n'
            "// SAMRAIHyperbolicLevelIntegrator in comment\n"
            "SAMRAIHyperbolicLevelIntegrator* p = nullptr;\n"
        )

        includes = MODULE._collect_alias_includes_from_text(
            text,
            {
                "SAMRAIPointer": "SAMRAIPointer.h",
                "SAMRAIHyperbolicLevelIntegrator": "SAMRAIHyperbolicLevelIntegrator.h",
            },
        )

        self.assertEqual(includes, {"SAMRAIHyperbolicLevelIntegrator.h"})


if __name__ == "__main__":
    unittest.main()
