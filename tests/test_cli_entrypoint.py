import subprocess
import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class CliEntrypointTests(unittest.TestCase):
    def run_python(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, *arguments],
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            check=False,
            text=True,
        )

    def test_package_module_displays_version(self):
        result = self.run_python("-m", "dorado_workflow.main", "--version")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Dorado Workflow v2.0", result.stdout)

    def test_direct_script_displays_version(self):
        result = self.run_python("dorado_workflow/main.py", "--version")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Dorado Workflow v2.0", result.stdout)

    def test_package_module_displays_help(self):
        result = self.run_python("-m", "dorado_workflow.main", "--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Dorado Workflow", result.stdout)
        self.assertIn("r-analysis", result.stdout)


if __name__ == "__main__":
    unittest.main()
