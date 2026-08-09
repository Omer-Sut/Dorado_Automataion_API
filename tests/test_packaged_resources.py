import importlib.util
import unittest
from pathlib import Path

import dorado_workflow
from dorado_workflow.managers.config_manager import ConfigManager


EXPECTED_R_SCRIPTS = {
    "main_analysis_pipeline.R",
    "batch_mapping_analysis.R",
    "batch_methylation_prep.R",
    "batch_nanotel_analysis.R",
}
EXPECTED_R_FUNCTIONS = {
    "mapping_functions.R",
    "methylation_functions.R",
    "nanotel_functions.R",
    "utils.R",
}
EXPECTED_ICONS = {
    "bam.png",
    "basecalling.png",
    "checked.png",
    "fastq.png",
    "folder.png",
    "nanotel.png",
    "pod5.png",
    "sub_basecalling.png",
    "sub_nanotel.png",
}


class PackagedResourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow_root = Path(dorado_workflow.__file__).resolve().parent

        gui_spec = importlib.util.find_spec("dorado_gui")
        if gui_spec is None or gui_spec.submodule_search_locations is None:
            raise AssertionError("dorado_gui package could not be located")
        cls.gui_root = Path(next(iter(gui_spec.submodule_search_locations))).resolve()

    def test_default_config_and_nanotel_are_available(self):
        config = ConfigManager()

        self.assertTrue(Path(config.config_path).is_file())
        self.assertTrue(Path(config.get_nanotel_script_path()).is_file())

    def test_runtime_r_scripts_are_available(self):
        r_root = self.workflow_root / "r_analysis"
        scripts = {path.name for path in r_root.glob("*.R")}
        functions = {path.name for path in (r_root / "functions").glob("*.R")}

        self.assertTrue(EXPECTED_R_SCRIPTS.issubset(scripts))
        self.assertTrue(EXPECTED_R_FUNCTIONS.issubset(functions))

    def test_r_dependency_loader_does_not_install_at_runtime(self):
        utils_source = (
            self.workflow_root / "r_analysis" / "functions" / "utils.R"
        ).read_text(encoding="utf-8")

        self.assertNotIn("install.packages(", utils_source)
        self.assertIn("Missing required R packages:", utils_source)
        self.assertIn("Reinstall or update the Conda environment.", utils_source)

    def test_gui_icons_are_available(self):
        icons = {path.name for path in (self.gui_root / "icons").glob("*.png")}

        self.assertEqual(icons, EXPECTED_ICONS)


if __name__ == "__main__":
    unittest.main()
