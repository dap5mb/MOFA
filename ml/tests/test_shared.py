import glob
import os
import re
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
import shared

ML_SCRIPTS = os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts"))


class TestSharedPreprocessing(unittest.TestCase):
    def test_exports_come_from_mofa_scripts(self):
        import plate_correction
        self.assertIs(shared.correct_plate, plate_correction.correct_plate)
        self.assertTrue(os.path.samefile(os.path.dirname(plate_correction.__file__), shared.MOFA_SCRIPTS))

    def test_no_copies_under_ml(self):
        # The shared functions are defined in mofa/scripts only
        for path in glob.glob(os.path.join(ML_SCRIPTS, "*.py")):
            src = open(path).read()
            for name in shared.__all__:
                self.assertIsNone(re.search(rf"^\s*def {name}\(", src, re.M),
                                  f"{name} is defined in {os.path.basename(path)}; import it from shared")

    def test_unknown_name_raises(self):
        with self.assertRaises(AttributeError):
            shared.not_a_function


if __name__ == "__main__":
    unittest.main()
