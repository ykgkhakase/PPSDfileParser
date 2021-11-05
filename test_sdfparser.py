import unittest
import sdfparser as sdf

class TestSDfileParser(unittest.TestCase):

	def test_PartiallyParseFromFile(self):
		file1 = 'data/PubChem_compound_text_egfr_records.sdf'
		mols1 = sdf.SDFileParser(file1, maxNumOfMol=10)
		self.assertEqual(len(mols1), 10)

		file2 = 'data/PubChem_compound_text_protease_records.sdf.gz'
		mols2 = sdf.SDFileParser(file2, maxNumOfMol=30)
		self.assertEqual(len(mols2), 30)

	def test_SkipFirstNMoleculesFromFile(self):
		file1 = 'data/PubChem_compound_text_egfr_records.sdf'
		mols1 = sdf.SDFileParser(file1, maxNumOfMol=10, numOfSkippedMol=5)
		self.assertEqual(len(mols1), 10)


if __name__ == "__main__":
    unittest.main()