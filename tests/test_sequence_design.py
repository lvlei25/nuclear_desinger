"""
序列设计模块测试
"""

import unittest
from src.sequence_design import NucleicAcidDesigner, SequenceAnalyzer

class TestNucleicAcidDesigner(unittest.TestCase):
    
    def setUp(self):
        self.designer = NucleicAcidDesigner()
    
    def test_complement(self):
        sequence = "ATCG"
        expected = "TAGC"
        self.assertEqual(self.designer.complement(sequence), expected)
    
    def test_reverse_complement(self):
        sequence = "ATCG"
        expected = "CGAT"
        self.assertEqual(self.designer.reverse_complement(sequence), expected)
    
    def test_gc_content(self):
        sequence = "ATCGCG"
        expected = 50.0
        self.assertEqual(self.designer.calculate_gc_content(sequence), expected)
    
    def test_is_valid_sequence(self):
        self.assertTrue(self.designer.is_valid_sequence("ATCG"))
        self.assertFalse(self.designer.is_valid_sequence("ATCGXYZ"))
    
    def test_design_sirna(self):
        mrna = "A" * 200
        sense, antisense = self.designer.design_sirna(mrna, position=100, length=21)
        
        self.assertEqual(len(sense), 21)
        self.assertEqual(len(antisense), 21)
        self.assertEqual(antisense, self.designer.reverse_complement(sense))


class TestSequenceAnalyzer(unittest.TestCase):
    
    def test_find_motifs(self):
        sequence = "ATATCGATAT"
        positions = SequenceAnalyzer.find_motifs(sequence, "ATAT")
        
        self.assertEqual(len(positions), 2)
        self.assertIn(1, positions)
        self.assertIn(7, positions)
    
    def test_calculate_tm(self):
        sequence = "ATCGCGATCG"
        tm = SequenceAnalyzer.calculate_tm(sequence)
        
        self.assertIsInstance(tm, float)
        self.assertGreater(tm, 0)
    
    def test_detect_palindrome(self):
        sequence = "ATCGCGAT"
        palindromes = SequenceAnalyzer.detect_palindrome(sequence)
        
        self.assertIsInstance(palindromes, list)


if __name__ == '__main__':
    unittest.main()