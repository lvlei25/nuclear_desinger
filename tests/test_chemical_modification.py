"""
化学修饰模块测试
"""

import unittest
from src.chemical_modification import ChemicalModifier, DeliverySystem

class TestChemicalModifier(unittest.TestCase):
    
    def setUp(self):
        self.modifier = ChemicalModifier()
    
    def test_add_modification(self):
        self.modifier.add_modification(2, '2OMe')
        
        modifications = self.modifier.get_modifications()
        self.assertEqual(len(modifications), 1)
        self.assertEqual(modifications[0]['position'], 2)
        self.assertEqual(modifications[0]['type'], '2OMe')
    
    def test_apply_modifications(self):
        self.modifier.add_modification(2, '2OMe')
        result = self.modifier.apply_modifications("ATCG")
        
        self.assertEqual(result, "A[2OMe]TCG")
    
    def test_suggest_modifications(self):
        suggestions = self.modifier.suggest_modifications('siRNA')
        
        self.assertIsInstance(suggestions, list)
        self.assertGreater(len(suggestions), 0)


class TestDeliverySystem(unittest.TestCase):
    
    def test_recommend_delivery(self):
        methods = DeliverySystem.recommend_delivery('肝脏')
        
        self.assertIsInstance(methods, list)
        self.assertIn('GalNAc', methods)
    
    def test_calculate_dosage(self):
        dosage = DeliverySystem.calculate_dosage(21, '肝脏')
        
        self.assertIsInstance(dosage, float)
        self.assertGreater(dosage, 0)


if __name__ == '__main__':
    unittest.main()