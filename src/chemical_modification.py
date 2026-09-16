"""
小核酸药物化学修饰模块
包含各种化学修饰类型和修饰策略
"""

from typing import List, Dict, Optional, Tuple

class ChemicalModifier:
    """
    化学修饰器
    """

    MODIFICATION_TYPES = {
        '2OMe': {'name': '2-O-甲基修饰', 'position': '糖环2\'位', 'effect': '增加稳定性，降低免疫原性', 'symbol': 'm', 'synth_format': 'm{N}'},
        'PS': {'name': '硫代磷酸酯修饰', 'position': '磷酸骨架', 'effect': '增加核酸酶抗性', 'symbol': '*', 'synth_format': '{N}p'},
        'LNA': {'name': '锁核酸修饰', 'position': '糖环', 'effect': '提高结合亲和力，增加稳定性', 'symbol': '+', 'synth_format': '+{N}'},
        'FANA': {'name': '2\'氟代阿拉伯糖核酸', 'position': '糖环2\'位', 'effect': '提高RNA亲和力', 'symbol': 'F', 'synth_format': 'F{N}'},
        'PMO': {'name': '肽核酸类似物', 'position': '骨架', 'effect': '高特异性，抗降解', 'symbol': 'P', 'synth_format': 'P{N}'},
        'cEt': {'name': '2\'氟代乙氧基修饰', 'position': '糖环2\'位', 'effect': '增强效力和稳定性', 'symbol': 'e', 'synth_format': 'e{N}'},
        'GalNAc': {'name': 'N-乙酰半乳糖胺修饰', 'position': '3\'端', 'effect': '靶向肝细胞', 'symbol': 'G', 'synth_format': 'GalNAc-{N}'},
        'Cholesterol': {'name': '胆固醇修饰', 'position': '末端', 'effect': '增加细胞膜通透性', 'symbol': 'C', 'synth_format': 'Chol-{N}'}
    }

    SUGAR_MODS = ['2OMe', 'FANA', 'LNA', 'cEt']
    BACKBONE_MODS = ['PS']

    def __init__(self):
        self.modifications = []

    def add_modification(self, position: int, mod_type: str) -> None:
        """
        添加修饰
        :param position: 修饰位置（1-based）
        :param mod_type: 修饰类型
        """
        if mod_type not in self.MODIFICATION_TYPES:
            raise ValueError(f"未知修饰类型: {mod_type}")

        self.modifications.append({
            'position': position,
            'type': mod_type,
            'info': self.MODIFICATION_TYPES[mod_type]
        })

    def remove_modification(self, position: int) -> None:
        """
        移除指定位置的修饰
        """
        self.modifications = [m for m in self.modifications if m['position'] != position]

    def get_modifications(self) -> List[Dict]:
        """
        获取所有修饰
        """
        return self.modifications

    def clear_modifications(self) -> None:
        """
        清除所有修饰
        """
        self.modifications = []

    def apply_modifications(self, sequence: str, format_type: str = 'mark') -> str:
        """
        将修饰应用到序列上
        :param sequence: 原始序列
        :param format_type: 输出格式 ('mark' - 标记形式, 'synth' - 合成订单格式)
        :return: 带修饰的序列
        """
        modified = list(sequence)

        for mod in self.modifications:
            pos = mod['position'] - 1
            if pos < len(modified):
                base = modified[pos]
                mod_type = mod['type']
                mod_info = self.MODIFICATION_TYPES.get(mod_type, {})

                if format_type == 'synth':
                    synth_format = mod_info.get('synth_format', '{N}*')
                    modified[pos] = synth_format.replace('{N}', base)
                else:
                    symbol = mod_info.get('symbol', '*')
                    modified[pos] = f"{base}{symbol}"

        return ''.join(modified)

    def suggest_modifications(self, sequence_type: str = 'siRNA') -> List[Dict]:
        """
        根据序列类型推荐修饰策略
        :param sequence_type: 序列类型 (siRNA, ASO, miRNA)
        :return: 推荐的修饰列表
        """
        suggestions = []

        if sequence_type == 'siRNA':
            suggestions = [
                {'position': 2, 'type': '2OMe', 'reason': '增强稳定性'},
                {'position': 14, 'type': '2OMe', 'reason': '增强稳定性'},
                {'position': 1, 'type': 'PS', 'reason': '增强核酸酶抗性'},
                {'position': 21, 'type': 'PS', 'reason': '增强核酸酶抗性'}
            ]
        elif sequence_type == 'ASO':
            suggestions = [
                {'position': 1, 'type': 'PS', 'reason': '增强核酸酶抗性'},
                {'position': 2, 'type': 'PS', 'reason': '增强核酸酶抗性'},
                {'position': 3, 'type': 'LNA', 'reason': '提高结合亲和力'},
                {'position': len(self.modifications)-2, 'type': 'LNA', 'reason': '提高结合亲和力'},
                {'position': len(self.modifications)-1, 'type': 'PS', 'reason': '增强核酸酶抗性'}
            ]
        elif sequence_type == 'miRNA':
            suggestions = [
                {'position': 2, 'type': '2OMe', 'reason': '增强稳定性，降低脱靶'},
                {'position': 3, 'type': '2OMe', 'reason': '增强稳定性'},
                {'position': 5, 'type': 'PS', 'reason': '增强核酸酶抗性'}
            ]

        return suggestions

    def design_siRNA_modification(self, sense: str, antisense: str, platform: str = 'GalNAc', duplex_config: str = None) -> Dict:
        """
        根据10步算法设计siRNA化学修饰

        Step 1: 平台选择 (tissue → platform GalNAc/LNP)
        Step 2: 双链构型 (自动识别 19+21, 21+23, 17+17...)
        Step 3: 初始化修饰图谱 - 创建Nuc数组: {base, pos, sugar=null, ps_after=false}
        Step 4: 全链2'-OMe - GalNAc:双链全覆盖; LNP:奇数位选择性覆盖
        Step 5: 2'-F覆盖写 - AS first 25% + last 25%, 排除切割位点中心15%
        Step 6: PS骨架 - SS: [1,2,3] + 50-65%区间; AS: [1,2,3] + last 3
        Step 7: 种子区验证 - AS pos 2-8 必须有sugar (兜底2'-OMe)
        Step 8: 末端标注 - AS 5'磷酸化 + SS 3'偶联(GalNAc/dTdT)
        Step 9: 序列化输出 - m=2'-OMe, f=2'-F, p=PS键
        Step 10: 校验 - 6项检查: 种子覆盖/切割位点/PS密度/5'磷酸/覆盖率/长度一致性

        :param sense: 正义链序列 (SS)
        :param antisense: 反义链序列 (AS)
        :param platform: 递送平台 ('GalNAc' 或 'LNP')
        :return: 修饰结果字典
        """
        sense_upper = sense.upper().replace('T', 'U')
        antisense_upper = antisense.upper().replace('T', 'U')

        ss_len = len(sense_upper)
        as_len = len(antisense_upper)

        duplex_config = f"{ss_len}+{as_len}"

        ss = [{'base': sense_upper[i], 'pos': i+1, 'sugar': None, 'ps_after': False} for i in range(ss_len)]
        as_seq = [{'base': antisense_upper[i], 'pos': i+1, 'sugar': None, 'ps_after': False} for i in range(as_len)]

        if platform == 'GalNAc':
            for i in range(ss_len):
                ss[i]['sugar'] = '2OMe'
            for i in range(as_len):
                as_seq[i]['sugar'] = '2OMe'
        else:
            for i in range(ss_len):
                if (i + 1) % 2 == 1:
                    ss[i]['sugar'] = '2OMe'
            for i in range(as_len):
                if (i + 1) % 2 == 1:
                    as_seq[i]['sugar'] = '2OMe'

        if as_len > 0:
            first_25 = int(as_len * 0.25)
            last_25 = int(as_len * 0.25)
            cut_center_start = int(as_len * 0.425)
            cut_center_end = int(as_len * 0.575)

            for i in range(first_25):
                if i < cut_center_start or i >= cut_center_end:
                    as_seq[i]['sugar'] = 'FANA'

            for i in range(as_len - last_25, as_len):
                if i < cut_center_start or i >= cut_center_end:
                    as_seq[i]['sugar'] = 'FANA'

        if ss_len > 0:
            for i in [0, 1, 2]:
                if i < ss_len:
                    ss[i]['ps_after'] = True

            start_50 = int(ss_len * 0.50)
            end_65 = int(ss_len * 0.65)
            for i in range(start_50, min(end_65, ss_len)):
                ss[i]['ps_after'] = True

        if as_len > 0:
            for i in [0, 1, 2]:
                if i < as_len:
                    as_seq[i]['ps_after'] = True

            for i in range(max(0, as_len - 3), as_len):
                as_seq[i]['ps_after'] = True

        for i in range(1, 8):
            if i <= as_len and as_seq[i-1]['sugar'] is None:
                as_seq[i-1]['sugar'] = '2OMe'

        as_5p = True
        ss_3conjugate = 'GalNAc' if platform == 'GalNAc' else 'dTdT'

        sense_output = self._map_to_output_new(ss, False, None)
        antisense_output = self._map_to_output_new(as_seq, as_5p, None)

        if platform == 'GalNAc':
            sense_output += '-GalNAc'
        else:
            sense_output += '-dTdT'

        validation = self._validate_modification_new(ss, as_seq, platform)

        return {
            'platform': platform,
            'duplex_config': duplex_config,
            'sense_seq': sense,
            'antisense_seq': antisense,
            'sense_modified': sense_output,
            'antisense_modified': antisense_output,
            'sense_map': ss,
            'antisense_map': as_seq,
            'validation': validation,
            'modifications': []
        }

    def _map_to_output_new(self, mod_map: List[Dict], five_prime_phos: bool, three_conjugate: Optional[str]) -> str:
        """将修饰图谱转换为合成订单格式输出"""
        output = []

        if five_prime_phos:
            output.append('p')

        for i, m in enumerate(mod_map):
            base = m['base']
            sugar = m.get('sugar')
            ps_after = m.get('ps_after', False)

            if sugar == '2OMe':
                output.append(f"m{base}")
            elif sugar == 'FANA':
                output.append(f"f{base}")
            elif sugar == 'LNA':
                output.append(f"+{base}")
            elif sugar == 'cEt':
                output.append(f"e{base}")
            else:
                output.append(base)

            if ps_after and i < len(mod_map) - 1:
                output[-1] += 'p'

        if three_conjugate:
            output.append(f"-{three_conjugate}")

        return ''.join(output)

    def _validate_modification_new(self, ss_map: List[Dict], as_map: List[Dict], platform: str) -> Dict:
        """校验修饰方案 - 6项检查"""
        validation = {
            'valid': True,
            'checks': [],
            'errors': [],
            'warnings': []
        }

        seed_coverage = sum(1 for i in range(1, 8) if i <= len(as_map) and as_map[i-1]['sugar'] is not None)
        seed_ok = seed_coverage == 7
        validation['checks'].append(f"种子区覆盖: {seed_coverage}/7 {'✓' if seed_ok else '✗'}")
        if not seed_ok:
            validation['valid'] = False
            validation['errors'].append(f"种子区覆盖不完整: {seed_coverage}/7")

        cut_center_start = int(len(as_map) * 0.425)
        cut_center_end = int(len(as_map) * 0.575)
        f_in_cut = sum(1 for i in range(cut_center_start, cut_center_end) if as_map[i].get('sugar') == 'FANA')
        cut_ok = f_in_cut == 0
        validation['checks'].append(f"切割位点2'-F排除: {f_in_cut}个 {'✓' if cut_ok else '✗'}")
        if not cut_ok:
            validation['valid'] = False
            validation['errors'].append(f"切割位点中心区域存在2'-F修饰: {f_in_cut}个")

        ss_ps_count = sum(1 for m in ss_map if m.get('ps_after'))
        as_ps_count = sum(1 for m in as_map if m.get('ps_after'))
        ss_ps_ok = ss_ps_count <= 8
        as_ps_ok = as_ps_count <= 8
        validation['checks'].append(f"SS链PS密度: {ss_ps_count}/8 {'✓' if ss_ps_ok else '✗'}")
        validation['checks'].append(f"AS链PS密度: {as_ps_count}/8 {'✓' if as_ps_ok else '✗'}")
        if not ss_ps_ok:
            validation['warnings'].append(f"SS链PS修饰过多: {ss_ps_count}/8")
        if not as_ps_ok:
            validation['warnings'].append(f"AS链PS修饰过多: {as_ps_count}/8")

        validation['checks'].append(f"AS 5'磷酸化: ✓")

        ss_sugar_coverage = sum(1 for m in ss_map if m.get('sugar') is not None) / len(ss_map) * 100 if ss_map else 0
        as_sugar_coverage = sum(1 for m in as_map if m.get('sugar') is not None) / len(as_map) * 100 if as_map else 0
        coverage_ok = ss_sugar_coverage >= 70 and as_sugar_coverage >= 70
        validation['checks'].append(f"糖修饰覆盖率: SS={ss_sugar_coverage:.1f}%, AS={as_sugar_coverage:.1f}% {'✓' if coverage_ok else '✗'}")
        if not coverage_ok:
            validation['warnings'].append(f"糖修饰覆盖率不足: SS={ss_sugar_coverage:.1f}%, AS={as_sugar_coverage:.1f}%")

        len_consistent = len(ss_map) == len(as_map) or abs(len(ss_map) - len(as_map)) <= 4
        validation['checks'].append(f"双链长度一致性: {len(ss_map)}+{len(as_map)} {'✓' if len_consistent else '✗'}")
        if not len_consistent:
            validation['warnings'].append(f"双链长度差异较大: {len(ss_map)} vs {len(as_map)}")

        return validation


class DeliverySystem:
    """
    递送系统设计
    """

    DELIVERY_METHODS = {
        'LNP': {'name': '脂质纳米颗粒', 'target': '全身性', 'advantage': '高效递送，可规模化生产'},
        'GalNAc': {'name': 'GalNAc偶联', 'target': '肝脏', 'advantage': '高特异性，无需载体'},
        'AAV': {'name': '腺相关病毒', 'target': '多种组织', 'advantage': '长效表达'},
        'PEI': {'name': '聚乙烯亚胺', 'target': '体外转染', 'advantage': '高效转染'},
        'GoldNP': {'name': '金纳米颗粒', 'target': '肿瘤', 'advantage': '高负载量，生物相容性'}
    }

    @staticmethod
    def recommend_delivery(target_tissue: str) -> Dict:
        """
        根据目标组织推荐递送系统
        """
        recommendations = {
            '肝脏': ['GalNAc', 'LNP'],
            '眼睛': ['LNP', 'AAV'],
            '中枢神经系统': ['AAV', 'LNP'],
            '肿瘤': ['LNP', 'GoldNP'],
            '肺部': ['LNP', 'AAV'],
            '肌肉': ['AAV', 'LNP']
        }

        return recommendations.get(target_tissue, ['LNP'])

    @staticmethod
    def calculate_dosage(sequence_length: int, target_tissue: str) -> float:
        """
        计算推荐剂量（mg/kg）
        """
        base_dosage = {
            '肝脏': 1.0,
            '眼睛': 0.5,
            '中枢神经系统': 2.0,
            '肿瘤': 1.5,
            '肺部': 1.0,
            '肌肉': 1.0
        }

        length_factor = min(sequence_length / 20, 2.0)

        return base_dosage.get(target_tissue, 1.0) * length_factor