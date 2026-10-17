"""
配置文件
"""

from dataclasses import dataclass


# 版本号（与 src/__init__.py、根目录 CHANGELOG.md 保持同步）
__version__ = '1.1.0'


@dataclass
class ScoringConfig:
    """siRNA 评分与过滤的全部可调参数

    默认值与 src/sequence_design.py 的历史硬编码值完全一致，
    外置目的仅为集中管理与调参，不改变默认计算结果。
    """

    # ---- 硬性过滤 (hard_filter) ----
    hf1_gc_min: float = 25.0            # HF1: GC含量下限 (%)
    hf1_gc_max: float = 65.0            # HF1: GC含量上限 (%)
    hf2_poly_run: int = 5               # HF2: 禁止 ≥N 连续相同碱基
    hf3_reynolds_min: int = 2           # HF3: Reynolds 评分下限

    # ---- 脱靶检测 (check_off_target / get_specificity_score) ----
    off_target_identity_min: float = 80.0    # 判定脱靶的最小相似度 (%，含)
    off_target_identity_max: float = 100.0   # 完全匹配不计入脱靶 (%，不含)
    off_target_penalty_per_match: float = 5.0    # 每处脱靶扣分
    off_target_penalty_cap: float = 50.0     # 脱靶扣分上限

    # ---- BLAST 比对 (blast_alignment) ----
    seed_size: int = 7                  # 种子区长度
    seed_flank: int = 5                 # 种子匹配后向两侧扩展的比对窗口

    # ---- 评分规则公共阈值 ----
    rule_gc_min: float = 30.0           # Reynolds/Amarzguioui/Ui-Tei 的最优 GC 下限 (%)
    rule_gc_max: float = 55.0           # 最优 GC 上限 (%)
    uitei_dg_threshold: float = -5.0    # UT5: 正义链5'端 ΔG 下限 (kcal/mol)

    # ---- 综合评分权重 (comprehensive_sirna_evaluation) ----
    hard_pass_bonus: float = 30.0       # 通过全部硬过滤的基础分
    reynolds_weight: float = 4.0        # Reynolds 每分权重
    amarzguioui_weight: float = 3.0     # Amarzguioui 每分权重
    seed_complexity_weight: float = 8.0     # 种子区复杂度（0~1）权重
    thermo_asymmetry_bonus: float = 3.0     # 热力学不对称性通过加分
    uitei_pass_bonus: float = 3.0       # Ui-Tei 全部通过加分
    optimal_gc_bonus: float = 5.0       # GC 落入最优区间加分
    max_score: float = 96.0             # 满分
    recommend_threshold: float = 70.0   # "推荐"分级下限
    good_threshold: float = 50.0        # "良好"分级下限

    # ---- 靶区域加分 (get_target_region_bonus) ----
    utr3_bonus: float = 10.0            # 3'-UTR 加分
    cds_bonus: float = 5.0              # CDS 加分
    utr5_bonus: float = 0.0             # 5'-UTR 加分
    start_codon_penalty: float = -5.0   # 起始密码子附近减分
    cds_start_fraction: float = 0.25    # CDS 起点占 mRNA 比例
    cds_end_fraction: float = 0.75      # CDS 终点占 mRNA 比例
    start_codon_window: int = 50        # 起始密码子前后判定窗口 (nt)

    # ---- 默认全局配置实例 ----
    @classmethod
    def default(cls) -> 'ScoringConfig':
        return cls()


# 序列设计参数
SEQUENCE_DESIGN = {
    'default_sirna_length': 21,
    'default_gc_content': 50.0,
    'min_gc_content': 30.0,
    'max_gc_content': 70.0,
    'target_region_start': 100,
    'target_region_length': 500
}

# 修饰参数
MODIFICATION = {
    'default_modifications': ['2OMe', 'PS', 'LNA'],
    'max_modifications': 5,
    'modification_cost': {
        '2OMe': 1.0,
        'PS': 0.5,
        'LNA': 2.0,
        'FANA': 1.5,
        'PMO': 3.0,
        'cEt': 2.5,
        'GalNAc': 4.0,
        'Cholesterol': 2.0
    }
}

# 数据库路径
PATHS = {
    'sequence_db': 'data/sequences',
    'target_genes': 'data/target_genes.csv',
    'design_results': 'data/design_results.csv',
    'logs': 'logs/'
}

# 输出设置
OUTPUT = {
    'format': 'csv',  # csv, json, excel
    'include_analysis': True,
    'include_modifications': True,
    'include_off_target': True
}
