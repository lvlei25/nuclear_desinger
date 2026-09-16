"""
配置文件
"""

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