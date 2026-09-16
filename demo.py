"""
小核酸药物设计演示脚本
"""

from src.sequence_design import NucleicAcidDesigner, SequenceAnalyzer
from src.chemical_modification import ChemicalModifier, DeliverySystem
from src.data_processing import DataHandler, SequenceDatabase

def main():
    print("=" * 60)
    print("        小核酸药物设计工具包演示")
    print("=" * 60)
    
    # 1. 序列设计演示
    print("\n【模块1】序列设计")
    print("-" * 40)
    
    designer = NucleicAcidDesigner()
    
    # 目标mRNA序列（示例）
    mrna_sequence = "AUGGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGC"
    
    # 设计siRNA
    sense, antisense = designer.design_sirna(mrna_sequence, position=5, length=21)
    
    print(f"目标mRNA序列 (前50nt): {mrna_sequence[:50]}...")
    print(f"\n设计的siRNA:")
    print(f"  正义链: {sense}")
    print(f"  反义链: {antisense}")
    
    # 计算GC含量
    gc_content = designer.calculate_gc_content(sense)
    print(f"\nGC含量: {gc_content:.2f}%")
    
    # 序列分析
    analyzer = SequenceAnalyzer()
    tm = analyzer.calculate_tm(sense)
    print(f"熔解温度(Tm): {tm:.2f}°C")
    
    # 2. 化学修饰演示
    print("\n\n【模块2】化学修饰")
    print("-" * 40)
    
    modifier = ChemicalModifier()
    
    # 添加修饰
    modifier.add_modification(2, '2OMe')
    modifier.add_modification(14, '2OMe')
    modifier.add_modification(1, 'PS')
    modifier.add_modification(21, 'PS')
    
    modified = modifier.apply_modifications(sense)
    print(f"修饰后序列: {modified}")
    
    # 推荐修饰策略
    suggestions = modifier.suggest_modifications('siRNA')
    print("\n推荐修饰策略:")
    for i, suggestion in enumerate(suggestions, 1):
        mod_info = modifier.MODIFICATION_TYPES[suggestion['type']]
        print(f"  {i}. 位置{suggestion['position']}: {mod_info['name']} - {suggestion['reason']}")
    
    # 3. 递送系统演示
    print("\n\n【模块3】递送系统")
    print("-" * 40)
    
    target_tissue = "肝脏"
    methods = DeliverySystem.recommend_delivery(target_tissue)
    
    print(f"目标组织: {target_tissue}")
    print(f"推荐递送方式:")
    for method in methods:
        info = DeliverySystem.DELIVERY_METHODS[method]
        print(f"  - {info['name']}: {info['advantage']}")
    
    dosage = DeliverySystem.calculate_dosage(len(sense), target_tissue)
    print(f"推荐剂量: {dosage:.2f} mg/kg")
    
    # 4. 数据保存演示
    print("\n\n【模块4】数据保存")
    print("-" * 40)
    
    db = SequenceDatabase()
    db.add_sequence(
        name="siRNA_Example",
        sequence=sense,
        description="示例siRNA序列，靶向mRNA的第5-25位"
    )
    
    print("序列已保存到数据库")
    print(f"数据库中的序列: {db.list_sequences()}")
    
    # 保存设计结果
    results = [{
        '名称': 'siRNA_Example',
        '正义链': sense,
        '反义链': antisense,
        'GC含量': f"{gc_content:.2f}%",
        'Tm': f"{tm:.2f}°C",
        '修饰数': len(modifier.get_modifications())
    }]
    
    DataHandler.write_csv(results, 'data/design_results.csv')
    print("设计结果已保存到 data/design_results.csv")
    
    print("\n" + "=" * 60)
    print("            演示完成！")
    print("=" * 60)

if __name__ == "__main__":
    main()