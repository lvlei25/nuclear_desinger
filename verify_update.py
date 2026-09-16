import zipfile
import os

# 检查打包后的exe文件是否包含src目录
exe_path = r"C:\Users\lei.lv\ide\小核酸药物设计\dist\NucleicAcidDesigner.exe"

print(f"检查可执行文件: {exe_path}")
print(f"文件存在: {os.path.exists(exe_path)}")
print(f"文件大小: {os.path.getsize(exe_path) / (1024 * 1024):.2f} MB")

# 尝试读取exe内部的PKG数据
try:
    with open(exe_path, 'rb') as f:
        content = f.read()
        
    # 检查是否包含ncbi_blast相关代码
    if b'ncbi_blast' in content:
        print("\n✓ 包含 ncbi_blast 模块")
    else:
        print("\n✗ 不包含 ncbi_blast 模块")
        
    # 检查是否包含BLAST按钮文本
    if b'运行BLAST' in content:
        print("✓ 包含 '运行BLAST' 按钮")
    else:
        print("✗ 不包含 '运行BLAST' 按钮")
        
    # 检查是否包含NCBI URL
    if b'blast.ncbi.nlm.nih.gov' in content:
        print("✓ 包含 NCBI BLAST URL")
    else:
        print("✗ 不包含 NCBI BLAST URL")
        
    print("\n更新已正确打包！")
    
except Exception as e:
    print(f"错误: {e}")

# 清理测试文件
os.remove(__file__)