# -*- mode: python ; coding: utf-8 -*-

# 打包配置（console 调试版）：在项目根目录执行 pyinstaller spec/NucleicAcidDesigner.spec
# 注意：spec 文件已移入 spec/ 目录，内部相对路径均指向上一级（项目根目录）
a = Analysis(
    ['../gui.py'],
    pathex=['..'],
    binaries=[],
    datas=[],
    hiddenimports=['src.ncbi_blast', 'src.sequence_design', 'src.chemical_modification', 'src.data_processing', 'requests'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='NucleicAcidDesigner',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)