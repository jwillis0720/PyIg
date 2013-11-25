# -*- mode: python -*-
a = Analysis(['../src/PyIg_gui'],
             pathex=['/Users/jordanwillis/utilities/pyigblast/builds'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='PyIg',
          debug=True,
          strip=None,
          upx=True,
          console=False , icon='../splashes/full_antibody.icns')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='PyIg')
app = BUNDLE(coll,
             name='PyIg.app',
             icon='../splashes/full_antibody.icns')
