# -*- mode: python -*-
def Datafiles(*filenames, **kw):
    import os
    
    def datafile(path, strip_path=True):
        parts = path.split('/')
        path = name = os.path.join(*parts)
        if strip_path:
            name = os.path.basename(path)
        return name, path, 'DATA'

    strip_path = kw.get('strip_path', True)
    return TOC(
        datafile(filename, strip_path=strip_path)
        for filename in filenames
        if os.path.isfile(filename))

docfiles = Datafiles(r'README_gui.txt',strip_path=True)


a = Analysis(['PyIg_gui.py'],
             pathex=['C:\\Users\\crowelab\\pyigblast\\src'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='test.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False , icon='splashes\\full_antibody.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               docfiles,
               strip=None,
               upx=True,
               name='test')
