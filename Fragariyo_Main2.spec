# -*- mode: python -*-
# ('TerminalFragments_template_example.csv', 'InternalFragsments_template_example.csv', 'ModificationsRepository.txt','.')

import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)

block_cipher = None

b = [
    ('C:\\Users\\caror\\AppData\\Roaming\\Python\\Python37\\site-packages\\sklearn\\.libs\\vcomp140.dll', '.\\sklearn\\.libs')
    ]
	
resource_files = [
('UI/*', 'UI'),
('tooltips.txt', '.'),
('C:\\Users\\caror\\AppData\\Roaming\\Python\\Python37\\site-packages\\PyQt5\\Qt5\\bin\\*','PyQt5\\Qt5\\bin')]

a = Analysis(['Fragariyo_Main.py'],
             pathex=['C:\\Users\\caror\\PycharmProjects\\fragmentor_build', 'C:\\Users\\caror\\AppData\\Roaming\\Python\\Python37\\site-packages\\PyQt5\\Qt5\bin'],
             binaries=b,
             datas=resource_files,
             hiddenimports=[],
             hookspath=['.\\extra_hooks'],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='Fragariyo',
          debug=False,
          strip=False,
          upx=True,
          icon='Fragariyo_win10.ico',
          console=True)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='Fragariyo')
			   

