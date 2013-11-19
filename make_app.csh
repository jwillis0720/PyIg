rm -rf /Applications/PyIg*
python ../pyinstaller-dev/pyinstaller.py --log-level=DEBUG --debug --distpath=/Applications/ -n PyIg --clean --window --workpath=/tmp/ -i splashes/full_antibody.icns gui_setup.py
rm -rf tmp/
cp -r junctional_data /Applications/PyIg.app/Contents/MacOS/.
cp -r database /Applications/PyIg.app/Contents/MacOS/.
cp -r internal_data /Applications/PyIg.app/Contents/MacOS/.
cp -r optional_file /Applications/PyIg.app/Contents/MacOS/.
cp -r igblast_source/darwin/ncbi-igblast-1.2.0/bin/igblastn /Applications/PyIg.app/Contents/MacOS/.
cp README_gui.txt /Applications/PyIg.app/Contents/MacOS/
rm -rf /Applications/PyIg
