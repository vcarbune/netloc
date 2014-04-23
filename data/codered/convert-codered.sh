#!/bash/bin

python convert-codered.py > tmp
sort -k2 -n tmp > codered_gt.txt
rm -rf tmp
