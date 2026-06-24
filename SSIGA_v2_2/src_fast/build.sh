mkdir -p src

make fast.x

cd src
rsync -avz --quiet "../../src/"*.cpp "./"
rsync -avz --quiet "../../src/"*.hpp "./"

rsync -avz --quiet "../../analysis/material/constant.ini" "./"
rm -f merged_output.cpp SSIGA.x ../SSIGA.x
make merged_output.cpp
make SSIGA.x
cd ../

rsync -avz --quiet "./src/SSIGA.x" ./
echo "build completed."