mkdir -p src

make fast.x

cd src
rsync -avz --quiet "../../src/"*.cpp "./"
rsync -avz --quiet "../../src/"*.hpp "./"

rsync -avz --quiet "../../analysis/material/constant.ini" "./"
rm -f merged_output.cpp SIGA.x ../SIGA.x
make merged_output.cpp
make SIGA.x
cd ../

rsync -avz --quiet "./src/SIGA.x" ./
echo "build completed."