mkdir -p src

make -j$(nproc) fast.x

cd src
rsync -avz --quiet "../../src/"*.cpp "./"
rsync -avz --quiet "../../src/"*.hpp "./"

rsync -avz --quiet "../../analysis/material/constant.ini" "./"
rm -f merged_output.cpp SIGA.x ../SIGA.x
make -j$(nproc) merged_output.cpp
make -j$(nproc) SIGA.x
cd ../

rsync -avz --quiet "./src/SIGA.x" ./
echo "build completed."