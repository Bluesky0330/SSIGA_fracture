インプットデータのフォーマットが以前のバージョンと違うので注意してください．
一行目にDIMENSIONを入力する箇所が変更されています．
MI3D -> MC3D の順で操作すると問題なく作成されるはずです．
詳しくはmake_input/make_input_3D/MI3D, MC3D の textfile/template...を参照してください．

./src/S_IGA_header.hpp で主要な解析条件を変更できるようになっています．適宜変更してください．

・コンパイル方法
$ cd (Makefileがあるディレクトリ)
$ make
(or $ make compile)

・実行方法
Makefile 内の run等 のコマンドライン引数のファイル名や数を適宜変更
$ cd (Makefileがあるディレクトリ)
$ make run

・Makefileを使用しない場合のコンパイル・実行方法
$ cd (path/to/src)
$ g++ -o S-IGA_PCG.x S-IGA_PCG_main_v2.cpp S-IGA_PCG_sub_pre_v2.cpp S-IGA_PCG_sub_solver_v2.cpp S-IGA_PCG_sub_tool_v2.cpp S-IGA_PCG_sub_post_v2.cpp S-IGA_PCG_sub_NR_v2.cpp -Ofast -lm -Wall -Wextra -pipe -std=c++17
$ cd (path/to/output_dir)
$ mkdir bin
IGAの場合
$ (path/to/S-IGA_PCG.x) (path/to/input.txt)
S-IGAの場合
$ (path/to/S-IGA_PCG.x) (path/to/input_global_patch.txt) (path/to/input_local_patch.txt)
括弧内のパスは適宜指定してください．
解析を実行する前にoutputファイル内にbinというファイルを作成しないと，
paraview用の出力ファイルのバイナリデータが作成されないので注意してください．

・データ確認方法
outputフォルダ内のdatファイルで確認してください．
datファイルの出力に関する設定は./src/S_IGA_header.hppで行うことができます．

出力されるファイル
_Displacement_overlay_at_ele_vertex.dat
_Displacement_overlay_at_GP.dat
_PhysicalCoordinate_overlay_at_ele_vertex.dat
_PhysicalCoordinate_overlay_at_GP.dat
_Strain_overlay_at_ele_vertex.dat
_Strain_overlay_at_GP.dat
_Stress_overlay_at_ele_vertex.dat
_Stress_overlay_at_GP.dat

local要素は重ね合わせを行った値，global要素はローカル領域でも重ね合わせを行っていない値となっています．

~at_ele_vertex.dat:   要素の各頂点での値
~at_GP.dat:           要素のガウス点での値

Displacement:       変位
PhysicalCoordinate: 物理座標
Strain:             ひずみ
Stress:             応力

・可視化方法
global_patch.xmf
global_patch_boundary_line.xmf
global_patch_control_point.xmf
local_patch.xmf
local_patch_boundary_line.xmf
local_patch_control_point.xmf
deformed_global_patch.xmf
deformed_global_patch_boundary_line.xmf
deformed_local_patch.xmf
deformed_local_patch_boundary_line.xmf
(localはS-IGAの場合のみ, deformedは変形後のメッシュ出力を有効にした場合のみ)

上記のファイルをparaviewのアプリケーションにドラッグしてXDMF Readerを選択すると可視化できます．
paraview:   https://www.paraview.org/download/

paraview の bin (例 C:\Program Files\ParaView 5.10.1-Windows-Python3.9-msvc2017-AMD64\bin)を環境変数に追加し，
パスを通すと，以下のコマンドで可視化できます．

IGA の場合
$ cd (Makefileがあるディレクトリ)
$ make view_0

S-IGA の場合
$ cd (Makefileがあるディレクトリ)
$ make view_1

変形後の出力, IGA の場合
$ cd (Makefileがあるディレクトリ)
$ make view_2

変形後の出力, S-IGA の場合
$ cd (Makefileがあるディレクトリ)
$ make view_3