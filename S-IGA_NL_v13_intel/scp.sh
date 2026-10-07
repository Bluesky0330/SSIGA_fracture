#!/bin/bash

# Convert script file to LF if necessary
[ "$(head -n 1 "$0" | grep -o $'\r')" ] && exec 2> /dev/null < <(tr -d '\r' < "$0") || true

# このスクリプトがあるディレクトリを取得
script_dir=$(dirname "$(realpath "$0")")

# ターゲットディレクトリをスクリプトのあるディレクトリに設定
target_dir_path="$script_dir"

# ターゲットディレクトリ名を取得
target_dir_name=$(basename "$target_dir_path")

# 引数チェック
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <hostname> <operation_number>"
    exit 1
fi

# 引数を変数に格納
HOSTNAME=$1
OPERATION_NUMBER=$2

# 操作番号に基づいてSCPでファイルを転送
case $OPERATION_NUMBER in
    1)
        echo "Operation 1: Copying all files to $HOSTNAME:~/"
        rsync -avzP --exclude='*output*' "$target_dir_path" "$HOSTNAME:~/"
        ;;
    2)
        echo "Operation 2: Copying src to $HOSTNAME:~/$target_dir_name"
        rsync -avzP "$target_dir_path/src" "$HOSTNAME:~/$target_dir_name"
        rsync -avzP "$target_dir_path/analysis/input/"*.txt "$HOSTNAME:~/$target_dir_name/analysis/input"
        rsync -avzP "$target_dir_path/analysis/material/"*.txt "$HOSTNAME:~/$target_dir_name/analysis/material"
        rsync -avzP "$target_dir_path/analysis/material/"*.ini "$HOSTNAME:~/$target_dir_name/analysis/material"
        rsync -avzP "$target_dir_path/analysis/material/"*.sh "$HOSTNAME:~/$target_dir_name/analysis/material"
        rsync -avzP "$target_dir_path/Makefile" "$HOSTNAME:~/$target_dir_name"
        ;;
    3)
        echo "Operation 3: Copying output"
        rm -rf "$target_dir_path/analysis/output"
        mkdir -p "$target_dir_path/analysis/output/bin"
        rsync -avzP "$HOSTNAME:~/$target_dir_name/analysis/output" "$target_dir_path/analysis"
        ;;
    *)
        echo "Invalid operation number. Please use 1, 2, or 3."
        exit 1
        ;;
esac

echo "Operation completed."
