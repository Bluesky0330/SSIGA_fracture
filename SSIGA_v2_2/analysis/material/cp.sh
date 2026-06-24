#!/bin/bash

OPTION_FILE=$1

# 宛先ディレクトリ（任意引数2 or 環境変数OUTPUT_DIR or ../output_* を自動検出）
DEST_ARG=${2-}

resolve_dest_dir() {
    local dest

    # 1) 明示指定があれば最優先
    if [ -n "$DEST_ARG" ]; then
        if [ -d "$DEST_ARG" ]; then
            echo "$DEST_ARG"
            return 0
        else
            echo "Destination directory $DEST_ARG not found." 1>&2
            return 1
        fi
    fi

    # 2) OUTPUT_DIR 環境変数があれば優先
    if [ -n "${OUTPUT_DIR-}" ] && [ -d "../$OUTPUT_DIR" ]; then
        echo "../$OUTPUT_DIR"
        return 0
    fi

    # 3) ../output_* を検出（mkdir_output後は1つ存在する想定）
    #   複数ある場合は更新時刻が最新のものを選ぶ
    local found_latest
    found_latest=$(
        find .. -maxdepth 1 -type d -name 'output_*' -printf '%T@ %p\n' 2>/dev/null | \
        sort -nr | head -n1 | awk '{print $2}'
    )

    if [ -n "$found_latest" ] && [ -d "$found_latest" ]; then
        echo "$found_latest"
        return 0
    fi

    # 4) 後方互換: ../output があれば使用
    if [ -d "../output" ]; then
        echo "../output"
        return 0
    fi

    echo "No output_* directory found near material/. Please ensure make has created it." 1>&2
    return 1
}

if [ ! -f "$OPTION_FILE" ];then
    echo "Option file $OPTION_FILE not found."
    exit 1
fi

# 宛先ディレクトリ決定
DEST_DIR=$(resolve_dest_dir) || exit 1
echo "Copy destination: $DEST_DIR"

cp "./$OPTION_FILE" "$DEST_DIR"

# CRLFに対応しつつ、コメント行と空行を飛ばして最後の行も含めて処理
while IFS= read -r file || [ -n "$file" ]; do
    # \r (CR) を削除
    file=$(echo "$file" | tr -d '\r')

    # '#' で始まる行と空行をスキップ
    if [[ "$file" =~ ^# ]] || [[ -z "$file" ]]; then
        continue
    fi

    # コピー
    if [ -f "./$file" ]; then
        cp "./$file" "$DEST_DIR"
    elif [ -f "../input/$file" ]; then
        cp "../input/$file" "$DEST_DIR"
    else
        echo "Material file $file not found in ./ or ../input/"
        exit 1
    fi
done < "$OPTION_FILE"
