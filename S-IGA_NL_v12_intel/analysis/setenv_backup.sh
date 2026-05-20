#!/bin/bash

# 環境内のjemallocを探してexportし，見つからなければエラーを出力する
JEMALLOC_PATH=$(find /usr/lib* -name "libjemalloc.so*" 2>/dev/null | head -n 1)
if [ -n "$JEMALLOC_PATH" ]; then
    export LD_PRELOAD="$JEMALLOC_PATH"
else
    echo "jemalloc not found"
    echo "Please install jemalloc and try again."
    exit 1
fi

# スレッド数，バインディング，配置を設定
# export OMP_NUM_THREADS=$(nproc)
export OMP_NUM_THREADS=40
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# OpenMPのKMP_AFFINITYとKMP_BLOCKTIMEを設定
# export KMP_AFFINITY=granularity=fine,compact,1,0 # OMP_PROC_BINDとOMP_PLACESが無視される
export KMP_BLOCKTIME=0

# このシェルで設定したすべての環境変数を色を付けて出力
printf "\n" >&2
printf "\tEnvironment variables set by setenv.sh:\n" >&2
printf "\t\tLD_PRELOAD\t= \e[35m%s\e[0m\n" "$LD_PRELOAD" >&2
printf "\t\tOMP_NUM_THREADS\t= \e[35m%s\e[0m\n" "$OMP_NUM_THREADS" >&2
printf "\t\tOMP_PROC_BIND\t= \e[35m%s\e[0m\n" "$OMP_PROC_BIND" >&2
printf "\t\tOMP_PLACES\t= \e[35m%s\e[0m\n" "$OMP_PLACES" >&2
printf "\t\tKMP_AFFINITY\t= \e[35m%s\e[0m\n" "$KMP_AFFINITY" >&2
printf "\t\tKMP_BLOCKTIME\t= \e[35m%s\e[0m\n" "$KMP_BLOCKTIME" >&2
printf "\n" >&2