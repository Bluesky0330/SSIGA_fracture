#!/bin/bash
set -euo pipefail

# --- 色・装飾定義 ---
CYAN='\033[36m'
YELLOW='\033[33m'
GREEN='\033[32m'
RED='\033[1;31m'
RESET='\033[0m'

# ローカル端末(Cygwin)側もUTF-8で扱う
export LANG=ja_JP.UTF-8
export LC_ALL=ja_JP.UTF-8

# プレフィックス定義
P_LOC="${CYAN}[Local]${RESET} "
P_REM="${CYAN}[Remote]${RESET}"

# --- ヘッダー関数 ---
print_section() {
    echo -e "\n${CYAN}============================================================${RESET}"
    echo -e "${CYAN}   $1${RESET}"
    echo -e "${CYAN}============================================================${RESET}"
}
print_sub() {
    echo -e "${YELLOW}--- $1 ---${RESET}"
}

# --- 設定 ---
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
DIR_NAME=$(basename "$SCRIPT_DIR")
REMOTE_DIR="~/$DIR_NAME"
REMOTE_LOG="$REMOTE_DIR/remote_log.txt"
REMOTE_LOG_COLOR="$REMOTE_DIR/remote_log_color.txt"
REMOTE_ENV="LANG=ja_JP.UTF-8 LC_ALL=ja_JP.UTF-8"

# --- 引数処理 (getoptsでオプション解析) ---
FORCE_RUN=false
ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--force)
            FORCE_RUN=true
            shift
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

# 引数を再割り当て
set -- "${ARGS[@]}"
HOST="${REMOTE_HOST:-}"
TARGET=""

if [ $# -eq 2 ]; then 
    HOST="$1"
    TARGET="$2"
elif [ $# -eq 1 ]; then 
    TARGET="$1"
fi

[ -z "$HOST" ] && { echo "Usage: $0 [user@host] <make_target> [-f|--force]"; exit 1; }

# tmuxセッション名にはピリオド(.)が使えないため、_に置換した変数を用意
SAFE_TARGET="${TARGET//./_}"
SESSION="build_${DIR_NAME}_${SAFE_TARGET}"

# --- 1. リモートの状態確認 ---
print_section "状態確認"

# 誤検知除外(psコマンド)でPID取得
CHECK_CMD="ps -ef | grep 'make $TARGET' | grep -v grep | grep -v bash | grep -v ssh | awk '{print \$2}' | head -n 1"
REMOTE_PID=$(ssh "$HOST" "$REMOTE_ENV $CHECK_CMD")

MODE=""
if [ -n "$REMOTE_PID" ]; then
    # A. 実行中 -> 強制フラグがあっても、実行中はアタッチを優先（安全のため）
    MODE="ATTACH"
    echo -e "${P_REM} 実行中のプロセスを検知 (PID: $REMOTE_PID)"
    if [ "$FORCE_RUN" = "true" ]; then
        echo -e "${YELLOW}※ -f が指定されましたが、実行中のためアタッチします。再実行したい場合はリモートで停止してください。${RESET}"
    fi
else
    # B. 停止中 -> 未回収チェック or 再実行チェック

    # 強制実行フラグがあれば、未回収チェックなどをスキップしてRERUNへ
    if [ "$FORCE_RUN" = "true" ]; then
        MODE="RERUN"
        echo -e "${P_LOC} ${RED}強制再実行モード (-f)${RESET}"
    else
        # 通常フロー: 未回収ファイルの確認
        UNSYNCED_COUNT=$(rsync -nai "$HOST:$REMOTE_DIR/analysis/output_$TARGET/" "$SCRIPT_DIR/analysis/output_$TARGET/" | grep '^>f' | wc -l || true)

        if [ "$UNSYNCED_COUNT" -gt 0 ]; then
            print_sub "未回収データの検出"
            echo -e "${YELLOW}解析は既に終了していますが、未回収の結果 ($UNSYNCED_COUNT files) があります。${RESET}"
            echo -e "${P_LOC} 上書き防止のため、結果を回収して終了します。"
            echo -e "${P_LOC} (再解析したい場合は回収後に再度実行するか、-f オプションを付けてください)"
            
            echo -e "${P_LOC} ダウンロード中..."
            mkdir -p "$SCRIPT_DIR/analysis/output_$TARGET"
            rsync -azP "$HOST:$REMOTE_DIR/analysis/output_$TARGET/" "$SCRIPT_DIR/analysis/output_$TARGET/" 2>/dev/null || true
            rsync -azP "$HOST:$REMOTE_DIR/remote_log.txt" "$SCRIPT_DIR/analysis/output_$TARGET/" 2>/dev/null || true
            
            echo -e "${GREEN}===== 回収完了 =====${RESET}"
            exit 0
        else
            # 未回収なし -> ソース変更チェック
            CHANGE_COUNT=$(rsync -nai --exclude='*output*' --exclude='.git' "$SCRIPT_DIR/" "$HOST:$REMOTE_DIR/" | grep '^<' | wc -l || true)
            
            if [ "$CHANGE_COUNT" -gt 0 ]; then
                MODE="RERUN"
                echo -e "${P_LOC} ソース変更あり ($CHANGE_COUNT files) -> 再実行"
            else
                MODE="SYNC"
                echo -e "${P_LOC} 変更なし -> 結果回収のみ"
            fi
        fi
    fi
fi

# --- 2. 分岐実行 (ATTACH / RERUN) ---

if [ -n "$TARGET" ] && [ "$MODE" != "RERUN" ]; then
    echo -e "${P_REM} makeターゲットの存在チェック (SYNC/ATTACH): $TARGET"
    if ! ssh "$HOST" "$REMOTE_ENV cd $REMOTE_DIR && make -n $TARGET >/dev/null 2>&1"; then
        echo -e "${RED}${P_REM} Error: make ターゲット '$TARGET' がリモートに存在しないか make 実行時にエラーが発生しました。中止します。${RESET}"
        ssh "$HOST" "cat $REMOTE_LOG" || true
        exit 1
    fi
fi

if [ "$MODE" = "RERUN" ]; then
    print_section "解析再実行 (RERUN)"

    print_sub "ファイル同期"
    echo -e "${P_LOC} rsync実行中..."
    rsync -azP --exclude='*output*' "$SCRIPT_DIR/" "$HOST:$REMOTE_DIR/"

    print_sub "リモート実行準備"
    echo -e "${P_REM} ログ初期化 & 古いセッションの削除"
    # 重要: 古いセッション(PIDなし)が残っているとnew-sessionが失敗するため、確実にkill
    ssh "$HOST" "$REMOTE_ENV mkdir -p $REMOTE_DIR && rm -f $REMOTE_LOG $REMOTE_LOG_COLOR && printf '\xEF\xBB\xBF' > $REMOTE_LOG && : > $REMOTE_LOG_COLOR && tmux kill-session -t $SESSION 2>/dev/null || true"

    echo -e "${P_REM} makeターゲットの存在チェック: $TARGET"
    if ! ssh "$HOST" "$REMOTE_ENV cd $REMOTE_DIR && make -n $TARGET >/dev/null 2>&1"; then
        echo -e "${RED}${P_REM} Error: make ターゲット '$TARGET' が見つからないか make 実行時にエラーが発生しました。中止します。${RESET}"
        ssh "$HOST" "cat $REMOTE_LOG" || true
        exit 1
    fi

    # Makeコマンド — bash -lc で実行し、端末はカラーのまま、ログはANSIシーケンス除去して保存
    MAKE_CMD="bash -lc 'cd $REMOTE_DIR && $REMOTE_ENV TERM=xterm-256color CLICOLOR_FORCE=1 FORCE_COLOR=1 stdbuf -oL -eL make $TARGET 2>&1 | tee -a $REMOTE_LOG_COLOR >(stdbuf -oL sed -r '\''s/\x1B\[[0-9;]*[A-Za-z]//g'\'' >> $REMOTE_LOG)'"
    
    echo -e "${P_REM} tmuxセッション作成 & make投入"
    ssh "$HOST" "$REMOTE_ENV tmux new-session -d -s $SESSION 2>/dev/null || true"
    ssh "$HOST" "tmux send-keys -t $SESSION \"$MAKE_CMD\" C-m"

    # 起動待機
    for i in {1..5}; do
        sleep 1
        REMOTE_PID=$(ssh "$HOST" "$REMOTE_ENV $CHECK_CMD")
        [ -n "$REMOTE_PID" ] && break
    done
    
    if [ -n "$REMOTE_PID" ]; then
        echo -e "${GREEN}${P_REM} 解析開始 (PID: $REMOTE_PID)${RESET}"
    else
        echo -e "${RED}${P_REM} Error: 解析即時終了。ログを表示します:${RESET}"
        ssh "$HOST" "cat $REMOTE_LOG"
        exit 1
    fi
fi

# --- 3. ログ監視 (ATTACH or RERUN) ---
if [ "$MODE" != "SYNC" ] && [ -n "$REMOTE_PID" ]; then
    print_section "ログモニタリング"
    echo -e "  - ${RED}'q' キー${RESET} : デタッチ"
    echo -e "  - ${RED}'x' キー${RESET} : Ctrl+C (SIGINT) 相当 — リモート解析を中断"
    echo ""

    ssh "$HOST" "$REMOTE_ENV stdbuf -oL tail -n +1 -f $REMOTE_LOG_COLOR --pid=$REMOTE_PID" < /dev/null &
    SSH_PID=$!
    trap 'kill $SSH_PID 2>/dev/null' EXIT

    while kill -0 $SSH_PID 2>/dev/null; do
        read -t 0.5 -n 1 key || true
        if [[ $key == "q" ]]; then
            echo -e "\n${YELLOW}=== デタッチしました ===${RESET}"
            exit 0
        fi

        if [[ $key == "x" || $key == "X" ]]; then
            echo -e "\n${YELLOW}=== SIGINT をリモートプロセスに送信します (Ctrl+C 相当) ===${RESET}"
            # tmux セッション内のプロセスに Ctrl-C を送信
            ssh "$HOST" "$REMOTE_ENV tmux send-keys -t $SESSION C-c" >/dev/null 2>&1 || true
            echo -e "${YELLOW}リモートプロセスに SIGINT を送信しました。${RESET}"
            # ローカルの tail プロセスを終了して、ログ監視ループを抜ける
            kill $SSH_PID 2>/dev/null || true
            break
        fi
    done
fi

# --- 4. 結果回収 ---
if [ -z "$REMOTE_PID" ] && [ "$MODE" != "SYNC" ]; then
    exit 0
fi

print_section "結果回収"
print_sub "Remote -> Local 同期"

mkdir -p "$SCRIPT_DIR/analysis/output_$TARGET"

echo -e "${P_LOC} rsync実行中..."
rsync -azP "$HOST:$REMOTE_DIR/analysis/output_$TARGET/" "$SCRIPT_DIR/analysis/output_$TARGET/" 2>/dev/null || true
rsync -azP "$HOST:$REMOTE_DIR/remote_log.txt" "$SCRIPT_DIR/analysis/output_$TARGET/" 2>/dev/null || true
rsync -azP "$HOST:$REMOTE_DIR/remote_log_color.txt" "$SCRIPT_DIR/analysis/output_$TARGET/" 2>/dev/null || true

echo -e "\n${GREEN}===== 全工程完了 =====${RESET}"