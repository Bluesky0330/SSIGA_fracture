#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
	echo "Usage: $0 <ssh_alias> [interval_seconds]" >&2
	exit 1
fi

HOST="$1"
INTERVAL="${2:-2}"

REMOTE_CMD='ps -eo user:20=,pcpu= --no-headers | awk '\''{cpu[$1]+=$2} END {for (u in cpu) printf "%-20s %.1f%%\n", u, cpu[u]}'\'' | sort -k2 -rn'
HEADER_CMD="ssh -T ${HOST} 'ps -eo user:20=,pcpu= --no-headers | awk ... | sort -k2 -rn'"

trap 'exit 0' INT TERM

clear_previous_block() {
	local lines="$1"
	if [[ "$lines" -gt 0 ]]; then
		printf '\033[%sA\033[J' "$lines"
	fi
}

previous_lines=0

while true; do
	rendered_output=$(ssh -T "$HOST" "$REMOTE_CMD")
	clear_previous_block "$previous_lines"
	printf 'Every %ss: %s\n\n%s\n' "$INTERVAL" "$HEADER_CMD" "$rendered_output"
	previous_lines=$((2 + $(printf '%s\n' "$rendered_output" | wc -l)))
	sleep "$INTERVAL"
done