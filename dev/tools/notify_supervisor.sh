#!/usr/bin/env bash
# Usage from any lane: bash /home/datahouse-raid5/wjt/COSMolKit/dev/tools/notify_supervisor.sh "report"
# Success means submitted to the supervisor's conversation/queue, not reviewed.
set -euo pipefail

if [[ $# != 1 || -z ${1//[[:space:]]/} ]]; then
    echo 'Usage: notify_supervisor.sh "nonempty report"' >&2
    exit 2
fi
if (( ${#1} > 1000 )); then
    echo 'Notification aborted: use at most 1000 characters; put detailed evidence in the receipt.' >&2
    exit 2
fi

exec 9>/tmp/ck-supervisor-notify.lock
flock -w 10 9 || { echo 'Notification failed: sender lock busy.' >&2; exit 1; }

screen() {
    zellij --session ck action dump-screen --pane-id terminal_0
}

composer_empty() {
    awk '
        /^[[:space:]]*›/ {
            last = $0
            sub(/^[[:space:]]*/, "", last)
            sub(/[[:space:]]*$/, "", last)
        }
        END { exit last == "› Ask Codex to do anything" ? 0 : 1 }
    ' <<< "$1"
}

before=$(screen)
if ! composer_empty "$before"; then
    echo 'Notification aborted: supervisor composer is not empty; nothing pasted.' >&2
    exit 1
fi

# A unique prefix distinguishes this submission from earlier reports and drafts.
marker="[CK-NOTIFY:$(date +%s%N)-$$]"
zellij --session ck action write-chars --pane-id terminal_0 "$marker $1"
# Zellij returns before the UI necessarily consumes the complete paste.
# Wait BEFORE Enter; otherwise it can arrive in the middle of the payload.
sleep 4
zellij --session ck action write --pane-id terminal_0 13
sleep 2
# A separate Enter is required when the first merely ends a pasted draft.
zellij --session ck action write --pane-id terminal_0 13

extra_enter=false
for attempt in 1 2 3 4 5 6; do
    sleep 0.5
    after=$(screen)
    if composer_empty "$after" && [[ $after == *"$marker"* ]]; then
        printf 'DELIVERED %s\n' "$marker"
        exit 0
    fi
    # Scroll/paste handling can consume an Enter. Retry Enter once, ONLY when
    # the last composer prompt still starts with our own unique marker.
    if [[ $extra_enter == false ]] && awk -v prefix="› $marker " '
        /^[[:space:]]*›/ { last = $0; sub(/^[[:space:]]*/, "", last) }
        END { exit index(last, prefix) == 1 ? 0 : 1 }
    ' <<< "$after"; then
        zellij --session ck action write --pane-id terminal_0 13
        extra_enter=true
    fi
done

echo "Notification unverified: $marker not confirmed submitted with an empty composer. Inspect terminal_0; do not blindly resend." >&2
exit 1
