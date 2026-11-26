#!/usr/bin/env bash
# Decide whether an assembler may resume from the files already in its work directory.
#
# usage: resume_guard.sh WORKDIR READS [PARAMS...]
#
# Prints "resume" or "fresh" on stdout. The run's signature (reads path, size and mtime plus
# PARAMS, i.e. the assembler options that change the result; never threads or memory) is
# kept in WORKDIR/.ouroassemble_resume:
#   - WORKDIR empty or missing          -> fresh
#   - signature unchanged               -> resume (e.g. the job was killed at walltime)
#   - signature changed                 -> WORKDIR is emptied, fresh (stale state discarded)
#   - files but no signature (a run started before resume support) -> adopted, resume
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: resume_guard.sh WORKDIR READS [PARAMS...]" >&2
    exit 2
fi
workdir=$1
reads=$2
shift 2
case $workdir in
    "" | "/" | "." | "..") echo "resume_guard: refusing work directory '$workdir'" >&2; exit 2 ;;
esac

stats=$(stat -L -c 'size=%s mtime=%Y' "$reads" 2>/dev/null || stat -L -f 'size=%z mtime=%m' "$reads")
sig="reads=$(cd "$(dirname "$reads")" && pwd -P)/$(basename "$reads") $stats params=$*"
sigfile=$workdir/.ouroassemble_resume

mkdir -p "$workdir"
if [[ -f $sigfile ]]; then
    if [[ $(cat "$sigfile") == "$sig" ]]; then
        echo resume
        exit 0
    fi
    echo "resume_guard: reads or parameters changed; clearing $workdir" >&2
    find "$workdir" -mindepth 1 -delete
elif [[ -n $(ls -A "$workdir") ]]; then
    printf '%s\n' "$sig" > "$sigfile"
    echo resume
    exit 0
fi
printf '%s\n' "$sig" > "$sigfile"
echo fresh
