#!/usr/bin/env bash
# Run the MESS executables on examples/ and compare their output with tests/reference/.
# The reference outputs come from the upstream static binaries.
# Usage: tests/run-examples.sh [RTOL]   (default relative tolerance 1e-6)
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
rtol="${1:-1e-6}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cp -r "$root/examples/." "$work/"
cd "$work"

# Compare two files token by token: numbers within rtol, everything else exactly.
compare() {
  awk -v rtol="$rtol" -v name="$1" '
    function isnum(x) { return x ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/ }
    NR == FNR { ref[FNR] = $0; nref = FNR; next }
    {
      n = split(ref[FNR], a); m = split($0, b)
      if (n != m) { printf "%s:%d: token count differs\n", name, FNR; bad = 1; next }
      for (i = 1; i <= n; i++) {
        if (isnum(a[i]) && isnum(b[i])) {
          d = a[i] - b[i]; if (d < 0) d = -d
          s = (a[i] < 0 ? -a[i] : a[i]); t = (b[i] < 0 ? -b[i] : b[i]); if (t > s) s = t
          if (d > rtol * s) { printf "%s:%d: %s vs %s\n", name, FNR, a[i], b[i]; bad = 1 }
        } else if (a[i] != b[i]) { printf "%s:%d: \"%s\" vs \"%s\"\n", name, FNR, a[i], b[i]; bad = 1 }
      }
    }
    END { if (FNR != nref) { printf "%s: line count differs\n", name; bad = 1 }; exit bad }
  ' "$root/tests/reference/$1" "$1"
}

status=0
run() {  # run <executable> <input> <output>
  echo "== $1 $2"
  if ! "$1" "$2" > "$1.stdout" 2>&1; then
    echo "FAIL: $1 exited with an error"; tail -5 "$1.stdout"; status=1; return
  fi
  if [ ! -s "$3" ]; then echo "FAIL: $3 is missing or empty"; status=1; return; fi
  if compare "$3"; then echo "ok"; else echo "FAIL: $3 differs from tests/reference/$3"; status=1; fi
}

run mess mess.inp mess.out
run messpf messpf.inp messpf.dat
run messabs messabs.inp messabs.out

if [ "$status" -ne 0 ] && [ -n "${MESS_TEST_KEEP:-}" ]; then
  mkdir -p "$MESS_TEST_KEEP" && cp -r "$work/." "$MESS_TEST_KEEP/"
fi
exit "$status"
