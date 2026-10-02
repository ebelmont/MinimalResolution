#!/usr/bin/env python3
"""Summarize mr_BP's DIAG_TIMING output as a tree of steps.

Usage:
    python3 diag_summary.py ERRFILE [--depth N] [--min-share P]
    python3 diag_summary.py ERRFILE --per-s PATH

ERRFILE is mr_BP's stderr from a run with DIAG_TIMING=1 set (see diag.h).
Each timed step is shown with its total wall time, its share of the
enclosing step, cpu time, utilization (cpu / (wall * threads), so a serial
step on 8 threads shows ~12%), how many times it ran, and the peak resident
memory seen when it finished. Per-homological-degree steps ("s=12") are
merged into one "s=*" row; --per-s PATH lists them one degree at a time,
e.g. --per-s mr_BP/resolve/s=*/resolvor_modeled/del_and_gaussian.
"(other)" rows are time in a step not covered by any timed sub-step.
Works on a run still in progress: steps begun but not finished are listed
at the end.
"""
import argparse
import re
import sys
from collections import defaultdict

LINE = re.compile(r"^DIAG path=(\S+) (.*)$")
BEGIN = re.compile(r"^DIAGBEGIN path=(\S+) t=(\S+)")
S_RE = re.compile(r"s=\d+")


def fmt_time(x):
    if x >= 3600:
        return f"{x/3600:.1f}h"
    if x >= 60:
        return f"{x/60:.1f}m"
    return f"{x:.1f}s"


def parse(fn):
    done = []
    open_steps = {}
    last_t = 0.0
    with open(fn, errors="replace") as f:
        for line in f:
            m = BEGIN.match(line)
            if m:
                open_steps[m.group(1)] = float(m.group(2))
                last_t = max(last_t, float(m.group(2)))
                continue
            m = LINE.match(line)
            if not m:
                continue
            path = m.group(1)
            kv = dict(x.split("=", 1) for x in m.group(2).split())
            rec = {k: float(v) for k, v in kv.items()}
            rec["path"] = path
            done.append(rec)
            open_steps.pop(path, None)
            last_t = max(last_t, rec.get("t", 0.0))
    return done, open_steps, last_t


def aggregate(done):
    agg = defaultdict(lambda: {"wall": 0.0, "cpu": 0.0, "threads": 0, "n": 0, "peak": 0.0})
    for r in done:
        a = agg[S_RE.sub("s=*", r["path"])]
        a["wall"] += r["wall"]
        a["cpu"] += r["cpu"]
        a["threads"] = max(a["threads"], int(r["threads"]))
        a["n"] += 1
        a["peak"] = max(a["peak"], r.get("peak_gb", 0.0))
    return agg


def print_tree(agg, depth, min_share):
    children = defaultdict(list)
    for p in agg:
        parent = p.rsplit("/", 1)[0] if "/" in p else ""
        children[parent].append(p)
    # a step whose enclosing step hasn't finished yet (run still going)
    # is shown at top level
    roots = [p for p in agg if (p.rsplit("/", 1)[0] if "/" in p else "") not in agg]
    total = sum(agg[r]["wall"] for r in roots) or 1.0
    print(f"{'step':<60} {'wall':>8} {'%parent':>8} {'cpu':>8} {'util':>6} {'count':>6} {'peakGB':>7}")

    def row(name, indent, a, parent_wall):
        share = 100 * a["wall"] / parent_wall if parent_wall else 100.0
        util = a["cpu"] / (a["wall"] * a["threads"]) if a["wall"] > 0 and a["threads"] else 0
        print(f"{'  '*indent + name:<60} {fmt_time(a['wall']):>8} {share:>7.1f}% "
              f"{fmt_time(a['cpu']):>8} {100*util:>5.0f}% {a['n']:>6} {a['peak']:>7.1f}")

    def walk(p, indent, parent_wall):
        a = agg[p]
        row(p.rsplit("/", 1)[-1] if indent else p, indent, a, parent_wall)
        if indent + 1 >= depth:
            return
        kids = sorted(children.get(p, []), key=lambda k: -agg[k]["wall"])
        covered = 0.0
        for k in kids:
            covered += agg[k]["wall"]
            if agg[k]["wall"] >= min_share / 100 * total:
                walk(k, indent + 1, a["wall"])
        other = a["wall"] - covered
        if kids and other >= min_share / 100 * total:
            row("(other)", indent + 1, {"wall": other, "cpu": 0, "threads": 0, "n": 0, "peak": 0}, a["wall"])

    for r in sorted(roots, key=lambda k: -agg[k]["wall"]):
        walk(r, 0, 0)


def print_per_s(done, pattern):
    rx = re.compile("^" + re.escape(pattern).replace(r"s=\*", r"s=(\d+)") + "$")
    rows = []
    for r in done:
        m = rx.match(r["path"])
        if m:
            rows.append((int(m.group(1)), r))
    if not rows:
        sys.exit(f"no lines match {pattern}")
    total = sum(r["wall"] for _, r in rows) or 1.0
    print(f"{'s':>4} {'wall':>8} {'%total':>7} {'cpu':>8} {'util':>6} {'rssGB':>7} {'peakGB':>7}")
    for s, r in sorted(rows, key=lambda x: x[0]):
        util = r["cpu"] / (r["wall"] * r["threads"]) if r["wall"] > 0 else 0
        print(f"{s:>4} {fmt_time(r['wall']):>8} {100*r['wall']/total:>6.1f}% {fmt_time(r['cpu']):>8} "
              f"{100*util:>5.0f}% {r.get('rss_gb', 0):>7.1f} {r.get('peak_gb', 0):>7.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("errfile")
    ap.add_argument("--depth", type=int, default=5, help="tree depth to show (default 5)")
    ap.add_argument("--min-share", type=float, default=0.5,
                    help="hide steps under this %% of total wall time (default 0.5)")
    ap.add_argument("--per-s", metavar="PATH", help="list one step degree by degree")
    args = ap.parse_args()
    done, open_steps, last_t = parse(args.errfile)
    if args.per_s:
        print_per_s(done, args.per_s)
        return
    # leave out sub-steps of a degree (s=N) that hasn't finished yet, so the
    # totals only cover finished degrees
    open_s = [p for p in open_steps if re.search(r"s=\d+$", p)]
    done_full = [r for r in done if not any(r["path"].startswith(p + "/") for p in open_s)]
    print_tree(aggregate(done_full), args.depth, args.min_share)
    if open_steps:
        print(f"\nstill running (last timestamp t={fmt_time(last_t)}):")
        for p, t in sorted(open_steps.items(), key=lambda x: x[1]):
            print(f"  {p}  (started t={fmt_time(t)}, {fmt_time(last_t - t)} ago at last timestamp)")


if __name__ == "__main__":
    main()
