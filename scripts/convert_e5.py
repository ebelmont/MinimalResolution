#!/usr/bin/env python3
"""Turn a directory of `<prefix>_BP*.txt` ANSS tables (e.g. `../E5_200`,
`../380_results`) into E5-style CSVs (factor1,factor2,result / s,f,rank /
factor1,factor2 failures) -- the format the ANSS_prime3 propagator in
`../yuqin` reads directly.

Usage:
    python3 convert_e5.py DIR [OUT_DIR] (--full-ext | --ext-mod-3-40)

DIR must contain files named `<prefix>_BPAANSS_table.txt`, `..._h0.txt`,
`..._a0.txt`, `..._BPB2A_table.txt`, `..._BPBocSS_table.txt`,
`..._BPBocSS_a0.txt`, and `..._BPBocSS_theta{2..7}.txt`. The prefix is
inferred automatically. Output CSVs are written to OUT_DIR (default: DIR
itself).

This is a memory-conscious reimplementation of beta.py's ANSSData, needed
because the theta tables can be tens of millions of lines: they are parsed
*filtered*, keeping only keys that can possibly be looked up (the image of
B2A_inv, bounded by the tiny B2A table), instead of materializing the full
tables. `three()`/`boc_h0table` are not computed (unused by the alpha1 +
beta-family product outputs generated here).

Ext(BP_*/3^40) input contains two copies of (almost) every Ext(BP_*) class:
a "first copy" (image of Ext(BP_*) -> Ext(BP_*/3^40)) and a "second copy"
("top-cell" class, one Ext-degree lower). The flags below select what to do
about this, and exactly one of them is required:

  --full-ext       Input is old-style Ext(BP_*) data (or you want the raw,
                    unprojected Ext(BP_*/3^40) output). Runs the pre-existing
                    code path unchanged. If second-copy (top-cell) classes are
                    found anyway, prints a prominent warning but still runs
                    (and the output is not filtered -- it will contain the
                    second copies).
  --ext-mod-3-40   Input is Ext(BP_*/3^40) data and you want it projected
                    onto the first copy (so it can be compared with/used like
                    old Ext(BP_*) data). Builds first-copy generators from the
                    Bockstein side and re-expresses everything in a fresh
                    basis for just that subspace.

With neither flag the script exits immediately, before writing anything.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from pathlib import Path

OUTPUTS = [
    {"method": "alpha1", "shift": (3, 1)},
    {"method": "beta1", "shift": (10, 2)},
    {"method": "beta2", "shift": (26, 2)},
    {"method": "beta33", "shift": (34, 2)},
    {"method": "beta4", "shift": (58, 2)},
    {"method": "beta5", "shift": (74, 2)},
    {"method": "beta63", "shift": (82, 2)},
]

BETA_THETA_FILE = {
    "beta1": "_BPBocSS_theta2.txt",
    "beta2": "_BPBocSS_theta3.txt",
    "beta33": "_BPBocSS_theta4.txt",
    "beta4": "_BPBocSS_theta5.txt",
    "beta5": "_BPBocSS_theta6.txt",
    "beta63": "_BPBocSS_theta7.txt",
}

# Strip a leading v0^k power off a BocSS name, e.g. "v0^37v1^9[0-0]" ->
# "v1^9[0-0]". Used to detect the top-cell classification rule below.
_V0_PREFIX_RE = re.compile(r"^v0\^\d+")

# Leading numeral of a bracketed name, e.g. "v1^9[0-0]" -> 0 (the ANSS/BocSS
# s-filtration; the second bracket number is a distinct index and not used
# here).
_BRACKET_S_RE = re.compile(r"\[(\d+)-")

# The literal "|deg=N" trailer on a BPBocSS_table.txt line.
_DEG_RE = re.compile(r"\|deg=(-?\d+)\s*$")

_A0_CLOSURE_CAP = 45


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def find_prefix(source_dir):
    for p in source_dir.glob("*_BPAANSS_table.txt"):
        return p.name[: -len("_BPAANSS_table.txt")]
    raise SystemExit(f"no *_BPAANSS_table.txt found in {source_dir}")


def parse_terms(raw):
    """`a+b+o` -> {a:1, b:1} (mod 3, o dropped, zero coefficients dropped)."""
    d = {}
    for t in raw.split("+"):
        if t == "o" or t == "":
            continue
        d[t] = (d.get(t, 0) + 1) % 3
    return {k: v for k, v in d.items() if v}


def parse_mult_file(path):
    d = {}
    with path.open() as fh:
        for line in fh:
            src, target = line.rstrip("\n").split("\t->\t", 1)
            d[src] = parse_terms(target)
    return d


def parse_mult_file_filtered(path, needed_keys):
    d = {}
    with path.open() as fh:
        for line in fh:
            src, target = line.rstrip("\n").split("\t->\t", 1)
            if src in needed_keys:
                d[src] = parse_terms(target)
    return d


def invert(pairs, exclude_fn):
    """Port of ANSSData.invert, operating on plain dict[str,int] linear combos."""
    inv = {}
    for src, v in pairs:
        targets = list(v.keys())
        if not targets:
            continue
        elif len(targets) == 1:
            t0 = targets[0]
            if t0 in inv:
                raise Exception(f"Failed to invert dict; {t0} is hit multiple times")
            # src -> c*t0 inverts to t0 -> c^-1*src, so this should be the
            # INVERSE of v[t0], not v[t0] itself. It is correct as written only
            # because these coefficients are units mod 3, where 1^-1=1 and
            # 2^-1=2. It would be wrong at any other prime. Same issue at
            # `scale` in the multi-target branch below.
            inv[t0] = {src: v[t0]}
        else:
            old, new = [], []
            for val in targets:
                if exclude_fn(val):
                    continue
                if val in inv:
                    old.append(val)
                else:
                    new.append(val)
            if not new:
                continue
            elif len(new) == 1:
                new_val = new[0]
                new_keys = {src: 1}
                for val in old:
                    coeff = (-v[val]) % 3
                    for k2, c2 in inv[val].items():
                        new_keys[k2] = (new_keys.get(k2, 0) + c2 * coeff) % 3
                # solving src = sum_i c_i t_i for the single unknown t_new
                # divides through by t_new's coefficient, so this should be
                # scale^-1; multiplying by scale is equivalent only mod 3 --
                # see the note in the single-target branch above.
                scale = v[new_val] % 3
                new_keys = {k: (c * scale) % 3 for k, c in new_keys.items()}
                new_keys = {k: c for k, c in new_keys.items() if c}
                inv[new_val] = new_keys
            else:
                # beta.py's original behaviour stores a "?" placeholder here
                # (an unresolvable multi-term ambiguity, expected only at the
                # edge of the computation); skip these rather than poison
                # downstream sums with a non-numeric coefficient.
                print(f"  (invert: stuck on {src}, skipping {len(new)} ambiguous targets)",
                      flush=True)
    return inv


def lc_map(lc, f):
    """Apply f (name -> dict[str,int]) linearly to lc, mod 3."""
    rtn = {}
    for name, c in lc.items():
        for k2, c2 in f(name).items():
            rtn[k2] = (rtn.get(k2, 0) + c * c2) % 3
    return {k: v for k, v in rtn.items() if v}


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert <prefix>_BP*.txt ANSS tables into E5-style CSVs.")
    p.add_argument("dir", type=Path, help="source directory")
    p.add_argument("out_dir", type=Path, nargs="?", default=None,
                    help="output directory (default: same as dir)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--full-ext", action="store_true",
                    help="input is Ext(BP_*)-style data (or: keep the raw, "
                         "unprojected Ext(BP_*/3^40) output). Runs the "
                         "original, unprojected code path; warns (but still "
                         "runs) if second-copy classes are found.")
    g.add_argument("--ext-mod-3-40", action="store_true",
                    help="input is Ext(BP_*/3^40) data; project onto the "
                         "first copy (comparable to old Ext(BP_*) data).")
    args = p.parse_args()
    if not (args.full_ext or args.ext_mod_3_40):
        p.exit(2, "convert_e5.py: nothing written -- pass exactly one of:\n"
               "  --ext-mod-3-40   Ext(BP_*/3^40) data (current mr_BP output), "
               "projected onto the first copy\n"
               "  --full-ext       old Ext(BP_*) data, or Ext(BP_*/3^40) data with "
               "the second copies kept\n")
    return args


def main():
    args = parse_args()
    source_dir = args.dir
    out_dir = args.out_dir if args.out_dir is not None else source_dir
    out_dir.mkdir(exist_ok=True)
    prefix = find_prefix(source_dir)
    log(f"prefix: {prefix!r}")

    log("Parsing degree table...")
    deg_table = {}  # name -> (s, f, adams, index)
    basis_by_degree = {}  # (s,f) -> list[name]
    max_adams = 0
    with (source_dir / f"{prefix}_BPAANSS_table.txt").open() as fh:
        for line in fh:
            if "<-" in line:
                continue
            name, rest = line.split("\t", 1)
            deg_start = rest.index("deg=(")
            deg_str = rest[deg_start + 5:rest.index(")", deg_start)]
            s_str, adams_str = deg_str.split(",")
            stem = int(s_str)
            adams = int(adams_str)
            # ANSS filtration = leading numeral in the bracketed name, e.g. v0^3[7-2] -> 7
            bracket = name[name.index("[") + 1:name.index("-", name.index("["))]
            filt = int(bracket)
            basis = basis_by_degree.setdefault((stem, filt), [])
            index = len(basis)
            basis.append(name)
            deg_table[name] = (stem, filt, adams, index)
            if adams > max_adams:
                max_adams = adams
    log(f"deg_table: {len(deg_table)} generators, max Adams filtration {max_adams}")

    def alias(name):
        s, f, a, i = deg_table[name]
        return f"{s}_{f}_{i}"

    log("Writing E5_rank.csv...")
    rank_rows = sorted((s, f, len(names)) for (s, f), names in basis_by_degree.items())
    with (out_dir / "E5_rank.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s", "f", "rank"])
        w.writerows(rank_rows)
    log(f"E5_rank.csv: {len(rank_rows)} rows")

    log("Parsing alpha1 (h0) table...")
    h0table = parse_mult_file(source_dir / f"{prefix}_BPAANSS_h0.txt")
    log(f"h0table: {len(h0table)} entries")

    log("Parsing algNSS a0 (multiply-by-3) table...")
    a0table = parse_mult_file(source_dir / f"{prefix}_BPAANSS_a0.txt")
    a0div_targets = {t for terms in a0table.values() for t in terms}
    log(f"a0table: {len(a0table)} entries, {len(a0div_targets)} a0-divisible targets")

    log("Parsing B2A table...")
    b2a_pairs = []
    with (source_dir / f"{prefix}_BPB2A_table.txt").open() as fh:
        for line in fh:
            src, target = line.rstrip("\n").split("\t->\t", 1)
            b2a_pairs.append((src, parse_terms(target)))
    B2A = dict(b2a_pairs)
    log(f"B2A: {len(B2A)} entries")

    def exclude_for_b2a_inv(x):
        if x not in deg_table:
            return True
        s, f, adams, i = deg_table[x]
        return adams >= max_adams or x in a0div_targets

    log("Inverting B2A -> B2A_inv...")
    B2A_inv = invert(iter(b2a_pairs), exclude_for_b2a_inv)
    log(f"B2A_inv: {len(B2A_inv)} entries")

    needed_keys = set()
    for lc in B2A_inv.values():
        needed_keys.update(lc.keys())
    log(f"needed theta keys (image of B2A_inv): {len(needed_keys)}")

    log("Parsing bockstein table...")
    bockstein = {}
    boc_deg = {}  # BocSS name -> literal deg=N off any line mentioning it

    def record_deg(name, deg):
        if deg is not None:
            boc_deg.setdefault(name, deg)

    with (source_dir / f"{prefix}_BPBocSS_table.txt").open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if "<-" in line:
                parts = line.split("\t")
                target, _, src, dlen = parts[0], parts[1], parts[2], parts[3]
                deg = None
                if len(parts) > 4:
                    m = _DEG_RE.search(parts[4])
                    if m:
                        deg = int(m.group(1))
                # The literal deg=N on a differential line is written for the
                # target, but (given source's s is always target's s - 1) the
                # same literal, combined with EITHER name's own bracket s via
                # the stem formula below, recovers that name's own stem --
                # see NOTES.md for the derivation. So it's safe to record it
                # for both names.
                record_deg(target, deg)
                record_deg(src, deg)
                if "d0" not in dlen:
                    bockstein[src] = target
                elif src not in bockstein:
                    # A d0 "differential" is search-bookkeeping noise, not a
                    # real page-1 Bockstein differential (README: "Bockstein
                    # d0's can be ignored ... treat [as] 0"); src has no
                    # other table entry of its own, so register it as a
                    # permanent cycle (delta(src) == 0) instead of leaving
                    # it unregistered, which raised a spurious KeyError.
                    bockstein[src] = ""
            else:
                parts = line.split("\t", 1)
                elt = parts[0]
                deg = None
                if len(parts) > 1:
                    m = _DEG_RE.search(parts[1])
                    if m:
                        deg = int(m.group(1))
                record_deg(elt, deg)
                bockstein[elt] = ""
    bockstein_values = set(bockstein.values())
    log(f"bockstein: {len(bockstein)} entries")

    # --- Classify BocSS survivors as first-copy / top-cell (second-copy). ---
    # Rule (verified on rl=94 data against an independent rank count in every
    # bidegree): a surviving BocSS class (bockstein[name] == "") is top-cell
    # iff its v0-free base is the SOURCE of some non-d0 differential. The
    # sources of non-d0 differentials are exactly the bockstein keys whose
    # value is non-empty (a d0 line never assigns a non-empty value).
    diff_source_bases = {
        _V0_PREFIX_RE.sub("", src) for src, tgt in bockstein.items() if tgt != ""
    }

    def boc_bidegree(name):
        m = _BRACKET_S_RE.search(name)
        s = int(m.group(1))
        deg = boc_deg[name]
        stem = (deg + s) // 2 - s
        return stem, s

    # The raw BocSS table's universe of names goes far beyond the (finite,
    # degree-truncated) range the AANSS table was actually resolved to (e.g.
    # very high v0-power towers) -- those extra names are not part of this
    # computation's scope at all, so both the classification counts below and
    # the per-bidegree first-copy count used for the rank cross-check are
    # restricted to bidegrees basis_by_degree actually has AANSS generators
    # for. (Verified: on rl=94 data this restriction takes 462598 raw
    # survivors down to exactly 1118 -- 539 top-cell + 579 first-copy --
    # matching an independent rank count in every bidegree.)
    survivors = [name for name, tgt in bockstein.items()
                 if tgt == "" and boc_bidegree(name) in basis_by_degree]
    topcell_survivors = [n for n in survivors if _V0_PREFIX_RE.sub("", n) in diff_source_bases]
    first_copy_survivors = [n for n in survivors if _V0_PREFIX_RE.sub("", n) not in diff_source_bases]
    n_topcell = len(topcell_survivors)
    log(f"BocSS survivors within the resolved AANSS degree range: {len(survivors)} total, "
        f"{n_topcell} top-cell (second-copy), {len(first_copy_survivors)} first-copy")

    first_copy_boc_counts = {}
    for n in first_copy_survivors:
        bideg = boc_bidegree(n)
        first_copy_boc_counts[bideg] = first_copy_boc_counts.get(bideg, 0) + 1

    if args.ext_mod_3_40:
        mode = "projected"
    elif args.full_ext:
        mode = "unprojected"
        if n_topcell:
            log(f"WARNING: {n_topcell} second-copy (top-cell) classes are present in this "
                f"data (this looks like Ext(BP_*/3^40), not Ext(BP_*)); running --full-ext "
                f"anyway, so the output WILL contain these second copies. Pass "
                f"--ext-mod-3-40 instead to project them out.")

    def exclude_for_boc_a0div(x):
        return not (x in bockstein or x in bockstein_values)

    log("Parsing BocSS a0 table and inverting -> boc_a0div (this is the big one)...")

    def boc_a0_pairs():
        with (source_dir / f"{prefix}_BPBocSS_a0.txt").open() as fh:
            for line in fh:
                src, target = line.rstrip("\n").split("\t->\t", 1)
                yield src, parse_terms(target)

    boc_a0div = invert(boc_a0_pairs(), exclude_for_boc_a0div)
    log(f"boc_a0div: {len(boc_a0div)} entries")

    if mode == "projected":
        # First-copy generator families are seeded from ALL B2A keys (not
        # just the ones hit by B2A_inv), so the theta tables must be filtered
        # to include those too.
        needed_keys.update(B2A.keys())
        log(f"[projected] needed theta keys expanded to include all B2A keys: {len(needed_keys)}")

    log("Parsing the six theta tables (filtered to needed keys)...")
    beta_table = {}
    for method, suffix in BETA_THETA_FILE.items():
        beta_table[method] = parse_mult_file_filtered(source_dir / f"{prefix}{suffix}", needed_keys)
        log(f"  {method}: {len(beta_table[method])} rows kept")

    # bockstein_value(name): the algNSS-terms value of a single Bockstein-
    # space name, i.e. what `name` literally equals (not "name / 3" -- that
    # is delta()'s job, below). Recurses through boc_a0div as many times as
    # `name` itself needs, multiplying the algNSS result by v0 (a0table --
    # B2A has already renamed to algNSS names by this point, so the
    # 3-multiples go back in through the sphere-side a0 table, not the
    # Bockstein one) at *every* level this recurses through.
    #
    # A single a0.txt relation is itself routinely multi-term (the README's
    # "quirks" section: v0^1[5-16] -> [5-16]+v2^7[5-1]+...), so one round of
    # dividing a whole linear combination by boc_a0div can produce a mix of
    # already-irreducible terms (like [5-16]) alongside still-divisible ones
    # (deeper v0-power terms in the same sum) *simultaneously*. Treating the
    # combination as one all-or-nothing unit -- as a single shared
    # `is_a0_divisible`/`a0_count` did previously -- stops dividing the whole
    # sum as soon as any one term becomes irreducible, then hands every term
    # to B2A regardless of its own depth: B2A is only defined for
    # non-3-divisible names, so every remaining divisible term raised a
    # KeyError. Recursing per term reduces each one exactly as many times as
    # it alone needs.
    #
    # A v0-power tower's root is itself a multi-term a0 relation (the same
    # README quirk), and invert() resolves each rung against whichever OTHER
    # rungs already had a one-step answer at the time it ran, without itself
    # recursing -- so this recursion can revisit a name already on the
    # current call stack (a "cycle" in the one-step boc_a0div graph), even
    # though the true, fully-unrolled value is perfectly well-defined and
    # finite. Z3 represents 3-adic coefficients as Z/3^40 (Z3_Op::
    # MAX_VALUATION=39 in the C++; see NOTES.md) -- v0 is nilpotent of order
    # 40, and each round of this recursion corresponds to one more factor of
    # v0, so a name that keeps recurring must genuinely vanish once revisited
    # more times than Z3 can represent; that's the same nilpotence that makes
    # nested a0-division terminate at all, not an approximation. Rather than
    # solve the linear system this implies over F_3 directly, count each
    # name's occurrences on the CURRENT path and truncate (treat as zero)
    # once any name recurs more than Z3 can represent -- comfortably capped
    # here, with margin, at 50 -- and log how often/deep this fires so a
    # real infinite regression (a bug, not nilpotence) would stand out.
    #
    # A name's fully-converged value doesn't actually depend on which
    # caller asks for it or how deep in some unrelated chain that caller
    # happens to be -- only on how many more genuine divisions *that name
    # itself* can still support before nilpotence. In particular, whether a
    # given name hits the revisit cap is itself a fixed property of the
    # one-step boc_a0div graph (every top-level call starts _occurrences
    # fresh and walks the same deterministic dependency structure), so a
    # name that truncates does so identically on every call -- recomputing
    # it never yields a more-exact answer, only the same truncated one at
    # full cost again. So every result is cached regardless of `exact`
    # (previously only exact results were cached, which meant any name
    # anywhere in a large shared subgraph that ever truncated -- expected
    # to happen periodically, per the nilpotence argument above -- poisoned
    # the cache for that whole subtree and forced it to be re-walked from
    # scratch on every subsequent top-level delta() call; on the full 185
    # dataset this turned an O(V) traversal into something that hadn't
    # finished after 7+ hours). `exact` is still tracked and returned so
    # callers/logging can tell a cached truncated value from a fully
    # resolved one.
    _REVISIT_CAP = 50
    _truncation_hits = []
    _value_cache = {}  # name -> (result, exact)

    def bockstein_value(name, _occurrences=None):
        if name in _value_cache:
            return _value_cache[name]
        if name not in boc_a0div:
            result = dict(B2A[name])
            _value_cache[name] = (result, True)
            return result, True
        if _occurrences is None:
            _occurrences = {}
        seen = _occurrences.get(name, 0)
        if seen >= _REVISIT_CAP:
            _truncation_hits.append(name)
            return {}, False
        _occurrences[name] = seen + 1
        result = {}
        exact = True
        for t, c in boc_a0div[name].items():
            sub, sub_exact = bockstein_value(t, _occurrences)
            exact = exact and sub_exact
            remult = lc_map(sub, lambda s: a0table.get(s, {}))
            for k2, c2 in remult.items():
                result[k2] = (result.get(k2, 0) + c * c2) % 3
        _occurrences[name] = seen
        result = {k: v for k, v in result.items() if v}
        _value_cache[name] = (result, exact)
        return result, exact

    def delta(elt):
        # delta(elt) == bockstein[elt] / 3 (README: "d_r(x)=y ... hence
        # delta(x) = y/3"), i.e. exactly ONE division must be "spent" and
        # NOT multiplied back in -- unlike every deeper division inside
        # bockstein_value(), which reconstructs a term's own true value and
        # so *does* put every division's v0 back. bockstein[elt] (a
        # Bockstein differential's target) is a0-divisible by construction;
        # if it weren't, delta() would be asking for y/3 with no recorded
        # way to divide y by 3 at all -- a genuine gap, not a bug in this
        # function, so it's reported the same way as any other missing
        # lookup rather than silently guessed at.
        if elt not in bockstein:
            raise KeyError(f"{elt} not found as a Bockstein element")
        target = bockstein[elt]
        if target == "":
            return {}
        if target not in boc_a0div:
            raise KeyError(f"{target} (bockstein target of {elt!r}) is not a0-divisible")
        result = {}
        for t, c in boc_a0div[target].items():
            v, _exact = bockstein_value(t)
            for k2, c2 in v.items():
                result[k2] = (result.get(k2, 0) + c * c2) % 3
        return {k: v for k, v in result.items() if v}

    def beta_n(theta_table, elt):
        rtn = B2A_inv.get(elt, {})
        rtn = lc_map(rtn, lambda s: theta_table.get(s, {}))
        rtn = lc_map(rtn, lambda s: delta(s))
        return rtn

    def alpha1(elt):
        return h0table.get(elt, {})

    all_failures = set()

    if mode == "unprojected":
        log("Computing product tables...")
        for config in OUTPUTS:
            method_name = config["method"]
            shift_s, shift_f = config["shift"]
            multiplier = basis_by_degree[(shift_s, shift_f)][0]
            multiplier_alias = alias(multiplier)
            multiplier_adams = deg_table[multiplier][2]

            if method_name == "alpha1":
                method = alpha1
            else:
                # Normalise theta_n so that beta_n * 1 is exactly the named
                # multiplier class (the downstream validator requires the
                # identity row to read 1*multiplier). Without this, the beta
                # tables come out as a uniform unit multiple of the intended one.
                unit = basis_by_degree[(0, 0)][0]
                unit_terms = {k: v % 3 for k, v in beta_n(beta_table[method_name], unit).items() if v % 3}
                c = unit_terms.get(multiplier, 0)
                if c == 0:
                    raise Exception(f"{method_name} * 1 has no {multiplier_alias} term: {unit_terms}")
                if set(unit_terms) != {multiplier}:
                    log(f"WARNING: {method_name} * 1 = {unit_terms}, not a multiple of "
                        f"{multiplier_alias} alone; normalising by its coefficient only")
                inv_c = pow(c, -1, 3)
                log(f"{method_name}: beta_n * 1 has coefficient {c}; scaling table by {inv_c}")

                def method(x, m=method_name, s=inv_c):
                    return {k: (v * s) % 3 for k, v in beta_n(beta_table[m], x).items()}

            product_rows = []
            zero_count = 0
            error_count = 0

            for source, (s1, f1, adams1, idx1) in deg_table.items():
                s2, f2 = s1 + shift_s, f1 + shift_f
                try:
                    terms = method(source)
                except KeyError:
                    minimum = adams1 + multiplier_adams
                    receivers = [
                        n for n in basis_by_degree.get((s2, f2), ())
                        if deg_table[n][2] >= minimum
                    ]
                    if receivers:
                        error_count += 1
                        all_failures.add((multiplier_alias, alias(source)))
                    else:
                        zero_count += 1
                    continue

                terms = {k: v for k, v in terms.items() if v}
                if not terms:
                    zero_count += 1
                    continue

                basis = basis_by_degree.get((s2, f2), [])
                expr_terms = [f"{terms[t]}*{alias(t)}" for t in basis if terms.get(t)]
                expr = "+".join(expr_terms) if expr_terms else "0"
                product_rows.append((s1, f1, alias(source), expr))

            product_rows.sort(key=lambda r: (r[0], r[1]))
            out = out_dir / f"E5_{method_name}_products.csv"
            with out.open("w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["factor1", "factor2", "result"])
                for _s1, _f1, factor2, expr in product_rows:
                    w.writerow([multiplier_alias, factor2, expr])
            log(f"{method_name}: {len(product_rows)} nonzero rows, {zero_count} zero omitted, "
                f"{error_count} undetermined")

    else:  # mode == "projected"
        log("[projected] Building first-copy generators from B2A keys...")
        METHODS = [c["method"] for c in OUTPUTS]

        def h0_extend(vec):
            return lc_map(vec, lambda n: h0table.get(n, {}))

        def bidegree_of_vec(vec, context):
            bideg = None
            for name in vec:
                if name not in deg_table:
                    raise SystemExit(f"[projected] {context}: name {name!r} not in AANSS deg_table")
                b = (deg_table[name][0], deg_table[name][1])
                if bideg is None:
                    bideg = b
                elif b != bideg:
                    log(f"WARNING: [projected] {context}: vec spans multiple bidegrees "
                        f"({bideg} vs {b}, term {name!r}) -- single-bidegree assumption violated")
            return bideg

        generator_rows = []  # list of {"vec": {...}, "payload": {method: dict|None}}
        cap_hit_seeds = 0
        for x, vec0 in B2A.items():
            bidegree_of_vec(vec0, f"B2A[{x!r}]")
            payload0 = {"alpha1": h0_extend(vec0)}
            for m in METHODS:
                if m == "alpha1":
                    continue
                seed_lc = beta_table[m].get(x, {})
                try:
                    payload0[m] = lc_map(seed_lc, delta)
                except KeyError:
                    payload0[m] = None

            vec, payload = vec0, payload0
            i = 0
            hit_cap = False
            while True:
                generator_rows.append({"vec": vec, "payload": payload})
                vec_next = lc_map(vec, lambda n: a0table.get(n, {}))
                if not vec_next:
                    break
                i += 1
                if i > _A0_CLOSURE_CAP:
                    hit_cap = True
                    break
                payload_next = {"alpha1": h0_extend(vec_next)}
                for m in METHODS:
                    if m == "alpha1":
                        continue
                    prev = payload[m]
                    payload_next[m] = None if prev is None else lc_map(prev, lambda n: a0table.get(n, {}))
                vec, payload = vec_next, payload_next
            if hit_cap:
                cap_hit_seeds += 1
                log(f"WARNING: [projected] seed {x!r} hit the {_A0_CLOSURE_CAP}-iteration "
                    f"a0-closure cap without emptying; stopped")
        log(f"[projected] {len(generator_rows)} generator rows from {len(B2A)} B2A seeds"
            + (f", {cap_hit_seeds} seed(s) hit the iteration cap" if cap_hit_seeds else ""))

        rows_by_bidegree = {}
        for row in generator_rows:
            bideg = bidegree_of_vec(row["vec"], "generator row")
            if bideg is None:
                continue
            rows_by_bidegree.setdefault(bideg, []).append(row)

        def combine(row_a, row_b, c):
            """row_a + c * row_b, mod 3; None payload propagates."""
            if c % 3 == 0:
                return row_a
            vec = dict(row_a["vec"])
            for k, v in row_b["vec"].items():
                vec[k] = (vec.get(k, 0) + c * v) % 3
            vec = {k: v for k, v in vec.items() if v}
            payload = {}
            for m in METHODS:
                pa, pb = row_a["payload"][m], row_b["payload"][m]
                if pa is None or pb is None:
                    payload[m] = None
                else:
                    d = dict(pa)
                    for k, v in pb.items():
                        d[k] = (d.get(k, 0) + c * v) % 3
                    payload[m] = {k: v for k, v in d.items() if v}
            return {"vec": vec, "payload": payload}

        def scale(row, c):
            vec = {k: (v * c) % 3 for k, v in row["vec"].items()}
            vec = {k: v for k, v in vec.items() if v}
            payload = {}
            for m in METHODS:
                p = row["payload"][m]
                if p is None:
                    payload[m] = None
                else:
                    payload[m] = {k: (v * c) % 3 for k, v in p.items() if (v * c) % 3}
            return {"vec": vec, "payload": payload}

        relation_violation_count = 0
        relation_violation_examples = []
        new_basis_order = {}  # (stem,filt) -> [pivot AANSS name, ...] in new-basis order
        new_basis_rows = {}   # (stem,filt) -> [row, ...] parallel to new_basis_order

        log("[projected] Gaussian-eliminating generator rows per bidegree...")
        for bideg, rows in rows_by_bidegree.items():
            stem, filt = bideg
            names = basis_by_degree.get((stem, filt), [])
            col_rank = {n: (deg_table[n][2], deg_table[n][3]) for n in names}
            pivots = {}  # AANSS name -> reduced row
            for row in rows:
                work = row
                # Clear EVERY already-pivoted column, not just the row's
                # lowest-ranked one: a row whose lowest column is new can
                # still carry a component on some other pivot column, and
                # leaving it there breaks the "coefficient = value at the
                # pivot column" readout used to express products below.
                while True:
                    nz = [c for c, v in work["vec"].items() if v]
                    piv_cols = [c for c in nz if c in pivots]
                    if not piv_cols:
                        break
                    leading = min(piv_cols, key=lambda c: col_rank.get(c, (10**9, 10**9)))
                    factor = work["vec"][leading]  # pivot's own coeff there is 1
                    work = combine(work, pivots[leading], (-factor) % 3)
                nz = [c for c, v in work["vec"].items() if v]
                if not nz:
                    # A relation among generators: vec collapses to 0, so any
                    # KNOWN nonzero payload here means the products computed
                    # from this generator family aren't actually consistent.
                    for m in METHODS:
                        p = work["payload"][m]
                        if p is not None and p:
                            relation_violation_count += 1
                            if len(relation_violation_examples) < 8:
                                relation_violation_examples.append((bideg, m, dict(p)))
                    continue
                leading = min(nz, key=lambda c: col_rank.get(c, (10**9, 10**9)))
                inv_lead = pow(work["vec"][leading], -1, 3)
                work = scale(work, inv_lead)
                for col, prow in list(pivots.items()):
                    coeff = prow["vec"].get(leading, 0)
                    if coeff:
                        pivots[col] = combine(prow, work, (-coeff) % 3)
                pivots[leading] = work
            order = sorted(pivots.keys(), key=lambda n: deg_table[n][3])
            new_basis_order[bideg] = order
            new_basis_rows[bideg] = [pivots[n] for n in order]

        if relation_violation_count:
            log(f"WARNING: [projected] {relation_violation_count} generator relation(s) not "
                f"respected by products (examples: {relation_violation_examples})")
        else:
            log("[projected] no generator-relation violations (all known payloads vanish "
                "where the generator vector does)")

        total_rank = sum(len(v) for v in new_basis_order.values())
        log(f"[projected] total first-copy rank: {total_rank}")

        mismatch_count = 0
        mismatch_examples = []
        for bideg in sorted(set(new_basis_order) | set(first_copy_boc_counts)):
            got = len(new_basis_order.get(bideg, []))
            want = first_copy_boc_counts.get(bideg, 0)
            if got != want:
                mismatch_count += 1
                if len(mismatch_examples) < 10:
                    mismatch_examples.append((bideg, got, want))
        if mismatch_count:
            log(f"WARNING: [projected] dimension check: {mismatch_count} bidegree(s) where "
                f"reconstructed rank != first-copy BocSS survivor count "
                f"(bideg, got, want): {mismatch_examples}")
        else:
            log("[projected] dimension check: reconstructed rank matches first-copy BocSS "
                "survivor counts in every bidegree")

        log("[projected] Writing (projected) E5_rank.csv...")
        rank_rows = sorted((s, f, len(order)) for (s, f), order in new_basis_order.items())
        with (out_dir / "E5_rank.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["s", "f", "rank"])
            w.writerows(rank_rows)
        log(f"[projected] E5_rank.csv: {len(rank_rows)} rows, total rank {total_rank}")

        log("[projected] Computing product tables...")
        residual_violation_count = 0
        residual_examples = []
        for config in OUTPUTS:
            method_name = config["method"]
            shift_s, shift_f = config["shift"]
            target_home = (shift_s, shift_f)
            mult_order = new_basis_order.get(target_home, [])
            if not mult_order:
                raise SystemExit(f"[projected] no first-copy basis at shift bidegree "
                                  f"{target_home} for {method_name}")
            multiplier_alias = f"{shift_s}_{shift_f}_0"
            multiplier_pivot = mult_order[0]
            multiplier_adams = deg_table[multiplier_pivot][2]

            if method_name == "alpha1":
                inv_c = 1
            else:
                unit_order = new_basis_order.get((0, 0), [])
                if "[0-0]" not in unit_order:
                    raise SystemExit("[projected] no [0-0] pivot in bidegree (0,0); "
                                      "cannot normalise beta tables")
                unit_idx = unit_order.index("[0-0]")
                unit_row = new_basis_rows[(0, 0)][unit_idx]
                unit_val = unit_row["payload"][method_name]
                if unit_val is None:
                    raise SystemExit(f"[projected] {method_name} * 1 is undetermined (None); "
                                      f"cannot normalise")
                coeffs = {j: unit_val.get(pn, 0) for j, pn in enumerate(mult_order)}
                nz_coeffs = {j: c for j, c in coeffs.items() if c}
                if not nz_coeffs:
                    raise SystemExit(f"[projected] {method_name} * 1 has no component on "
                                      f"{target_home}'s new basis: {unit_val}")
                residual = dict(unit_val)
                for j, c in nz_coeffs.items():
                    for k, v in new_basis_rows[target_home][j]["vec"].items():
                        residual[k] = (residual.get(k, 0) - c * v) % 3
                residual = {k: v for k, v in residual.items() if v}
                if residual:
                    log(f"WARNING: [projected] {method_name}*1 normalisation residual "
                        f"nonzero: {residual}")
                c0 = nz_coeffs.get(0, 0)
                if c0 == 0:
                    raise Exception(f"{method_name} * 1 has no {multiplier_alias} term "
                                     f"(new-basis coeffs: {nz_coeffs})")
                if set(nz_coeffs) != {0}:
                    log(f"WARNING: {method_name} * 1 = {nz_coeffs} (new-basis indices), not "
                        f"a multiple of {multiplier_alias} alone; normalising by its "
                        f"coefficient only")
                inv_c = pow(c0, -1, 3)
                log(f"{method_name}: beta_n * 1 has coefficient {c0}; scaling table by {inv_c}")

            product_rows = []
            zero_count = 0
            error_count = 0
            undetermined_count = 0

            for (s1, f1), order in new_basis_order.items():
                rows = new_basis_rows[(s1, f1)]
                s2, f2 = s1 + shift_s, f1 + shift_f
                target_order = new_basis_order.get((s2, f2), [])
                target_rows = new_basis_rows.get((s2, f2), [])
                for idx1, pivot_name in enumerate(order):
                    row = rows[idx1]
                    val = row["payload"][method_name]
                    source_alias = f"{s1}_{f1}_{idx1}"
                    if val is None:
                        minimum = deg_table[pivot_name][2] + multiplier_adams
                        receivers = [n for n in target_order if deg_table[n][2] >= minimum]
                        if receivers:
                            error_count += 1
                            undetermined_count += 1
                            all_failures.add((multiplier_alias, source_alias))
                        else:
                            zero_count += 1
                        continue

                    val = {k: (v * inv_c) % 3 for k, v in val.items() if (v * inv_c) % 3}
                    if not val:
                        zero_count += 1
                        continue

                    coeffs = [(val.get(pn, 0)) for pn in target_order]
                    nz = [(c, j) for j, c in enumerate(coeffs) if c]
                    if not nz:
                        zero_count += 1
                        continue

                    residual = dict(val)
                    for c, j in nz:
                        for k, v in target_rows[j]["vec"].items():
                            residual[k] = (residual.get(k, 0) - c * v) % 3
                    residual = {k: v for k, v in residual.items() if v}
                    if residual:
                        residual_violation_count += 1
                        if len(residual_examples) < 8:
                            residual_examples.append((method_name, source_alias, dict(residual)))

                    expr_terms = [f"{c}*{s2}_{f2}_{j}" for c, j in sorted(nz, key=lambda t: t[1])]
                    expr = "+".join(expr_terms)
                    product_rows.append((s1, f1, source_alias, expr))

            product_rows.sort(key=lambda r: (r[0], r[1]))
            out = out_dir / f"E5_{method_name}_products.csv"
            with out.open("w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["factor1", "factor2", "result"])
                for _s1, _f1, factor2, expr in product_rows:
                    w.writerow([multiplier_alias, factor2, expr])
            log(f"{method_name}: {len(product_rows)} nonzero rows, {zero_count} zero omitted, "
                f"{undetermined_count} undetermined")

        if residual_violation_count:
            log(f"WARNING: [projected] {residual_violation_count} product(s) leave the first "
                f"copy (nonzero residual after projecting to the new basis) -- examples: "
                f"{residual_examples}")
        else:
            log("[projected] no residual violations: every product landed exactly in the "
                "first-copy span")

    fail_out = out_dir / "E5_product_failures.csv"
    with fail_out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["factor1", "factor2"])
        for factor1, factor2 in sorted(all_failures):
            w.writerow([factor1, factor2])
    log(f"E5_product_failures.csv: {len(all_failures)} rows")
    if _truncation_hits:
        log(f"bockstein_value: nilpotence truncation fired {len(_truncation_hits)} times "
            f"({len(set(_truncation_hits))} distinct names) -- see NOTES.md")
    log("Done.")


if __name__ == "__main__":
    # bockstein_value() recurses around each v0-tower cycle in boc_a0div up
    # to _REVISIT_CAP times, so its depth is ~_REVISIT_CAP x (tower length);
    # at 270 towers reach ~150 names, far past Python's default limit of
    # 1000. Run in a thread with a large C stack so deep recursion can't
    # overflow it.
    import threading
    sys.setrecursionlimit(1_000_000)
    threading.stack_size(1 << 30)
    # An exit or exception inside a thread doesn't set the process exit
    # status, so carry it back to the main thread explicitly.
    _status = [0]

    def _run():
        try:
            main()
        except SystemExit as e:
            if e.code is not None and not isinstance(e.code, int):
                print(e.code, file=sys.stderr)
                _status[0] = 1
            else:
                _status[0] = e.code or 0
        except BaseException:
            import traceback
            traceback.print_exc()
            _status[0] = 1

    _t = threading.Thread(target=_run)
    _t.start()
    _t.join()
    sys.exit(_status[0])
