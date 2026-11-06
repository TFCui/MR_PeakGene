#!/usr/bin/env python3
"""Manuscript checks (standard library only).

Run from Manuscript/ (or anywhere; paths are resolved from this file):

  python tools/ms_check.py invariants [--base REF | --snapshot DIR] [--summary]
      Compare sections/*.tex and tables/*.tex in the working copy with the
      version at git REF (default HEAD), or, with --snapshot, with the copies
      in DIR (paths mirrored, e.g. DIR/sections/results.tex; only the files
      present in DIR are compared). Numbers, \\cite keys, \\ref/\\label/
      \\nameref targets and math contents are compared as multisets.
      Prints added / removed items per file; exit 1 if anything differs.
  python tools/ms_check.py placeholders
      List \\TODO, XXX, (citation), (reference), (version), (Later),
      [CITE: ...] and author self-notes "(we need ...)", "(add ...)",
      "(later ...)" as file:line.
  python tools/ms_check.py consistency
      Stub (not implemented yet).

What counts as a number: integers, decimals, thousands separators
(40,087), scientific notation (1e-8), percentages (73\\%), genomic sizes (1MB, 500kb, 20 Mb), a leading
+/-/\\pm/±, plus SNP and region identifiers (rs705704,
chr11:117185715-117186216). Digits inside words (CD4, RPL11, hg38), lengths
(1.5pt, 0.7\\textheight), LaTeX comments, and arguments of \\cite / \\ref /
\\label / \\includegraphics / \\input are not counted as numbers. Numbers
inside math are covered by the math comparison, not the number list.
"""

import argparse
import collections
import re
import subprocess
import sys
from pathlib import Path

MS = Path(__file__).resolve().parent.parent  # Manuscript/
GLOBS = ("sections/*.tex", "tables/*.tex")

# ---------------------------------------------------------------- parsing

COMMENT = re.compile(r"(?<!\\)%.*")

MATH_ENVS = ("equation", "align", "gather", "multline", "eqnarray",
             "displaymath", "math", "alignat", "flalign")
MATH_PATTERNS = [
    re.compile(r"\\begin\{(%s)(\*?)\}(.*?)\\end\{\1\2\}" % "|".join(MATH_ENVS), re.S),
    re.compile(r"\\\[(.*?)\\\]", re.S),
    re.compile(r"\\\((.*?)\\\)", re.S),
    re.compile(r"(?<!\\)\$\$(.*?)(?<!\\)\$\$", re.S),
    re.compile(r"(?<!\\)\$((?:[^$\\]|\\.)+?)(?<!\\)\$", re.S),
]

CITE = re.compile(r"\\[cC]ite[a-zA-Z]*\*?(?:\s*\[[^\]]*\]){0,2}\s*\{([^}]*)\}")
REF = re.compile(r"\\(ref|eqref|pageref|autoref|cref|Cref|nameref|label)\*?\s*\{([^}]*)\}")
# command arguments that are file names, keys or layout, not numbers
STRIP_ARGS = re.compile(
    r"\\(?:[cC]ite[a-zA-Z]*|ref|eqref|pageref|autoref|cref|Cref|nameref|label|"
    r"includegraphics|input|include|graphicspath|usepackage|documentclass|"
    r"begin|end|setlength|addtolength|hspace|vspace|rule)\*?"
    r"(?:\s*\[[^\]]*\])*(?:\s*\{[^{}]*\})+")
# lengths: "1.5pt" (unit attached, not followed by a letter, so "428 individuals"
# is not "428 in"), or "0.7\textheight"
DIMEN = re.compile(r"[-+]?\d*\.?\d+(?:(?:pt|in|cm|mm|em|ex|sp|pc)(?![A-Za-z])|\s*\\(?:"
                   r"linewidth|textwidth|textheight|columnwidth|paperwidth|baselineskip)\b)")

IDENT = re.compile(r"(?<![A-Za-z0-9_])(?:rs\d+|chr[0-9XYM]+:\d+(?:-\d+)?)(?![A-Za-z0-9_])")
NUMBER = re.compile(
    r"(?<![A-Za-z0-9_.])"
    r"((?:\\pm\s*|±\s*|(?<![0-9])[-+−]\s?)?"
    r"(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?"
    r"(?<=\d)(?:[eE][-+−]?\d+)?"
    r"(?:\s*\\?%|\s?(?:kb|Kb|Mb|MB|Gb|bp)(?![A-Za-z]))?)"
    r"(?![A-Za-z0-9_])")


def strip_comments(text):
    return "\n".join(COMMENT.sub("", line) for line in text.split("\n"))


def norm_space(s):
    return re.sub(r"\s+", " ", s).strip()


def norm_number(s):
    s = s.replace("−", "-").replace("\\pm", "±").replace("\\%", "%")
    return re.sub(r"\s+", "", s)


def extract(text):
    """Return {category: Counter} for one .tex source."""
    text = strip_comments(text)
    out = {k: collections.Counter() for k in ("numbers", "cite", "refs", "math")}

    for m in CITE.finditer(text):
        for key in m.group(1).split(","):
            if key.strip():
                out["cite"][key.strip()] += 1
    for m in REF.finditer(text):
        out["refs"]["%s:%s" % (m.group(1), m.group(2).strip())] += 1

    # math first, so numbers in math are not double counted
    prose = text
    for pat in MATH_PATTERNS:
        def take(m):
            body = m.group(m.lastindex)
            if pat is MATH_PATTERNS[0]:
                body = m.group(1) + m.group(2) + ": " + m.group(3)
            out["math"][norm_space(body)] += 1
            return " "
        prose = pat.sub(take, prose)

    prose = STRIP_ARGS.sub(" ", prose)
    prose = DIMEN.sub(" ", prose)
    for m in IDENT.finditer(prose):
        out["numbers"][m.group(0)] += 1
    prose = IDENT.sub(" ", prose)
    for m in NUMBER.finditer(prose):
        tok = norm_number(m.group(1))
        if re.search(r"\d", tok):
            out["numbers"][tok] += 1
    return out


# ---------------------------------------------------------------- git

def git(*args):
    # Python 3.6 on the cluster: no capture_output / text
    return subprocess.run(["git", "-C", str(MS)] + list(args),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)


def base_files(ref):
    r = git("ls-tree", "-r", "--name-only", ref, "--", "sections", "tables")
    if r.returncode != 0:
        sys.exit("git ls-tree failed for ref %r: %s" % (ref, r.stderr.strip()))
    return {p for p in r.stdout.split() if p.endswith(".tex")}


def base_text(ref, rel):
    r = git("show", "%s:./%s" % (ref, rel))
    return r.stdout if r.returncode == 0 else ""


def work_files():
    return {str(p.relative_to(MS)) for g in GLOBS for p in MS.glob(g)}


# ---------------------------------------------------------------- commands

def snapshot_files(snap):
    return {str(p.relative_to(snap)) for p in snap.rglob("*.tex")}


def cmd_invariants(args):
    if args.snapshot:
        snap = Path(args.snapshot).resolve()
        if not snap.is_dir():
            sys.exit("snapshot directory not found: %s" % snap)
        files = sorted(snapshot_files(snap))
        if not files:
            sys.exit("no .tex files in snapshot %s" % snap)
        label = str(snap)
        old_text = lambda rel: (snap / rel).read_text(encoding="utf-8")
    else:
        files = sorted(base_files(args.base) | work_files())
        label = args.base
        old_text = lambda rel: base_text(args.base, rel)
    changed = 0
    for rel in files:
        path = MS / rel
        new = extract(path.read_text(encoding="utf-8") if path.exists() else "")
        old = extract(old_text(rel))
        if args.summary:
            print("%-42s " % rel + "  ".join(
                "%s=%d" % (k, sum(new[k].values())) for k in new))
        lines = []
        for cat in ("numbers", "cite", "refs", "math"):
            added = new[cat] - old[cat]
            removed = old[cat] - new[cat]
            for sign, items in (("+", added), ("-", removed)):
                for item, n in sorted(items.items()):
                    lines.append("  %s %-7s %s%s" % (sign, cat, item,
                                                     "  (x%d)" % n if n > 1 else ""))
        if lines:
            changed += 1
            print("== %s" % rel)
            print("\n".join(lines))
    if changed:
        print("invariants: %d file(s) differ from %s" % (changed, label))
        return 1
    print("invariants: no changes vs %s (%d files)" % (label, len(files)))
    return 0


PLACEHOLDERS = [
    ("TODO", re.compile(r"\\TODO\b")),
    ("XXX", re.compile(r"XXX")),
    ("(citation)", re.compile(r"\(\s*citation\b[^)]*\)", re.I)),
    ("(reference)", re.compile(r"\(\s*(?:[A-Za-z0-9]+\s+)?reference\b[^)]*\)", re.I)),
    # "(version)", "(hg38 East Asian, version)", "(Version later)", "[version]";
    # not a real version such as "(version 1.0.1)"
    ("(version)", re.compile(r"\((?:[^()]*,\s*)?version(?:\s+later)?\s*\)|\[version\]", re.I)),
    ("(Later)", re.compile(r"\(\s*later\s*\)", re.I)),
    ("[CITE:]", re.compile(r"\[CITE:[^\]]*\]")),
    ("self-note", re.compile(r"\(\s*(?:we need|add|later)\b[^)]*\)", re.I)),
]


def cmd_placeholders(args):
    counts = collections.Counter()
    hits = []
    for rel in sorted(work_files()):
        lines = (MS / rel).read_text(encoding="utf-8").split("\n")
        for i, line in enumerate(lines, 1):
            line = COMMENT.sub("", line)
            spans = []  # count each piece of text once, e.g. "(citation, version)"
            for name, pat in PLACEHOLDERS:
                for m in pat.finditer(line):
                    if any(a <= m.start() and m.end() <= b for a, b in spans):
                        continue
                    spans.append(m.span())
                    counts[name] += 1
                    snippet = norm_space(m.group(0))
                    hits.append("%s:%d: [%s] %s" % (rel, i, name,
                                                    snippet[:90] + ("…" if len(snippet) > 90 else "")))
    print("\n".join(hits))
    print("== placeholders: " + ", ".join(
        "%s %d" % (name, counts[name]) for name, _ in PLACEHOLDERS) +
        "; total %d" % sum(counts.values()))
    return 0


def cmd_consistency(args):
    print("consistency: not implemented yet")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd")  # required= needs Python 3.7
    p = sub.add_parser("invariants", help="numbers/cites/refs/math vs a git ref")
    p.add_argument("--base", default="HEAD", help="git ref to compare with (default HEAD)")
    p.add_argument("--snapshot", metavar="DIR",
                   help="compare with .tex copies in DIR instead of a git ref")
    p.add_argument("--summary", action="store_true", help="also print item counts per file")
    p.set_defaults(func=cmd_invariants)
    sub.add_parser("placeholders", help="list TODOs and placeholders").set_defaults(func=cmd_placeholders)
    sub.add_parser("consistency", help="(stub)").set_defaults(func=cmd_consistency)
    args = ap.parse_args()
    if not args.cmd:
        ap.error("choose a subcommand: invariants, placeholders, consistency")
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
