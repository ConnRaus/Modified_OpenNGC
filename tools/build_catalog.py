#!/usr/bin/env python3
"""
Builds dist/catalog.sqlite (schema v2) from database_files/: NGC.csv, addendum.csv and links.csv.

Schema:
  objects       one row per real object. id is permanent (it is the image key and the saved-state key).
  designations  every catalog id of an object. (catalog, code) is unique, so an id names exactly one object.
                Display text = catalogs.prefix + code, e.g. "NGC " + "2237". rank orders them within an object.
  names         common names of an object; rank 0 is the one to show.
  catalogs      the catalogs designations come from: display prefix, full name, priority.
  types         object type codes and their labels.
  relations     links between objects that stay separate, e.g. ("IC1396A", "IC1396", "part_of").
  meta          build info and schema_version.

Title of an object = names rank 0, else designations rank 0. Ids of objects folded into another object
(OpenNGC Dup rows, links.csv "same" rows) live on as designations of the object they were folded into.

A report of everything the build dropped or had to decide is written to dist/catalog.report.txt.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import sqlite3
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA_VERSION = "2"

# NGC.csv rows kept: canonical and lettered NGC/IC ids. NED component rows ("IC0080 NED01") are skipped.
RE_NGC_IC = re.compile(r"^(NGC|IC)\d+[A-Z]*$")


@dataclass(frozen=True)
class Catalog:
    id: str
    name: str
    prefix: str
    pattern: str
    code: str = "{0}"


# Display priority order: first = most prominent. Identifiers matching none of these are dropped.
CATALOGS: list[Catalog] = [
    Catalog("M", "Messier", "M", r"M\s*(\d+)"),
    Catalog("NGC", "New General Catalogue", "NGC ", r"NGC\s*(\d+)([A-Z]*)", "{0}{1}"),
    Catalog("IC", "Index Catalogue", "IC ", r"IC\s*(\d+)([A-Z]*)", "{0}{1}"),
    Catalog("C", "Caldwell", "C ", r"C\s*(\d+)"),
    Catalog("Sh2", "Sharpless", "Sh2-", r"(?:SH\s*2\s*-|Sharpless\s+)\s*(\d+)"),
    Catalog("B", "Barnard", "B ", r"B\s*(\d+)"),
    Catalog("RCW", "Rodgers, Campbell & Whiteoak", "RCW ", r"RCW\s*(\d+)"),
    Catalog("Gum", "Gum", "Gum ", r"Gum\s*(\d+)"),
    Catalog("LBN", "Lynds Bright Nebulae", "LBN ", r"LBN\s*(\d+)"),
    Catalog("LDN", "Lynds Dark Nebulae", "LDN ", r"LDN\s*(\d+)"),
    Catalog("SNR", "Green Galactic SNRs", "SNR G", r"SNR\s*G\s*(\d+\.\d)([+-])0*(\d+\.\d)", "{0}{1}{2}"),
    Catalog("Mel", "Melotte", "Mel ", r"Mel\s*(\d+)"),
    Catalog("Cr", "Collinder", "Cr ", r"(?:Cl|Cr|Collinder)\s*(\d+)"),
    Catalog("H", "Harvard", "H ", r"H\s*(\d+)"),
    Catalog("W", "Westerhout", "W ", r"(?:Westerhout\s*|W\s*)(\d+)"),
    Catalog("Arp", "Arp", "Arp ", r"Arp\s*(\d+)"),
    Catalog("VV", "Vorontsov-Velyaminov", "VV ", r"VV\s*(\d+)"),
    Catalog("HCG", "Hickson Compact Group", "HCG ", r"HCG\s*(\d+)"),
    Catalog("UGC", "Uppsala General Catalogue", "UGC ", r"UGC\s*(\d+)"),
    Catalog("UGCA", "UGC Appendix", "UGCA ", r"UGCA\s*(\d+)"),
    Catalog("ESO", "ESO/Uppsala Survey", "ESO ", r"ESO\s*(\d+)-(\d+)", "{0}-{1:0>3}"),
    Catalog("PGC", "Principal Galaxies Catalogue", "PGC ", r"PGC\s*(\d+)"),
    Catalog("MWSC", "Milky Way Star Clusters", "MWSC ", r"MWSC\s*(\d+)"),
    Catalog("HIP", "Hipparcos", "HIP ", r"HIP\s*(\d+)"),
]
# Every catalog also takes one optional trailing letter for lettered parts ("Gum 34a", "RCW 38A")
_COMPILED = [(c, re.compile(rf"^{c.pattern}([A-Z]?)$", re.IGNORECASE)) for c in CATALOGS]
PRIORITY = {c.id: i for i, c in enumerate(CATALOGS)}

TYPES: dict[str, str] = {
    "G": "Galaxy",
    "GPair": "Galaxy pair",
    "GTrpl": "Galaxy triplet",
    "GGroup": "Galaxy group",
    "GCl": "Globular cluster",
    "OCl": "Open cluster",
    "*Ass": "Stellar association",
    "Cl+N": "Cluster with nebula",
    "Neb": "Nebula",
    "EmN": "Emission nebula",
    "HII": "HII region",
    "RfN": "Reflection nebula",
    "DrkN": "Dark nebula",
    "PN": "Planetary nebula",
    "SNR": "Supernova remnant",
    "Nova": "Nova",
    "*": "Star",
    "**": "Double star",
    "Other": "Other object",
    "NonEx": "Non-existent object",
    "Dup": "Duplicate object",
}
# Dup rows fold into their target; only a Dup with no findable target stays an object.
DUP = "Dup"

LINK_KINDS = {"same", "part_of", "member_of"}

# Claim strength on a designation, lower wins: an object's own id, a folded row's own id, a cross-id.
OWN, FOLDED, CROSS = 0, 1, 2


@dataclass
class Report:
    lines: dict[str, list[str]] = field(default_factory=lambda: defaultdict(list))

    def add(self, section: str, line: str) -> None:
        self.lines[section].append(line)


def parse_designation(raw: str) -> tuple[Catalog, str] | None:
    s = (raw or "").strip()
    for cat, rx in _COMPILED:
        m = rx.match(s)
        if m:
            *groups, letter = m.groups()
            parts = [str(int(g)) if g.isdigit() else g.upper() for g in groups]
            return cat, cat.code.format(*parts) + letter.upper()
    return None


def key_of(raw: str) -> str | None:
    """Search key of an id in its canonical display form, so "Cl399" and "Cr 399" agree."""
    d = parse_designation(raw)
    if d:
        return search_key(d[0].prefix + d[1])
    # NC001 etc: objects with no catalog number. Keyed for links, but not a designation, so the app shows only the name
    return search_key(raw) if re.fullmatch(r"NC\d+", (raw or "").strip()) else None


def search_key(text: str) -> str:
    """Uppercase, no spaces, no leading zeros in digit runs. The app must normalise the same way."""
    s = re.sub(r"\s+", "", text.upper())
    return re.sub(r"\d+", lambda m: str(int(m.group(0))), s)


def parse_coord(coord: str, is_ra: bool) -> float | None:
    s = (coord or "").strip()
    if not s:
        return None
    neg = s.startswith("-")
    parts = s.lstrip("+-").split(":")
    try:
        if len(parts) != 3:
            return None
        v = float(parts[0]) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0
    except ValueError:
        return None
    if is_ra:
        v *= 15.0
    return -v if neg else v


def fnum(s: str | None) -> float | None:
    try:
        return float((s or "").strip())
    except ValueError:
        return None


def read_rows(csv_path: Path, ngc_only: bool) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            name = (r.get("Name") or "").strip()
            if not name or (ngc_only and not RE_NGC_IC.match(name)):
                continue
            r["Name"] = name
            r["Type"] = (r.get("Type") or "").strip()
            out.append(r)
    return out


def read_links(path: Path, report: Report) -> list[tuple[str, str, str]]:
    if not path.exists():
        return []
    out: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            obj, tgt, kind = (r.get("object") or "").strip(), (r.get("target") or "").strip(), (r.get("kind") or "").strip()
            if not obj or not tgt or kind not in LINK_KINDS:
                report.add("links.csv rows ignored (need object;target;kind with a known kind)", str(r))
                continue
            out.append((obj, tgt, kind))
    return out


def cross_ids(r: dict[str, str], ned_base: bool = False) -> list[str]:
    """Raw ids from the M, NGC and IC columns (comma lists allowed). NED parts are skipped, or with
    ned_base mapped to their parent ("3058 NED02" -> NGC3058)."""
    out: list[str] = []
    m = (r.get("M") or "").strip()
    if m:
        out.append(f"M{m}")
    for col in ("NGC", "IC"):
        for v in (r.get(col) or "").split(","):
            v = v.strip()
            if v and "NED" in v:
                v = v.split("NED")[0].strip() if ned_base else ""
            if v:
                out.append(f"{col}{v}")
    return out


@dataclass
class Obj:
    row: dict[str, str]
    ra: float
    dec: float
    folded: list[dict[str, str]] = field(default_factory=list)


def build(ngc_csv: Path, add_csv: Path, links_csv: Path, report: Report):
    rows: dict[str, dict[str, str]] = {}
    for r in read_rows(ngc_csv, ngc_only=True) + read_rows(add_csv, ngc_only=False):
        k = key_of(r["Name"])
        if k is None:
            report.add("Rows whose name matches no catalog (dropped, add the catalog to CATALOGS)", r["Name"])
            continue
        if k in rows:
            report.add("Repeated row names (first kept)", r["Name"])
            continue
        rows[k] = r

    def make_obj(r: dict[str, str]) -> Obj | None:
        ra, dec = parse_coord(r.get("RA", ""), True), parse_coord(r.get("Dec", ""), False)
        if ra is None or dec is None or not (0 <= ra <= 360 and -90 <= dec <= 90):
            report.add("Rows without usable coordinates (dropped, cannot be placed)", f"{r['Name']} ({r['Type']})")
            return None
        if r["Type"] not in TYPES:
            report.add("Unknown type codes (kept, label = code; fix the CSV)", f"{r['Name']}: {r['Type']}")
        return Obj(r, ra, dec)

    # Every non-Dup row with coordinates is an object (non-existent ones included, typed NonEx).
    objs: dict[str, Obj] = {}
    for r in rows.values():
        if r["Type"] != DUP:
            o = make_obj(r)
            if o:
                objs[r["Name"]] = o

    # Lookup from any id key to the object that owns it as its own id or Messier number.
    by_key: dict[str, str] = {key_of(n): n for n in objs}
    for n, o in objs.items():
        m = (o.row.get("M") or "").strip()
        if m:
            by_key.setdefault(key_of(f"M{m}"), n)

    def resolve(name: str, seen: set[str], ned_base: bool) -> str | None:
        r = rows.get(key_of(name))
        if r is None:
            return None
        for raw in cross_ids(r, ned_base):
            k = key_of(raw)
            if k is None:
                continue
            if k in by_key:
                return by_key[k]
            nxt = rows.get(k)
            if nxt is not None and nxt["Type"] == DUP and nxt["Name"] not in seen:
                hit = resolve(nxt["Name"], seen | {nxt["Name"]}, ned_base)
                if hit:
                    return hit
        return None

    for r in rows.values():
        if r["Type"] != DUP:
            continue
        # A Dup of a NED part falls back to the part's parent object.
        tgt = resolve(r["Name"], {r["Name"]}, False) or resolve(r["Name"], {r["Name"]}, True)
        if tgt:
            objs[tgt].folded.append(r)
            if not resolve(r["Name"], {r["Name"]}, False):
                report.add("Dup rows of a NED part, folded into the parent object", f"{r['Name']} -> {tgt}")
            continue
        report.add("Dup rows with no target object (kept as their own Duplicate object)", r["Name"])
        o = make_obj(r)
        if o:
            objs[r["Name"]] = o

    relations: list[tuple[str, str, str]] = []
    for obj, tgt, kind in read_links(links_csv, report):
        if obj not in objs or tgt not in objs:
            report.add("links.csv rows naming unknown objects (ignored)", f"{obj};{tgt};{kind}")
            continue
        if kind == "same":
            o = objs.pop(obj)
            objs[tgt].folded.extend([o.row, *o.folded])
        else:
            relations.append((obj, tgt, kind))
    for r in relations:
        if r[0] not in objs or r[1] not in objs:
            report.add("Relations dropped because an end was folded away", ";".join(r))
    relations = [r for r in relations if r[0] in objs and r[1] in objs]

    # Designation claims, then each (catalog, code) goes to its strongest claim.
    claims: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    order = {n: i for i, n in enumerate(objs)}
    for n, o in objs.items():
        for strength, r in [(OWN, o.row)] + [(FOLDED, f) for f in o.folded]:
            # The row's own name and Messier number are its own ids; the rest are cross-ids.
            m = (r.get("M") or "").strip()
            for raw in [r["Name"]] + ([f"M{m}"] if m else []):
                own = parse_designation(raw)
                if own:
                    claims[(own[0].id, own[1])].append((strength, order[n], n))
            # Common names that are really catalog ids ("RCW 130") count as cross-ids too.
            for raw in cross_ids(r) + (r.get("Identifiers") or "").split(",") + (r.get("Common names") or "").split(","):
                d = parse_designation(raw)
                if d:
                    claims[(d[0].id, d[1])].append((CROSS, order[n], n))
    designations: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for (cat, code), cs in claims.items():
        best: dict[str, tuple[int, int, str]] = {}
        for c in cs:
            if c[2] not in best or c < best[c[2]]:
                best[c[2]] = c
        ranked = sorted(best.values())
        winner = ranked[0]
        if len(ranked) > 1:
            prefix = CATALOGS[PRIORITY[cat]].prefix
            losers = ", ".join(c[2] for c in ranked[1:])
            section = (
                "Designations claimed by two objects at equal strength (first kept, review these)"
                if ranked[1][0] == winner[0]
                else "Cross-ids dropped because another object owns that id"
            )
            report.add(section, f"{prefix}{code}: kept on {winner[2]}, dropped from {losers}")
        designations[winner[2]].append((winner[0], cat, code))

    names: dict[str, list[str]] = {}
    for n, o in objs.items():
        seen: set[str] = set()
        out: list[str] = []
        for r in [o.row, *o.folded]:
            for nm in (r.get("Common names") or "").split(","):
                nm = nm.strip()
                if nm and nm.lower() not in seen and parse_designation(nm) is None:
                    seen.add(nm.lower())
                    out.append(nm)
        names[n] = out

    def natural(code: str) -> tuple:
        return tuple(int(t) if t.isdigit() else t for t in re.findall(r"\d+|\D+", code))

    for n in objs:
        designations[n].sort(key=lambda d: (PRIORITY[d[1]], d[0], natural(d[2])))

    return objs, designations, names, relations


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_report(report: Report, path: Path) -> None:
    out: list[str] = []
    for section, lines in report.lines.items():
        out.append(f"== {section} ({len(lines)})")
        out.extend(lines)
        out.append("")
    path.write_text("\n".join(out), encoding="utf-8")


SCHEMA = """
CREATE TABLE catalogs (
  id TEXT PRIMARY KEY,
  name TEXT NOT NULL,
  prefix TEXT NOT NULL,
  priority INTEGER NOT NULL
) WITHOUT ROWID;

CREATE TABLE types (
  code TEXT PRIMARY KEY,
  label TEXT NOT NULL
) WITHOUT ROWID;

CREATE TABLE objects (
  id TEXT PRIMARY KEY,
  type TEXT NOT NULL REFERENCES types(code),
  ra_deg REAL NOT NULL,
  dec_deg REAL NOT NULL,
  constellation TEXT,
  v_mag REAL,
  b_mag REAL,
  surface_brightness REAL,
  major_arcmin REAL,
  minor_arcmin REAL,
  pa_deg REAL
) WITHOUT ROWID;

CREATE TABLE designations (
  catalog TEXT NOT NULL REFERENCES catalogs(id),
  code TEXT NOT NULL,
  object_id TEXT NOT NULL REFERENCES objects(id),
  rank INTEGER NOT NULL,
  PRIMARY KEY (catalog, code)
) WITHOUT ROWID;

CREATE TABLE names (
  object_id TEXT NOT NULL REFERENCES objects(id),
  rank INTEGER NOT NULL,
  name TEXT NOT NULL,
  PRIMARY KEY (object_id, rank)
) WITHOUT ROWID;

CREATE TABLE relations (
  object_id TEXT NOT NULL REFERENCES objects(id),
  related_id TEXT NOT NULL REFERENCES objects(id),
  kind TEXT NOT NULL,
  PRIMARY KEY (object_id, related_id, kind)
) WITHOUT ROWID;

CREATE TABLE meta (
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
) WITHOUT ROWID;
"""


def content_hash(paths: list[Path]) -> str:
    """Identifies the catalog's content: same inputs give the same hash on any machine."""
    h = hashlib.sha256()
    for p in paths:
        h.update(p.name.encode("utf-8") + b"\0")
        h.update(p.read_bytes().replace(b"\r\n", b"\n"))
    return h.hexdigest()


def main(argv: list[str]) -> int:
    repo_root = Path(__file__).resolve().parents[1]
    data = repo_root / "database_files"
    out_dir = repo_root / "dist"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_db = out_dir / "catalog.sqlite"
    out_meta = out_dir / "catalog.meta.json"
    out_report = out_dir / "catalog.report.txt"

    source_sha = os.environ.get("GITHUB_SHA", "").strip() or "unknown"
    built_at = dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    # The app keys result freshness on this, so local and CI builds of the same data match.
    source_hash = content_hash(
        [data / "NGC.csv", data / "addendum.csv", data / "links.csv", Path(__file__).resolve()]
    )

    report = Report()
    objs, designations, names, relations = build(
        data / "NGC.csv", data / "addendum.csv", data / "links.csv", report
    )

    if out_db.exists():
        out_db.unlink()
    con = sqlite3.connect(str(out_db))
    try:
        con.execute("PRAGMA journal_mode=OFF;")
        con.execute("PRAGMA synchronous=OFF;")
        con.executescript(SCHEMA)
        used_types = {o.row["Type"] for o in objs.values()}
        with con:
            con.executemany(
                "INSERT INTO catalogs VALUES (?, ?, ?, ?);",
                [(c.id, c.name, c.prefix, i) for i, c in enumerate(CATALOGS)],
            )
            con.executemany(
                "INSERT INTO types VALUES (?, ?);",
                [(t, TYPES.get(t, t)) for t in sorted(used_types | set(TYPES))],
            )
            con.executemany(
                "INSERT INTO objects VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
                [
                    (
                        n,
                        o.row["Type"],
                        o.ra,
                        o.dec,
                        (o.row.get("Const") or "").strip() or None,
                        fnum(o.row.get("V-Mag")),
                        fnum(o.row.get("B-Mag")),
                        fnum(o.row.get("SurfBr")),
                        fnum(o.row.get("MajAx")),
                        fnum(o.row.get("MinAx")),
                        fnum(o.row.get("PosAng")),
                    )
                    for n, o in objs.items()
                ],
            )
            con.executemany(
                "INSERT INTO designations VALUES (?, ?, ?, ?);",
                [(cat, code, n, i) for n, ds in designations.items() for i, (_, cat, code) in enumerate(ds)],
            )
            con.executemany(
                "INSERT INTO names VALUES (?, ?, ?);",
                [(n, i, nm) for n, ns in names.items() for i, nm in enumerate(ns)],
            )
            con.executemany("INSERT INTO relations VALUES (?, ?, ?);", relations)
            counts = {
                t: con.execute(f"SELECT COUNT(*) FROM {t};").fetchone()[0]
                for t in ("objects", "designations", "names", "relations")
            }
            con.executemany(
                "INSERT INTO meta VALUES (?, ?);",
                [
                    ("schema_version", SCHEMA_VERSION),
                    ("built_at", built_at),
                    ("source_sha", source_sha),
                    ("source_hash", source_hash),
                    *[(f"{t}_count", str(c)) for t, c in counts.items()],
                ],
            )
        con.execute("VACUUM;")
    finally:
        con.close()

    meta = {
        "schema_version": int(SCHEMA_VERSION),
        "built_at": built_at,
        "source_sha": source_sha,
        "source_hash": source_hash,
        "row_count": counts["objects"],
        "sqlite_asset": "catalog.sqlite",
        "sha256": sha256_file(out_db),
    }
    out_meta.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    write_report(report, out_report)
    print(f"Wrote {out_db}: " + ", ".join(f"{c} {t}" for t, c in counts.items()))
    for section, lines in report.lines.items():
        print(f"  {section}: {len(lines)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
