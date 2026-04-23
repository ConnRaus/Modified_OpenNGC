#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import sqlite3
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


RE_NGC_IC = re.compile(r"^(NGC|IC)\d+$")


def parse_coord(coord: str, is_ra: bool) -> float:
    s = (coord or "").strip()
    if not s:
        return float("nan")
    neg = s.startswith("-")
    s = s.lstrip("+-")
    parts = s.split(":")
    try:
        if len(parts) == 3:
            h = float(parts[0])
            m = float(parts[1])
            sec = float(parts[2])
            dec = h + m / 60.0 + sec / 3600.0
            if is_ra:
                dec *= 15.0
            if neg:
                dec *= -1.0
            return dec
        return float(s) * (-1.0 if neg else 1.0)
    except Exception:
        return float("nan")


def normalize_messier_id(name: str) -> str:
    # Mirrors StarTeller's `normalizeMessierId` behavior for inputs like "M 31", "M31".
    s = (name or "").strip()
    m = re.match(r"^M\s*(\d+)$", s, re.IGNORECASE)
    if m:
        return f"M{int(m.group(1))}"
    return s


def clean_name(name: str) -> str:
    s = (name or "").strip()
    m = re.match(r"^(NGC|IC)(\d+)$", s)
    if m:
        return f"{m.group(1)} {int(m.group(2))}"
    m = re.match(r"^C(\d+)$", s)
    if m:
        return f"C {int(m.group(1))}"
    m = re.match(r"^B(\d+)$", s)
    if m:
        return f"B {int(m.group(1))}"
    m = re.match(r"^H(\d+)$", s)
    if m:
        return f"H {int(m.group(1))}"
    m = re.match(r"^Mel(\d+)$", s)
    if m:
        return f"Mel {int(m.group(1))}"
    m = re.match(r"^UGC(\d+)$", s)
    if m:
        return f"UGC {int(m.group(1))}"
    m = re.match(r"^PGC(\d+)$", s)
    if m:
        return f"PGC {int(m.group(1))}"
    m = re.match(r"^MWSC(\d+)$", s)
    if m:
        return f"MWSC {int(m.group(1))}"
    m = re.match(r"^HCG(\d+)$", s)
    if m:
        return f"HCG {int(m.group(1))}"
    m = re.match(r"^ESO(\d+)-(\d+)$", s)
    if m:
        return f"ESO {m.group(1)}-{m.group(2)}"
    m = re.match(r"^Cl(\d+)$", s)
    if m:
        return f"Cl {int(m.group(1))}"
    m = re.match(r"^M(\d+)$", s)
    if m:
        return f"M {int(m.group(1))}"
    return s


def format_messier(m_val: str) -> str:
    s = (m_val or "").strip()
    if not s:
        return ""
    try:
        n = int(float(s))
        return f"M{n}"
    except Exception:
        return ""


TYPE_EXPANSIONS: dict[str, str] = {
    "G": "Galaxy",
    "SNR": "Supernova remnant",
    "GCl": "Globular cluster",
    "GCI": "Globular cluster",
    "OCl": "Open cluster",
    "Neb": "Nebula",
    "HII": "HII region",
    "PN": "Planetary nebula",
    "RfN": "Reflection nebula",
    "DrkN": "Dark nebula",
    "**": "Double star",
    "*": "Star",
    "*Ass": "Stellar association",
    "GPair": "Galaxy pair",
    "GGroup": "Galaxy group",
    "GTrpl": "Galaxy triplet",
    "EmN": "Emission nebula",
    "Nova": "Nova",
    "Dup": "Duplicate object",
    "Other": "Other object",
    "Cl+N": "Cluster with nebula",
    "NonEx": "Non-existent object",
}


@dataclass(frozen=True)
class CatalogRow:
    Object: str
    Catalog_Name: str
    Name: str
    Right_Ascension: float
    Declination: float
    Type: str
    Magnitude: float
    Common_Name: str
    Messier: str
    Constellation: str
    V_Mag: float
    SurfBr: float
    Major_Axis_arcmin: float
    Minor_Axis_arcmin: float
    Position_Angle_deg: float


def fnum(s: str) -> float:
    s = (s or "").strip()
    if not s:
        return float("nan")
    try:
        return float(s)
    except Exception:
        return float("nan")


def iter_catalog_rows(csv_path: Path, enforce_ngc_ic: bool) -> Iterable[CatalogRow]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=";")
        for r in reader:
            name_raw = (r.get("Name") or "").strip()
            if not name_raw:
                continue
            if enforce_ngc_ic and not RE_NGC_IC.match(name_raw):
                continue

            ra_str = (r.get("RA") or "").strip()
            dec_str = (r.get("Dec") or "").strip()
            if not ra_str or not dec_str:
                continue

            ra = parse_coord(ra_str, True)
            de = parse_coord(dec_str, False)
            if not (ra == ra and de == de):
                continue
            if ra < 0 or ra > 360 or de < -90 or de > 90:
                continue

            obj_id = normalize_messier_id(name_raw)
            catalog_name = clean_name(name_raw)
            common_raw = (r.get("Common names") or "").strip()
            display_name = common_raw if common_raw else catalog_name

            t0 = (r.get("Type") or "").strip()
            t = TYPE_EXPANSIONS.get(t0, t0)

            v_mag = fnum(r.get("V-Mag") or "")
            b_mag = fnum(r.get("B-Mag") or "")
            magnitude = v_mag if (v_mag == v_mag) else b_mag

            yield CatalogRow(
                Object=obj_id,
                Catalog_Name=catalog_name,
                Name=display_name,
                Right_Ascension=ra,
                Declination=de,
                Type=t,
                Magnitude=magnitude if magnitude == magnitude else float("nan"),
                Common_Name=common_raw,
                Messier=format_messier(r.get("M") or ""),
                Constellation=(r.get("Const") or "").strip(),
                V_Mag=v_mag,
                SurfBr=fnum(r.get("SurfBr") or ""),
                Major_Axis_arcmin=fnum(r.get("MajAx") or ""),
                Minor_Axis_arcmin=fnum(r.get("MinAx") or ""),
                Position_Angle_deg=fnum(r.get("PosAng") or ""),
            )


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str]) -> int:
    repo_root = Path(__file__).resolve().parents[1]
    ngc_csv = repo_root / "database_files" / "NGC.csv"
    add_csv = repo_root / "database_files" / "addendum.csv"
    out_dir = repo_root / "dist"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_db = out_dir / "catalog.sqlite"
    out_meta = out_dir / "catalog.meta.json"

    source_sha = os.environ.get("GITHUB_SHA", "").strip() or "unknown"
    built_at = dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )

    if out_db.exists():
        out_db.unlink()

    con = sqlite3.connect(str(out_db))
    try:
        con.execute("PRAGMA journal_mode=OFF;")
        con.execute("PRAGMA synchronous=OFF;")
        con.execute("PRAGMA temp_store=MEMORY;")
        con.execute("PRAGMA locking_mode=EXCLUSIVE;")

        con.execute(
            """
            CREATE TABLE objects (
              Object TEXT NOT NULL,
              Catalog_Name TEXT NOT NULL,
              Name TEXT NOT NULL,
              Right_Ascension REAL NOT NULL,
              Declination REAL NOT NULL,
              Type TEXT NOT NULL,
              Magnitude REAL,
              Common_Name TEXT,
              Messier TEXT,
              Constellation TEXT,
              V_Mag REAL,
              SurfBr REAL,
              Major_Axis_arcmin REAL,
              Minor_Axis_arcmin REAL,
              Position_Angle_deg REAL
            );
            """
        )
        con.execute("CREATE INDEX idx_objects_object ON objects(Object);")
        con.execute("CREATE INDEX idx_objects_messier ON objects(Messier);")

        insert_sql = """
          INSERT INTO objects (
            Object, Catalog_Name, Name, Right_Ascension, Declination, Type, Magnitude,
            Common_Name, Messier, Constellation, V_Mag, SurfBr, Major_Axis_arcmin,
            Minor_Axis_arcmin, Position_Angle_deg
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """

        def batched(it: Iterable[CatalogRow], n: int = 2000) -> Iterable[list[CatalogRow]]:
            batch: list[CatalogRow] = []
            for x in it:
                batch.append(x)
                if len(batch) >= n:
                    yield batch
                    batch = []
            if batch:
                yield batch

        total = 0
        with con:
            for batch in batched(iter_catalog_rows(ngc_csv, enforce_ngc_ic=True)):
                con.executemany(insert_sql, [tuple(vars(r).values()) for r in batch])
                total += len(batch)
            for batch in batched(iter_catalog_rows(add_csv, enforce_ngc_ic=False)):
                con.executemany(insert_sql, [tuple(vars(r).values()) for r in batch])
                total += len(batch)

        con.execute(
            "CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);"
        )
        with con:
            con.executemany(
                "INSERT INTO meta(key, value) VALUES (?, ?);",
                [
                    ("built_at", built_at),
                    ("source_sha", source_sha),
                    ("row_count", str(total)),
                ],
            )

        con.execute("VACUUM;")
    finally:
        con.close()

    meta = {
        "built_at": built_at,
        "source_sha": source_sha,
        "row_count": None,
        "sqlite_asset": "catalog.sqlite",
        "sha256": sha256_file(out_db),
    }
    # read row count back from db for accuracy
    con2 = sqlite3.connect(str(out_db))
    try:
        row_count = con2.execute("SELECT COUNT(*) FROM objects;").fetchone()[0]
        meta["row_count"] = int(row_count)
    finally:
        con2.close()

    out_meta.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out_db} ({meta['row_count']} rows) and {out_meta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

