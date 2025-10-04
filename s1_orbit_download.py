#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sentinel-1 Precise Orbit (AUX_POEORB) downloader — JSON config version

Source directory: https://s1qc.asf.alaska.edu/aux_poeorb/

- Fetch directory listing
- Match filenames with regex
- Skip download if already exists locally
- Default output directory: S1_orbits/
- No login required; if in the future host requires auth, use ~/.netrc

Extra features:
- Select platform: S1A / S1B / BOTH
- Filter by validity period from filename (_VYYYYMMDDThhmmss_YYYYMMDDThhmmss)
- Retry logic and optional tqdm progress bar

Config:
- Reads all parameters from s1_config.json
- No CLI arguments, everything is controlled via JSON

Usage:
  python s1_orbit_download.py
"""

import os
import re
import sys
import json
from pathlib import Path
from datetime import datetime, timezone
from time import sleep

BASE_URL = "https://s1qc.asf.alaska.edu/aux_poeorb/"
REGEX = re.compile(
    r"(S1[A-B]_OPER_AUX_POEORB_OPOD_\d{8}T\d{6}_V\d{8}T\d{6}_\d{8}T\d{6}\.EOF)"
)

# ------------------ Config Loader ------------------
def load_config(cfg_path="s1_config.json"):
    p = Path(cfg_path)
    if not p.exists():
        print(f"Config file {cfg_path} not found", file=sys.stderr)
        sys.exit(1)
    with p.open("r", encoding="utf-8") as f:
        return json.load(f) or {}

# ------------------ Utilities ------------------
def dt_utc(s):
    """Parse string like 20240327T100224 into datetime UTC"""
    return datetime.strptime(s, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)

def valid_window_from_name(name):
    """Extract validity window from filename."""
    m = re.search(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$", name)
    if not m:
        return None, None
    return dt_utc(m.group(1)), dt_utc(m.group(2))

def filter_by_platform(names, platform):
    if platform.upper() == "BOTH":
        return names
    return [n for n in names if n.startswith(platform.upper())]

def filter_by_date(names, start_str, end_str):
    if not start_str and not end_str:
        return names
    start = dt_utc(start_str.replace("-", "") + "T000000") if start_str else None
    end   = dt_utc(end_str.replace("-", "") + "T235959")   if end_str   else None
    kept = []
    for n in names:
        vstart, vend = valid_window_from_name(n)
        if vstart is None:
            continue
        ok = True
        if start and vend < start:
            ok = False
        if end and vstart > end:
            ok = False
        if ok:
            kept.append(n)
    return kept

def get_session():
    """Return requests.Session with headers."""
    try:
        import requests
    except ImportError:
        print("requests library required: conda install requests OR pip install requests",
              file=sys.stderr)
        sys.exit(1)
    s = requests.Session()
    s.headers.update({"User-Agent":"s1-poeorb-downloader/1.0"})
    return s

def fetch_listing(session):
    r = session.get(BASE_URL, timeout=30)
    r.raise_for_status()
    return r.text

def parse_names(html):
    """Parse AUX_POEORB filenames from HTML index page."""
    return sorted(set(REGEX.findall(html)))

def download_one(session, name, outdir, retry=3):
    """Download one orbit file with retry logic and optional tqdm progress bar."""
    from requests.exceptions import RequestException
    url = BASE_URL + name
    dst = Path(outdir) / name
    if dst.exists():
        return "skip"
    for attempt in range(1, retry+1):
        try:
            with session.get(url, stream=True, timeout=60) as resp:
                resp.raise_for_status()
                total = int(resp.headers.get("Content-Length", 0))
                # tqdm optional
                try:
                    from tqdm import tqdm
                    bar = tqdm(total=total, unit="B", unit_scale=True, desc=name, leave=False)
                    with open(dst, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=1024*1024):
                            if chunk:
                                f.write(chunk)
                                bar.update(len(chunk))
                    bar.close()
                except Exception:
                    with open(dst, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=1024*1024):
                            if chunk:
                                f.write(chunk)
            return "ok"
        except RequestException as e:
            if attempt < retry:
                sleep(2*attempt)
            else:
                if dst.exists():
                    try: dst.unlink()
                    except Exception: pass
                return f"fail: {e}"

# ------------------ Main ------------------
def main():
    cfg = load_config("s1_config.json")
    orbit_cfg = cfg.get("orbit", {})

    outdir = Path(orbit_cfg.get("dir", "S1_orbits"))
    platform = orbit_cfg.get("platform", "BOTH")
    start = orbit_cfg.get("start")
    end = orbit_cfg.get("end")
    retry = int(orbit_cfg.get("retry", 3))

    outdir.mkdir(parents=True, exist_ok=True)

    s = get_session()
    print(f"Fetching listing from {BASE_URL} ...")
    html = fetch_listing(s)
    names = parse_names(html)
    if not names:
        print("Empty listing, server may be unavailable. Try again later.")
        sys.exit(2)

    names = filter_by_platform(names, platform)
    names = filter_by_date(names, start, end)
    print(f"There are {len(names)} files (platform={platform}, start={start}, end={end})")

    cnt_ok = cnt_skip = cnt_fail = 0
    for i, name in enumerate(names, 1):
        status = download_one(s, name, outdir, retry=retry)
        if status == "ok":
            cnt_ok += 1
            print(f"[{i}/{len(names)}] downloaded {name}")
        elif status == "skip":
            cnt_skip += 1
            print(f"[{i}/{len(names)}] already exists, skip: {name}")
        else:
            cnt_fail += 1
            print(f"[{i}/{len(names)}] failed: {name}\n    {status}")

    print(f"\nSummary: downloaded {cnt_ok} ; skipped {cnt_skip} ; failed {cnt_fail}")

if __name__ == "__main__":
    main()
