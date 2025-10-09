#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- User configurable parameters are stored in s1_config.json
- Program reads the config file as default values; CLI args can override them.
- Critical parameters like 'event_date' and 'default_aoi_wkt' are mandatory.
- This version handles multiple frames per acquisition for full AOI coverage.
Usage:
  python s1_data_search.py --config s1_config.json
"""

import json, os, argparse, re, sys
from pathlib import Path
from datetime import datetime, timedelta, timezone
from collections import defaultdict
import asf_search as asf

# ------------------ Load JSON config ------------------
def load_config(cfg_path: str):
    """
    Load JSON config file.
    If file does not exist, return empty dict.
    """
    p = Path(cfg_path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        return json.load(f) or {}

def cfg_get(dct, path, default=None):
    """
    Nested dictionary reader:
    e.g. cfg_get(cfg, "search.flight", "ASC")
    """
    cur = dct
    for k in path.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

# ---------- Utilities ----------
def iso2dt(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z","+00:00")).astimezone(timezone.utc)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def get_track_id(props: dict):
    return props.get("relativeOrbit") or props.get("pathNumber")

def get_zipname(p: dict) -> str:
    return (
        p.get("granuleName")
        or p.get("fileName")
        or (p.get("fileID", "").split("/")[-1] if p.get("fileID") else None)
        or (p.get("url", "").split("/")[-1] if p.get("url") else None)
        or ""
    )

def to_safe(name: str) -> str:
    if not name:
        return ""
    if name.endswith(".SAFE"):
        return name
    if name.endswith(".zip"):
        return name[:-4] + ".SAFE"
    return name + ".SAFE"

def get_bbox_from_wkt(wkt_string: str):
    """Extracts the bounding box from a WKT POLYGON string."""
    try:
        # Find all numbers within the first parenthesis
        coords_str = re.findall(r"\(([^()]+)\)", wkt_string)
        if not coords_str:
            return None
        
        points = [float(p) for p in re.findall(r"[-+]?\d*\.?\d+", coords_str[0])]
        
        longitudes = points[0::2] # Even indices
        latitudes = points[1::2]  # Odd indices
        
        if not longitudes or not latitudes:
            return None

        return (min(longitudes), min(latitudes), max(longitudes), max(latitudes))
    except (IndexError, ValueError):
        return None

def pick_best_pairs(results, event_dt):
    """
    Revised logic:
    1. Group images by track, then by date.
    2. For each track, find the pre/post dates closest to the event.
    3. Select the track with the shortest temporal baseline.
    4. Return all images for that best pre-date and post-date pair.
    """
    if not results:
        return [], {}

    tracks = defaultdict(list)
    for r in results:
        tid = get_track_id(r.properties)
        if tid is not None:
            tracks[tid].append(r)

    candidates = []
    for tid, track_images in tracks.items():
        acquisitions = defaultdict(list)
        for img in track_images:
            acq_date = iso2dt(img.properties["startTime"]).date()
            acquisitions[acq_date].append(img)

        pre_dates = sorted([d for d in acquisitions if d < event_dt.date()], reverse=True)
        post_dates = sorted([d for d in acquisitions if d >= event_dt.date()])

        if not pre_dates or not post_dates:
            continue

        best_pre_date = pre_dates[0]
        best_post_date = post_dates[0]
        
        pre_group = sorted(acquisitions[best_pre_date], key=lambda r: r.properties["startTime"])
        post_group = sorted(acquisitions[best_post_date], key=lambda r: r.properties["startTime"])
        
        temporal_baseline = (best_post_date - best_pre_date).total_seconds()
        
        candidates.append({
            "tid": tid,
            "pre": pre_group,
            "post": post_group,
            "baseline": temporal_baseline,
            "orbit_images": track_images
        })

    if not candidates:
        all_tracks_sorted = {tid: sorted(g, key=lambda r: r.properties["startTime"]) for tid, g in tracks.items()}
        return [], all_tracks_sorted

    best_candidate = min(candidates, key=lambda x: x["baseline"])
    
    return [best_candidate], {}

def make_topsapp_xml(pre_zips, post_zips, dem_path="~/dem/dem.tif", pol="VV",
                     swaths=(1,2,3), az_looks=7, rg_looks=2, filter_strength=0.5, do_unwrap=False,
                     region_of_interest=None):
    """
    Generate topsApp.xml template, handling lists of SAFE files and region of interest.
    """
    pre_safes = [f"./REF/{to_safe(name)}" for name in pre_zips]
    post_safes = [f"./SEC/{to_safe(name)}" for name in post_zips]

    pre_safe_str = json.dumps(pre_safes) if len(pre_safes) > 1 else f'"{pre_safes[0]}"'
    post_safe_str = json.dumps(post_safes) if len(post_safes) > 1 else f'"{post_safes[0]}"'

    roi_property_line = ""
    if region_of_interest:
        roi_property_line = f'    <property name="regionOfInterest">{str(region_of_interest)}</property>'

    return f"""<topsApp>
  <component name="topsApp">
    <property name="sensor name">SENTINEL1</property>
    <property name="swaths">{list(swaths)}</property>
    <property name="azimuth looks">{az_looks}</property>
    <property name="range looks">{rg_looks}</property>
    <property name="filter strength">{filter_strength}</property>
{roi_property_line}
    <property name="do unwrap">{str(do_unwrap)}</property>
    <component name="reference">
      <property name="safe">{pre_safe_str}</property>
      <property name="orbit directory">./orbits</property>
      <property name="output directory">./ref</property>
    </component>
    <component name="secondary">
      <property name="safe">{post_safe_str}</property>
      <property name="orbit directory">./orbits</property>
      <property name="output directory">./sec</property>
    </component>
  </component>
</topsApp>
"""

def write_footprints_png(results, out_png="footprints.png"):
    """
    Optional: plot footprints using pygmt.
    """
    try:
        import pygmt
    except Exception:
        print("pygmt not available, skip footprint plot.")
        return
    def wesn_of(item):
        xy = item.geometry["coordinates"][0]
        xs = [p[0] for p in xy]; ys = [p[1] for p in xy]
        return [min(xs), max(xs), min(ys), max(ys)]
    def merge(a,b):
        return [min(a[0],b[0]), max(a[1],b[1]), min(a[2],b[2]), max(a[3],b[3])]
    wesn = wesn_of(results[0])
    for it in results[1:]:
        wesn = merge(wesn, wesn_of(it))
    tmp = "_fp.txt"
    with open(tmp,"w") as f:
        for im in results:
            f.write(">\n")
            for x,y in im.geometry["coordinates"][0]:
                f.write(f"{x} {y}\n")
    fig = pygmt.Figure()
    fig.basemap(region=wesn, projection="M6i", frame='+t"Sentinel-1 SLC footprints"')
    fig.coast(region=wesn, borders="1", shorelines="0.5p,black", water="white")
    fig.plot(data=tmp, pen="0.4p,red")
    os.remove(tmp)
    fig.savefig(out_png)
    print(f"Footprint saved to {out_png}")

# ------------------ Main ------------------
def main():
    ap0 = argparse.ArgumentParser(add_help=False)
    ap0.add_argument("--config", default="s1_config.json", help="JSON config file")
    args0, _ = ap0.parse_known_args()
    cfg = load_config(args0.config)

    event_date_str = cfg_get(cfg, "event_date")
    if event_date_str is None:
        print(f"Error: Missing required parameter 'event_date' in config file '{args0.config}'.")
        sys.exit(1)
    default_aoi_wkt = cfg_get(cfg, "default_aoi_wkt")
    if default_aoi_wkt is None:
        print(f"Error: Missing required parameter 'default_aoi_wkt' in config file '{args0.config}'.")
        sys.exit(1)

    defcfg = {
        "flight": cfg_get(cfg, "search.flight", "ASC"),
        "pol": cfg_get(cfg, "search.pol", "VV"),
        "days": int(cfg_get(cfg, "search.days", 21)),
        "max": int(cfg_get(cfg, "search.max_results", 1000)),
        "plot": bool(cfg_get(cfg, "search.plot", False)),
        "aoi": cfg_get(cfg, "search.aoi") or default_aoi_wkt,
        "dem": cfg_get(cfg, "topsapp.dem", "~/dem/dem.tif"),
        "swaths": tuple(cfg_get(cfg, "topsapp.swaths", [1,2,3])),
        "az_looks": int(cfg_get(cfg, "topsapp.az_looks", 7)),
        "rg_looks": int(cfg_get(cfg, "topsapp.rg_looks", 2)),
        "filter_strength": float(cfg_get(cfg, "topsapp.filter_strength", 0.5)),
        "do_unwrap": bool(cfg_get(cfg, "topsapp.do_unwrap", False)),
    }
    ap = argparse.ArgumentParser(description="DInSAR pair picker for Sentinel-1 SLC/IW data", parents=[ap0])
    ap.add_argument("--flight", choices=["ASC","DES","BOTH"], default=defcfg["flight"], help="Orbit direction")
    ap.add_argument("--pol", default=defcfg["pol"], help="Polarization (default VV)")
    ap.add_argument("--days", type=int, default=defcfg["days"], help="Search days before/after event")
    ap.add_argument("--max", type=int, default=defcfg["max"], help="Max results")
    ap.add_argument("--plot", action="store_true" if not defcfg["plot"] else "store_false")
    ap.add_argument("--dem", default=defcfg["dem"], help="DEM path for topsApp.xml")
    ap.add_argument("--aoi", default=defcfg["aoi"], help="AOI polygon WKT")
    ap.add_argument("--swaths", type=str, default=",".join(map(str, defcfg["swaths"])), help="e.g. 1,2,3")
    ap.add_argument("--az-looks", type=int, default=defcfg["az_looks"])
    ap.add_argument("--rg-looks", type=int, default=defcfg["rg_looks"])
    ap.add_argument("--filter-strength", type=float, default=defcfg["filter_strength"])
    ap.add_argument("--do-unwrap", action="store_true" if not defcfg["do_unwrap"] else "store_false")

    args = ap.parse_args()
    swaths = tuple(int(x) for x in re.split(r"[,\s]+", args.swaths.strip()) if x)
    event_dt = datetime.fromisoformat(event_date_str).replace(tzinfo=timezone.utc)
    start = (event_dt - timedelta(days=args.days)).strftime("%Y-%m-%d")
    end   = (event_dt + timedelta(days=args.days)).strftime("%Y-%m-%d")

    fd = "ASCENDING" if args.flight == "ASC" else "DESCENDING" if args.flight == "DES" else None
    opts = dict(platform=asf.PLATFORM.SENTINEL1, processingLevel="SLC", beamMode="IW",
                start=f"{start}T00:00:00Z", end=f"{end}T23:59:59Z", maxResults=args.max)
    if fd: opts["flightDirection"] = fd

    print(f"[Config] {args0.config}")
    print(f"Search window: {start} ~ {end}; flight: {args.flight}; pol: {args.pol}")
    res = asf.geo_search(intersectsWith=args.aoi, **opts)
    if not res:
        print("No results found, try larger AOI or wider time window.")
        return

    def has_pol(pval, want):
        if pval is None: return False
        s = str(pval).upper().replace(" ","").replace(",","/").replace("+","/")
        parts = [x for x in s.split("/") if x]
        return want.upper() in parts or want.upper() == s

    res = [r for r in res if has_pol(r.properties.get("polarization"), args.pol)]
    if not res:
        print("No results with requested polarization. Try adjusting --pol.")
        return

    pairs, debug = pick_best_pairs(res, event_dt)
    if not pairs:
        print("No matching pre/post pairs in the same orbit. Orbit summary:")
        for tid, lst in sorted(debug.items(), key=lambda kv: kv[0]):
            times = [iso2dt(x.properties["startTime"]) for x in lst]
            print(f"  Orbit {tid}: {len(lst)} images, {times[0].date()} ~ {times[-1].date()}")
        print("Suggestion: try --flight DES, or increase --days (e.g. 60~400), or expand --aoi.")
        return

    best_pair = pairs[0]
    tid, pre_group, post_group, orbit_images = best_pair["tid"], best_pair["pre"], best_pair["post"], best_pair["orbit_images"]
    pre_names = [get_zipname(p.properties) for p in pre_group]
    post_names = [get_zipname(p.properties) for p in post_group]

    print(f"Selected orbit = {tid}")
    print(f"PRE  ({pre_group[0].properties['startTime'].split('T')[0]}):")
    for name in pre_names: print(f"  {name}")
    print(f"POST ({post_group[0].properties['startTime'].split('T')[0]}):")
    for name in post_names: print(f"  {name}")

    roi_bbox = get_bbox_from_wkt(default_aoi_wkt)
    
    roi_for_xml = None
    if roi_bbox:
        min_lon, min_lat, max_lon, max_lat = roi_bbox
        roi_for_xml = [min_lat, max_lat, min_lon, max_lon]
        print(f"Region of Interest from AOI (S,N,W,E): {roi_for_xml}")

    outdir = Path("S1_pick"); ensure_dir(outdir)

    def slim(x):
        p = x.properties
        return {"name": get_zipname(p), "startTime": p.get("startTime"), "url": p.get("url"),
                "relativeOrbit": p.get("relativeOrbit"), "pathNumber": p.get("pathNumber"),
                "flightDirection": p.get("flightDirection")}

    with (outdir/"pairs.json").open("w") as f:
        json.dump([
            {"track": tid, "pre": [slim(x) for x in pre_group], "post": [slim(x) for x in post_group],
             "orbit_images": [slim(x) for x in orbit_images]}
        ], f, indent=2)

    with (outdir/"download.txt").open("w") as f:
        for img in orbit_images:
            url = img.properties.get("url", "")
            if url: f.write(url + "\n")

    xml_txt = make_topsapp_xml(
        pre_names, post_names, dem_path=args.dem, pol=args.pol, swaths=swaths,
        az_looks=args.az_looks, rg_looks=args.rg_looks,
        filter_strength=args.filter_strength, do_unwrap=args.do_unwrap,
        region_of_interest=roi_for_xml)
    
    Path("topsApp.xml").write_text(xml_txt)
    print("Wrote topsApp.xml (ready to run ISCE2 topsApp)")

    if args.plot or defcfg["plot"]:
        write_footprints_png(res, out_png="footprints.png")

    print("\nNext Steps:")
    print("  # Step 1: Setup Earthdata Login credentials")
    print("  echo 'machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD>' > ~/.netrc")
    print("  chmod 600 ~/.netrc")
    print("  # Step 2: Download Sentinel-1 data")
    print("  aria2c -i S1_pick/download.txt -x 8 -s 8 -k 1M -d raw")
    print("  mkdir -p REF SEC orbits dem")
    print("  # Unzip all PRE images into REF/ and all POST images into SEC/")
    print("  # e.g., for f in raw/S1A_..._PRE_DATE_...zip; do unzip \"$f\" -d REF/; done")
    print("  # e.g., for f in raw/S1A_..._POST_DATE_...zip; do unzip \"$f\" -d SEC/; done")
    print("  python s1_orbit_download.py # (Or a similar script to get orbits for all images)")
    print("  # Place corresponding *.EOF into orbits/, and your DEM into", args.dem)
    print(" You can run step by step: python -m isce.applications.topsApp --start=preprocess --end=preprocess topsApp.xml")
    print(" Or just run all process: python -m isce.applications.topsApp topsApp.xml --steps")

if __name__ == "__main__":
    main()