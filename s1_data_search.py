#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
- User configurable parameters are stored in s1_config.json
- Program reads the config file as default values; CLI args can override them.
Usage:
  python s1_data_search.py --config s1_config.json
"""

import json, os, argparse, re
from pathlib import Path
from datetime import datetime, timedelta, timezone
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

# ------------------ Default constants ------------------
EVENT_DATE = "2024-04-03"
DEFAULT_AOI_WKT = "POLYGON((121.2 23.7, 121.8 23.7, 121.8 24.3, 121.2 24.3, 121.2 23.7))"

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

def pick_best_pairs(results, event_dt):
    """
    Group results by relative orbit; 
    pick one pre-event and one post-event image that overlap.
    """
    from collections import defaultdict
    groups = defaultdict(list)
    for r in results:
        tid = get_track_id(r.properties)
        if tid is not None:
            groups[tid].append(r)
    debug = {tid: sorted(g, key=lambda rr: iso2dt(rr.properties["startTime"])) for tid, g in groups.items()}
    candidates = []
    for tid, lst in debug.items():
        pre  = [rr for rr in lst if iso2dt(rr.properties["startTime"]) < event_dt]
        post = [rr for rr in lst if iso2dt(rr.properties["startTime"]) >= event_dt]
        if not pre or not post:
            continue
        pre_top  = sorted(pre,  key=lambda rr: abs((event_dt - iso2dt(rr.properties["startTime"])).total_seconds()))[:4]
        post_top = sorted(post, key=lambda rr: abs((iso2dt(rr.properties["startTime"]) - event_dt).total_seconds()))[:4]
        best = None
        for a in pre_top:
            for b in post_top:
                if not bbox_intersect(bbox_of(a), bbox_of(b), min_overlap_frac=0.02):
                    continue
                same_sf = same_slice_or_frame(a, b)
                dt = abs((iso2dt(b.properties["startTime"]) - iso2dt(a.properties["startTime"])).total_seconds())
                score = (0 if same_sf is True else (1 if same_sf is None else 2), dt)
                if (best is None) or (score < best[0]):
                    best = (score, (tid, a, b))
        if best is not None:
            candidates.append(best[1])
    candidates.sort(key=lambda t: (iso2dt(t[2].properties["startTime"]) - iso2dt(t[1].properties["startTime"])).total_seconds())
    return candidates, debug

def make_topsapp_xml(pre_zip, post_zip, dem_path="~/dem/dem.tif", pol="VV",
                     swaths=(1,2,3), az_looks=7, rg_looks=2, filter_strength=0.5, do_unwrap=False):
    """
    Generate topsApp.xml template.
    """
    pre_safe  = to_safe(pre_zip)
    post_safe = to_safe(post_zip)
    return f"""<topsApp>
  <component name="topsApp">
    <property name="sensor name">SENTINEL1</property>
    <property name="swaths">{list(swaths)}</property>
    <property name="azimuth looks">{az_looks}</property>
    <property name="range looks">{rg_looks}</property>
    <property name="filter strength">{filter_strength}</property>
    <property name="do unwrap">{str(do_unwrap)}</property>
    <component name="reference">
      <property name="safe">./REF/{pre_safe}</property>
      <property name="orbit directory">./orbits</property>
      <property name="output directory">./ref</property>
    </component>
    <component name="secondary">
      <property name="safe">./SEC/{post_safe}</property>
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

# -------------- bbox and slice/frame check --------------
def bbox_of(res):
    ring = res.geometry["coordinates"][0]
    xs = [p[0] for p in ring]; ys = [p[1] for p in ring]
    return (min(xs), min(ys), max(xs), max(ys))

def bbox_intersect(b1, b2, min_overlap_frac=0.02):
    minx1, miny1, maxx1, maxy1 = b1
    minx2, miny2, maxx2, maxy2 = b2
    ixmin, iymin = max(minx1, minx2), max(miny1, miny2)
    ixmax, iymax = min(maxx1, maxx2), min(maxy1, maxy2)
    if ixmin >= ixmax or iymin >= iymax:
        return False
    iarea = (ixmax - ixmin) * (iymax - iymin)
    a1 = (maxx1 - minx1) * (maxy1 - miny1)
    a2 = (maxx2 - minx2) * (maxy2 - miny2)
    ref = max(min(a1, a2), 1e-9)
    return (iarea / ref) >= min_overlap_frac

def same_slice_or_frame(a, b):
    pa, pb = a.properties, b.properties
    keys = ["sliceNumber", "frameNumber", "frame", "sceneFrameNumber"]
    for k in keys:
        va, vb = pa.get(k), pb.get(k)
        if (va is not None) and (vb is not None) and (va == vb):
            return True
    return None

# ------------------ Main ------------------
def main():
    # Parse --config first
    ap0 = argparse.ArgumentParser(add_help=False)
    ap0.add_argument("--config", default="s1_config.json", help="JSON config file")
    args0, _ = ap0.parse_known_args()

    cfg = load_config(args0.config)

    global EVENT_DATE, DEFAULT_AOI_WKT
    EVENT_DATE = cfg_get(cfg, "event_date", EVENT_DATE)
    DEFAULT_AOI_WKT = cfg_get(cfg, "default_aoi_wkt", DEFAULT_AOI_WKT)

    # Default values from config
    defcfg = {
        "flight": cfg_get(cfg, "search.flight", "ASC"),
        "pol": cfg_get(cfg, "search.pol", "VV"),
        "days": int(cfg_get(cfg, "search.days", 21)),
        "max": int(cfg_get(cfg, "search.max_results", 1000)),
        "plot": bool(cfg_get(cfg, "search.plot", False)),
        "aoi": cfg_get(cfg, "search.aoi", None) or DEFAULT_AOI_WKT,
        "dem": cfg_get(cfg, "topsapp.dem", "~/dem/dem.tif"),
        "swaths": tuple(cfg_get(cfg, "topsapp.swaths", [1,2,3])),
        "az_looks": int(cfg_get(cfg, "topsapp.az_looks", 7)),
        "rg_looks": int(cfg_get(cfg, "topsapp.rg_looks", 2)),
        "filter_strength": float(cfg_get(cfg, "topsapp.filter_strength", 0.5)),
        "do_unwrap": bool(cfg_get(cfg, "topsapp.do_unwrap", False)),
    }

    # argparse with defaults from config
    ap = argparse.ArgumentParser(
        description="0403 Hualien DInSAR pair picker (S1 SLC / IW)",
        parents=[ap0]
    )
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

    event_dt = datetime.fromisoformat(EVENT_DATE).replace(tzinfo=timezone.utc)
    start = (event_dt - timedelta(days=args.days)).strftime("%Y-%m-%d")
    end   = (event_dt + timedelta(days=args.days)).strftime("%Y-%m-%d")

    fd = None
    if args.flight == "ASC": fd = "ASCENDING"
    if args.flight == "DES": fd = "DESCENDING"

    opts = dict(
        platform=asf.PLATFORM.SENTINEL1,
        processingLevel="SLC",
        beamMode="IW",
        start=f"{start}T00:00:00Z",
        end=f"{end}T23:59:59Z",
        maxResults=args.max,
    )
    if fd: opts["flightDirection"] = fd

    print(f"[Config] {args0.config}")
    print(f"Search window: {start} ~ {end}; flight: {args.flight}; pol: {args.pol}")
    res = asf.geo_search(intersectsWith=args.aoi, **opts)
    if len(res)==0:
        print("No results found, try larger AOI or wider time window.")
        return

    def has_pol(pval, want):
        if pval is None: return False
        s = str(pval).upper().replace(" ","").replace(",","/").replace("+","/")
        parts = [x for x in s.split("/") if x]
        return want.upper() in parts or want.upper()==s

    res = [r for r in res if has_pol(r.properties.get("polarization"), args.pol)]
    if len(res)==0:
        print("No results with requested polarization. Try adjusting --pol.")
        return

    pairs, debug = pick_best_pairs(res, event_dt)
    if len(pairs)==0:
        print("No matching pre/post pairs in same orbit. Orbit summary:")
        for tid, lst in sorted(debug.items(), key=lambda kv: kv[0]):
            times = [iso2dt(x.properties["startTime"]) for x in lst]
            print(f"  Orbit {tid}: {len(lst)} images, {times[0].date()} ~ {times[-1].date()}")
        print("Suggestion: try --flight DES, or increase --days (e.g. 60~400), or expand --aoi.")
        return

    tid, pre, post = pairs[0]
    pre_name  = get_zipname(pre.properties)
    post_name = get_zipname(post.properties)
    print("Selected orbit =", tid)
    print("PRE :", pre_name,  pre.properties.get("startTime"))
    print("POST:", post_name, post.properties.get("startTime"))

    outdir = Path("S1_pick"); ensure_dir(outdir)

    def slim(x):
        p = x.properties
        return {
            "name": get_zipname(p),
            "startTime": p.get("startTime"),
            "url": p.get("url"),
            "relativeOrbit": p.get("relativeOrbit"),
            "pathNumber": p.get("pathNumber"),
            "flightDirection": p.get("flightDirection"),
        }
    with (outdir/"pairs.json").open("w") as f:
        json.dump([(tid, slim(pre), slim(post))], f, indent=2)

    with (outdir/"download.txt").open("w") as f:
        f.write(pre.properties.get("url","") + "\n")
        f.write(post.properties.get("url","") + "\n")
    print("Wrote S1_pick/download.txt and S1_pick/pairs.json")

    xml_txt = make_topsapp_xml(
        pre_name, post_name,
        dem_path=args.dem, pol=args.pol,
        swaths=swaths,
        az_looks=args.az_looks, rg_looks=args.rg_looks,
        filter_strength=args.filter_strength,
        do_unwrap=args.do_unwrap
    )
    Path("topsApp.xml").write_text(xml_txt)
    print("Wrote topsApp.xml (ready to run ISCE2 topsApp)")

    if args.plot or defcfg["plot"]:
        write_footprints_png(res, out_png="footprints.png")

    print("  # Step 1: Setup Earthdata Login credentials")
    print("  echo 'machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD>' > ~/.netrc")
    print("  chmod 600 ~/.netrc")

    print("  # Step 2: Download Sentinel-1 data")
    print("  aria2c -i S1_pick/download.txt -x 8 -s 8 -k 1M -d raw")
    print("  mkdir -p REF SEC orbits dem")
    print("  unzip raw/<PRE>.zip  -d REF/   &&  unzip raw/<POST>.zip -d SEC/")
    print("  # Place corresponding *.EOF into orbits/, DEM into", args.dem)
    print("  topsApp.py topsApp.xml --steps   # Verify steps, then run: topsApp.py topsApp.xml")

if __name__ == "__main__":
    main()
