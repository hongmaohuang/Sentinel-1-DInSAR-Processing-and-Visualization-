import numpy as np
import xarray as xr
import pygmt
import math
import json
from cubbie.read_write_insar_utilities import isce_read_write

# 讀取設定檔
with open("s1_config.json", "r") as f:
    cfg = json.load(f)

pygmt.config(COLOR_NAN="white")

##### PHASE DATA #####
xarr, yarr, data = isce_read_write.read_scalar_data(cfg["phase_file"], band=cfg["phase_band"])
phase = np.arctan2(np.imag(data), np.real(data))
data = phase

z = xr.DataArray(
    np.asarray(data, dtype=np.float32),
    coords={"lat": np.asarray(yarr, dtype=np.float64),
            "lon": np.asarray(xarr, dtype=np.float64)},
    dims=("lat", "lon"),
    name="intf"
)

# PLOT PHASE
pygmt.makecpt(cmap=cfg["phase_cmap"], series=cfg["phase_series"], continuous=True, output="phase.cpt")
W, E = float(np.nanmin(z.lon)), float(np.nanmax(z.lon))
S, N = float(np.nanmin(z.lat)), float(np.nanmax(z.lat))
with pygmt.config(MAP_FRAME_TYPE="plain",
                  MAP_FRAME_PEN="1p,black",
                  FORMAT_GEO_MAP="ddd.xx",
                  MAP_DEGREE_SYMBOL="none"):
    fig = pygmt.Figure()
    fig.grdimage(grid=z, region=[W, E, S, N],
                 projection=cfg["projection"], cmap="phase.cpt",
                 frame=cfg["frame"])
    fig.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
    fig.colorbar(frame=f"af+l\"{cfg['phase_colorbar_label']}\"")
    fig.savefig("phase.png", dpi=300)

##### UNW DATA #####
xarr, yarr, data = isce_read_write.read_scalar_data(cfg["unw_file"], band=cfg["unw_band"])

z = xr.DataArray(
    np.asarray(data, dtype=np.float32),
    coords={"lat": np.asarray(yarr, dtype=np.float64),
            "lon": np.asarray(xarr, dtype=np.float64)},
    dims=("lat", "lon"),
    name="intf"
)

# PLOT UNW
pygmt.makecpt(cmap=cfg["unw_cmap"], series=cfg["unw_series"], continuous=True, output="unw.cpt")
W, E = float(np.nanmin(z.lon)), float(np.nanmax(z.lon))
S, N = float(np.nanmin(z.lat)), float(np.nanmax(z.lat))
with pygmt.config(MAP_FRAME_TYPE="plain",
                  MAP_FRAME_PEN="1p,black",
                  FORMAT_GEO_MAP="ddd.xx",
                  MAP_DEGREE_SYMBOL="none"):
    fig = pygmt.Figure()
    fig.grdimage(grid=z, region=[W, E, S, N],
                 projection=cfg["projection"], cmap="unw.cpt",
                 frame=cfg["frame"])
    fig.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
    fig.colorbar(frame=f"af+l\"{cfg['unw_colorbar_label']}\"")
    fig.savefig("unw.png", dpi=300)
