# %%
import numpy as np
import xarray as xr
import pygmt
import math
import json
from cubbie.read_write_insar_utilities import isce_read_write

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
    #fig.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
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
    #fig.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
    fig.colorbar(frame=f"af+l\"{cfg['unw_colorbar_label']}\"")
    fig.savefig("unw.png", dpi=300)


##### LOS DATA #####
wavelength = 0.0555
los_def = z * (-wavelength / (4 * np.pi))
los_def.name = "los_deformation"  
los_def_cm = los_def * 100
los_def_cm.name = "los_deformation_cm"

vmax = np.nanmax(np.abs(los_def))
pygmt.makecpt(cmap="polar", series=cfg["los_series"], continuous=True, output="los.cpt")

# PLOT LOS
with pygmt.config(MAP_FRAME_TYPE="plain",
                  MAP_FRAME_PEN="1p,black",
                  FORMAT_GEO_MAP="ddd.xx",
                  MAP_DEGREE_SYMBOL="none"):
    fig = pygmt.Figure()
    fig.grdimage(grid=los_def_cm, region=[W, E, S, N],
                 projection=cfg["projection"], cmap="los.cpt",
                 frame=cfg["frame"])
    #fig.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
    fig.colorbar(frame=f"af+l\"LOS Displacement (cm)\"")
    fig.savefig("los.png", dpi=300)


##### COHERENCE DATA #####
xarr_coh, yarr_coh, coh_data = isce_read_write.read_scalar_data(cfg["coherence_file"], band=cfg["coherence_band"]) 

coh = xr.DataArray(
    np.asarray(coh_data, dtype=np.float32),
    coords={"lat": np.asarray(yarr_coh, dtype=np.float64),
            "lon": np.asarray(xarr_coh, dtype=np.float64)},
    dims=("lat", "lon"),
    name="coherence"
)

pygmt.makecpt(cmap="gray", series=[0, 1, 0.1], continuous=True, output="coh.cpt")

W_coh, E_coh = float(np.nanmin(coh.lon)), float(np.nanmax(coh.lon))
S_coh, N_coh = float(np.nanmin(coh.lat)), float(np.nanmax(coh.lat))

# PLOT COHERENCE
with pygmt.config(MAP_FRAME_TYPE="plain",
                  MAP_FRAME_PEN="1p,black",
                  FORMAT_GEO_MAP="ddd.xx",
                  MAP_DEGREE_SYMBOL="none"):
    fig_coh = pygmt.Figure()
    fig_coh.grdimage(grid=coh, region=[W_coh, E_coh, S_coh, N_coh],
                     projection=cfg["projection"], cmap="coh.cpt",
                     frame=cfg["frame"])
    fig_coh.colorbar(frame=f"af+l\"Coherence\"")
    fig_coh.savefig("coherence_map.png", dpi=300)


##### MASKED LOS #####
coherence_threshold = cfg["coherence_threshold"]

masked_los_def_cm = los_def_cm.copy(deep=True) 
masked_los_def_cm = masked_los_def_cm.where(coh >= coherence_threshold, np.nan)

with pygmt.config(MAP_FRAME_TYPE="plain",
                  MAP_FRAME_PEN="1p,black",
                  FORMAT_GEO_MAP="ddd.xx",
                  MAP_DEGREE_SYMBOL="none"):
    fig_masked_los = pygmt.Figure()
    fig_masked_los.grdimage(grid=masked_los_def_cm, region=[W, E, S, N],
                            projection=cfg["projection"], cmap="los.cpt", 
                            frame=cfg["frame"])
    #fig_masked_los.coast(shorelines=cfg["shorelines"], resolution=cfg["coast_resolution"])
    fig_masked_los.colorbar(frame=f"af+l\"LOS Displacement (cm)\"")
    fig_masked_los.savefig("los_co_masked.png", dpi=300)
