

# Set the name of the event
# Note: it will be used to deduce the data and media file storage paths
EVENT = "CURRENT_2H_colores_windy"

# Set the satellite (origin) identifier
# Note: it must be a string (e.g. "G16").
# Available origins:
#     "G16" to "G19": GOES 4th generation (R to U) series,
ORIGINS = "G19"

# Set the band identifier
# Note: Channel 13 - TOA Brightness Temperature - 10.3 μm for this example.
#     it can be a string (e.g. "C08") or a sequence of strings (e.g.
#     ["G07", "G15"]).
CHANNEL = "C13"

# Set the scene identifier
# Available scenes:
#     "F": "Full Disk",
#     "C": "CONUS (Contiguous United States)".
#    "M1": "Mesoscale", n = 1 or 2.
SCENE = "C"

# Set the central longitude and latitude for the region of interest (ROI)
LON_CEN = -99.91  # para Acapulco
LAT_CEN = 18   # para Acapulco
# Coordenadas de Acapulco
lon_cen, lat_cen = -99.91, 16.85

# Definir rango en grados para cubrir la región
lon_min, lon_max = lon_cen - 3, lon_cen + 3
lat_min, lat_max = lat_cen - 3, lat_cen + 3

# Set the extent of the ROI (longitude and latitude field of view in degrees)
ROI_EXTENT = (11, 11)


from pathlib import Path

# Set the path of the datasets repository
REPOSITORY_ROOT = Path("./data")
REPOSITORY_PATH = REPOSITORY_ROOT / EVENT
REPOSITORY_PATH.mkdir(parents=True, exist_ok=True)

# Set the path of the media repository
MEDIA_ROOT = Path("./media")
MEDIA_PATH = MEDIA_ROOT / EVENT
MEDIA_PATH.mkdir(parents=True, exist_ok=True)

# Set the media suffix (extension) for the generated images
# Note: this define the output file format.
MEDIA_SUFFIX = ".png"

# Set the media suffix for the generated animation
# Note: this define the output file format.
# Supported animation formats:
#     ".gif": Animated GIF, up to 256 colours palette,
#     ".png": Animated PNG, support true colours,
#     ".webp": Animated WebP, support true colours.
ANIMATION_SUFFIX = ".gif"

# Miscellaneous constants
DATE_FORMAT_DT = "%Y/%m/%d %H:%M:%S%:z"
DATE_FORMAT_IN = "%Y-%m-%dT%H:%M%z"
DATE_FORMAT_FG = "%Y/%m/%d %H:%M %Z"
ROI_CENTRE = (LON_CEN, LAT_CEN)

import datetime

# Get the current time in UTC (Coordinated Universal Time)
# The 'Z' (Zulu time) format is a military designation for UTC.
#
# Format the time in the required ISO 8601 format, truncating seconds and microseconds,
# and adding the literal 'Z' to indicate UTC.
# %Y: 4-digit year
# %m: Month as a number (01-12)
# %d: Day of the month as a number (01-31)
# T: Literal separator between date and time
# %H: Hour (24-hour clock) as a number (00-23)
# %M: Minute as a number (00-59)
# Z: Literal 'Z' character to indicate UTC

# 1. Get the current time in UTC
now_utc = datetime.datetime.now(datetime.timezone.utc)

# 2. Create a timedelta object to represent 2 hours 10 minutes
two_hours_delta = datetime.timedelta(minutes=10)

# 3. Subtract 2 hours from the current UTC time
past_time_utc = now_utc - two_hours_delta

# 4. Format the resulting times into the required format
# YYYY-MM-DDTHH:mmZ
TIME_START = past_time_utc.strftime(DATE_FORMAT_IN)
TIME_END = now_utc.strftime(DATE_FORMAT_IN)



from importlib import metadata, util



from cartopy import crs as ccrs

# Create the target (plotting) projection

target_globe = ccrs.Globe(ellipse="WGS84")

target_crs = ccrs.PlateCarree(
    central_longitude=0.0,
    globe=target_globe,
)


# Import the locator, datasource and downloader
from goesdl.goesr import GOESProductLocatorCMIP
from goesdl.datasource import DatasourceAWS
from goesdl.downloader import Downloader

# Initialize the product locator for GridSat-GOES
#
# Available scenes:
#     "F": "Full Disk",
#     "C": "CONUS (Contiguous United States)".
#
# Available origins:
#     "G08" to "G12": GOES 2nd generation (I to M) series,
#     "G13" to "G15": GOES 3rd generation (N to P) series.
#
# Note: The `origin` parameter can take a list of origins
#       as argument, e.g. origins=["G11", "G12"].
locator = GOESProductLocatorCMIP(scene=SCENE, channels=CHANNEL, origin=ORIGINS)

# GOES-R database is available from NOAA's AWS archive
datasource = DatasourceAWS(locator)

# Initialize the downloader with the locator and datasource
downloader = Downloader(
    datasource=datasource,
    locator=locator,
    repository=REPOSITORY_PATH,
    date_format=DATE_FORMAT_IN,
)

# --- Optional: Show both times for comparison ---
print(f"Coverage UTC start time : {TIME_START}")
print(f"Coverage UTC end time   : {TIME_END or TIME_START}\n")

# List the datasets within a given date range
#
# Note: To download for a specific date and time you might do
#       `downloader.download_files(start="2012-08-23T00:00Z")`.
files = downloader.list_files(
    start=TIME_START,
    end=TIME_END,
)

# Filter the resulting list, get the last twelve files, at most
files = files[-1:]

# Download the datasets within the filtered range
downloader.get_files(file_paths=files)

if not files:
    print("Unable to acquire files: no file found in the specified range\n")
else:
    print(f"Successfully acquired {len(files)} files")

# Note: The returned list of file paths are relative to the
#     repository root, to match the absolute path in the
#     datasource file syste, so, we need to reconstruct
#     the dataset paths to access our local copy.
dataset_paths = [REPOSITORY_PATH / file for file in files]


from goesdl.geodesy import RectangularRegion

# Calculate the data domain and the plot limits
region = RectangularRegion(domain=((lon_min, lon_max), (lat_min, lat_max)))


import numpy as np
from netCDF4 import Dataset
from goesdl.goesr import GOESLatLonGrid, GOESImage, GOESDatasetInfo

data_seq = []
metadata_seq = []

if isinstance(ORIGINS, str):
    ORIGINS = [ORIGINS]

# Set the path and file name of the dataset and extract required data
for i, dataset_path in enumerate(dataset_paths):
    with Dataset(dataset_path, "r") as dataframe:

        # Get all required data by entity in separate calls: lat. and
        # lon. grids, image data, coverage time*, metadata*).
        # *) These data are not required to plot the image, they are
        #    jus used to extract information for labeling the image.
        grid = GOESLatLonGrid(dataframe, region)  #, corners=True

        data = GOESImage(dataframe, grid)
        # Convertir de Kelvin a Celsius in-place
        data.raster.data[:] = data.raster.data - 273.15
        
        data_seq.append(data)


        metadata = GOESDatasetInfo(dataframe, CHANNEL)
        metadata_seq.append(metadata)

    print(f"\nDataset {i + 1}/{len(dataset_paths)}")
    print("-----------------------------------------------------")
    print(f"Record:    {metadata.dataset_name}")
    print("-----------------------------------------------------")
    print(f"Database                  : {metadata.database_name}")
    print(f"Time coverage start       : {metadata.coverage_start:{DATE_FORMAT_DT}}")
    print(f"Time coverage end         : {metadata.coverage_end:{DATE_FORMAT_DT}}")
    print(f"Platform                  : {metadata.plaform_name}")
    print(f"Measurement               : {metadata.measurement_name}")
    print(f"Band                      : Band {metadata.band_id} — {metadata.band_wavelength:.1f} μm")
    print(f"Minimum value             : {np.nanmin(data.image):.2f} [{metadata.units}]")
    print(f"Maximum value             : {np.nanmax(data.image):.2f} [{metadata.units}]")
    print(f"ROI shape (rows×columns)  : {data.image.shape[0]}×{data.image.shape[1]} px")
    print("-----------------------------------------------------\n")



from goesdl.enhancement import EnhancementScale, show_colormap, cmap

# Get a stock palette and configure the ticks for the colorbar, for this example
enhancement = cmap["IRCOLOR"]

# Create a preview of the color map
measurement = f"{metadata.measurement_name} [{metadata.units}]"

#show_colormap(enhancement, measurement)  #, offset=-269.4


from goesdl.plotting import GSPlot, GSPlotParameter

# Create the basic plotter
plotter = GSPlot()

# Set the plotting projection
plotter.crs = target_crs

# Create list of generated image files (optional)
# Note: this will be used to generate the animation file
media_files = []

# Generate the image files (or visualise)
for i, dataset_path in enumerate(dataset_paths):
    # Get the relevant plotting data
    data = data_seq[i]
    metadata = metadata_seq[i]


    # Set the title
    title_left = "Elaboro: C-31 de la CGPCyB de Acapulco"
    
    
    
    from datetime import timezone, timedelta
    
    # Convertir a UTC-6
    dt_local = metadata.coverage_end.astimezone(timezone(timedelta(hours=-6)))
    
    # Crear string formateado en español (día, mes, hora)
    dias = ["lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo"]
    meses = ["enero", "febrero", "marzo", "abril", "mayo", "junio",
             "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre"]
    
    dia_semana = dias[dt_local.weekday()]
    dia = dt_local.day
    mes = meses[dt_local.month - 1]
    hora = dt_local.strftime("%I:%M %p").lower()  # ejemplo: "10:19 am"
    
    title_right = f"{dia_semana}, {dia} de {mes}, {hora}"

    # Set XY-axis label
    x_label = ""
    y_label = ""

    # Set the colorbar caption
    # f"{cmi_entity} {cmi_measure} [{cmi_units}]" -> TOA Brightness Temperature
    # cb_label = f"{image.metadata.long_name} @ {metadata.wavelength:.1f} μm [{image.metadata.units}]"
    cb_label = "GOES-19 Canal 13 (IR) Temperatura de brillo (°C)"

    # Create the plot parameters
    param = GSPlotParameter((title_left, title_right), (x_label, y_label), cb_label)
    param.fig_dpi = 300  # Aumenta la resolución
    param.fig_width_px = 1400  # Aumenta la resolución
    param.fig_height_px = 1400  # Aumenta la resolución


    # Save the media file (ensure the destination path does exist)
    file_path = MEDIA_PATH/ dataset_path.with_suffix(MEDIA_SUFFIX).name
    

    # Create the figure and setup the axes
    if not file_path.exists():
        plotter.plot(data, param, file_path,show= False,enhancement=cmap["IRCOLOR2"])

    # Add the generate image file path to the media file list
    media_files.append(file_path)


import imageio.v2 as imageio  # ← usar la versión 2 explícitamente
from pathlib import Path

MEDIA_PATH = Path("./media/"+EVENT)

from pathlib import Path

# Carpeta de medios
MEDIA_PATH = Path("./media/" + EVENT)

# Obtener todos los PNG ordenados por fecha de modificación
all_png = sorted(MEDIA_PATH.glob("*.png"), key=lambda f: f.stat().st_mtime)

# Conservar solo los últimos 18
last_18_png = all_png[-38:]

# Archivos a eliminar: todos los PNG que no estén en last_18_png + todos los GIF
for file in MEDIA_PATH.glob("*"):
    if file.suffix in [".png"] and file not in last_18_png:
        file.unlink()
        print(f"Eliminado: {file.name}")



