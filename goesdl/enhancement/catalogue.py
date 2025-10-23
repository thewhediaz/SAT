import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from .scale import EnhancementScale

# 1. Colores de CONAGUA (12 colores)
CONAGUA_COLORS = [

    "#fcb0c8",  # -100
    "#fd5df0",  # -85
    "#e62f8c",  # -80
    "#dd0d0d",  # -75
    "#fe2530",  # -70
    "#a10a06",  # -65
    "#e95c14",  # -60
    "#e95c14",  # -55
    "#e95c14",  # -50
    "#fef756",  # -45
    "#fef756",  # -40
    "#fef756",  # -35
    "#0daa31",  # -30
    "#08721f",  # -25
    "#069bb6",  # -20
    "#1064da",  # -15
    "#01007f"  # -10

]

CONAGUA_BOUNDS = [-90,-85, -80, -75, -70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10]


# 2. Bordes de los rangos

# 3. Creamos el colormap discreto y la normalización
conagua_cmap = ListedColormap(CONAGUA_COLORS)
conagua_norm = BoundaryNorm(CONAGUA_BOUNDS, conagua_cmap.N)

# --- Función original ---
def _ircolor2() -> EnhancementScale:
    key_points = (-90, -5, 60)  # Se mantiene igual que antes
    color_maps = (conagua_cmap, "gray_r")  # Pasamos el objeto directamente
    cmap_name = f"IRCOLOR2 {int(key_points[0])}-{int(key_points[1])}K"
    return EnhancementScale.from_colormap(cmap_name, color_maps, key_points)
def _ircolor() -> EnhancementScale:
    # Define the properties of a stock palette
    key_points = (-73, -33, 60)  # (180, 240, 330)
    color_maps = ("jet_r", "gray_r")
    cmap_name = f"IRCOLOR {int(key_points[0])}-{int(key_points[1])}K"
    
    # Get a stock palette and configure the ticks for the colorbar
    return EnhancementScale.from_colormap(cmap_name, color_maps, key_points)

# Diccionario final
cmap: dict[str, EnhancementScale] = {
    "IRCOLOR": _ircolor(),
    "IRCOLOR2": _ircolor2()
}

