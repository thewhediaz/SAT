from typing import cast

from cartopy import crs as ccrs
from cartopy.crs import Globe, PlateCarree, Projection
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from matplotlib import pyplot as plt
from matplotlib.collections import QuadMesh
from matplotlib.ticker import MultipleLocator

from ..enhancement import EnhancementScale, cmap
from ..protocols import GeodeticRegion, SatImageData


class GSPlotParameter:

    title: str | tuple[str, str] = ""
    axis_label: tuple[str, str] = "", ""
    cbar_label: str = ""

    fig_width_px = 800
    fig_height_px = 800
    top_margin_px = 120
    bottom_margin_px = 110
    left_margin_px = 5
    right_margin_px = 5
    cbar_bottom_px = 80
    cbar_height_px = 16
    watermark_bottom_px = 7
    fig_dpi = 200
    

    def __init__(
        self,
        title: str | tuple[str, str] = "",
        axis_label: tuple[str, str] = ("", ""),
        cbar_label: str = "",
    ) -> None:
        self.title = title
        self.axis_label = axis_label
        self.cbar_label = cbar_label


class GSPlot:

    # The figure/image resolution (in number of dots per inch)
    fig_dpi: int = 200

    # The figure dimensions (in inches)
    fig_size: tuple[float, float] = 2.0, 2.0

    axes_box: tuple[float, float, float, float] = 0.0, 0.0, 0.0, 0.0
    cbar_box: tuple[float, float, float, float] = 0.0, 0.0, 0.0, 0.0
    plot_edges: tuple[float, float, float, float] = 0.0, 0.0, 0.0, 0.0
    watermark_loc: tuple[float, float] = 0.0, 0.0

    near_earth_scale: str = "10m"

    enhancement: EnhancementScale

    crs: Projection

    def __init__(self, enhancement: EnhancementScale | None = None) -> None:
        # Si se pasa enhancement lo usamos; si no, usamos IRCOLOR por defecto
        self.enhancement = enhancement or cmap["IRCOLOR"]

        target_globe = Globe(ellipse="WGS84")
        self.crs = PlateCarree(central_longitude=0.0, globe=target_globe)

    def plot(
        self,
        image: SatImageData,
        param: GSPlotParameter,
        save_path: str = "",
        show: bool = True,
        enhancement: EnhancementScale | None = None,  # NUEVO argumento
    ) -> None:
        # Si se pasa enhancement al plot, lo usamos; si no, usamos el actual
        if enhancement is not None:
            self.enhancement = enhancement
        # Setup the figure and axes boxes
        self._setup_boxes(param)

        # Create the figure and setup the axes
        fig = plt.figure("map", figsize=self.fig_size, dpi=self.fig_dpi)
        # 👉 Título principal fijo arriba de los otros dos
        # --- Título principal fijo arriba ---
        fig.subplots_adjust(top=0.95)  # deja espacio arriba para el título
        fig.suptitle("Sistema de Alerta Temprana Municipal de Acapulco SIATM-ACA", fontsize=9, y=0.98)


        # Create the axes ith the required projection, i.e., `target_crs` (see definition above)
        ax = fig.add_axes(self.axes_box, projection=self.crs)
        # --- Agregar geometrías del KML (sin alterar la proyección ni datos) ---
        import geopandas as gpd
        
        kml_path = "acapulco.kml"  # tu archivo KML
        gdf = gpd.read_file(kml_path, driver="KML")
        geoms = gdf.geometry.values
        
        # Agregar al GeoAxes de Cartopy
        ax.add_geometries(
            geoms,
            crs=ccrs.PlateCarree(),  # asegura que se proyecte correctamente
            facecolor='none',
            edgecolor='red',
            linewidth=0.5,
            zorder=2  # más alto que la malla de datos, pero menor que logos si quieres
        )

        # Adjust the subplots margins
        fig.subplots_adjust(
            left=self.plot_edges[0],
            bottom=self.plot_edges[1],
            right=self.plot_edges[2],
            top=self.plot_edges[3],
        )

        self._add_grid(ax, image.region, param.axis_label)

        self._add_admin_info(ax)

        #self._add_crosshair(ax)

        # Plot the data (in `gcrs`, see definition above) with a color map from the stock
        # shading;
        #   - center:  lon.shape == image.shape => "nearest" == "auto" == None | "gouraud"
        #   - corners: lon.shape != image.shape => "flat" == "auto" == None
        mesh = self._plot_data(ax, image)

        self._add_colorbar(mesh, fig, param.cbar_label)

        self._add_title(ax, param.title)

        
        from matplotlib.offsetbox import OffsetImage, AnnotationBbox
        import matplotlib.image as mpimg
        
        logo_img = mpimg.imread("LOGO.jpg")
        ab = AnnotationBbox(OffsetImage(logo_img, zoom=0.20), (0.90, 0.95), frameon=False, xycoords='axes fraction')
        ax.add_artist(ab)
        
        logo_img2 = mpimg.imread("SIATM.jpg")
        ab = AnnotationBbox(OffsetImage(logo_img2, zoom=0.08), (0.92, 0.85), frameon=False, xycoords='axes fraction')
        ax.add_artist(ab)




        # Add a watermark to the plot
        self._add_watermark(fig)

        # Save the media file (ensure the destination path does exist)
        if save_path:
            plt.savefig(save_path, dpi=self.fig_dpi, bbox_inches=None)

        # Show the plot
        if show:
            plt.show()
        else:
            plt.close()

    def save(
        self, save_path: str, image: SatImageData, param: GSPlotParameter
    ) -> None:
        self.plot(image, param, save_path, False)

    def show(self, image: SatImageData, param: GSPlotParameter) -> None:
        self.plot(image, param)

    def _add_admin_info(self, ax: plt.Axes) -> None:  # type: ignore
        # Create the Natural Earth projection (Plate-Carrée projection
        # on WGS84 ellipsoid)

        natearth_globe = ccrs.Globe(ellipse="WGS84")

        natearth_crs = ccrs.PlateCarree(
            central_longitude=0.0,
            globe=natearth_globe,
        )

        # Create political boundaries, in `natearth_crs` (see definition above)
        provinces = NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces",
            scale=self.near_earth_scale,
            facecolor="none",
            edgecolor='black',
            transform=natearth_crs,
        )

        countries = NaturalEarthFeature(
            category="cultural",
            name="admin_0_countries",
            scale=self.near_earth_scale,
            facecolor="none",
            edgecolor='black',
            transform=natearth_crs,
        )

        # Add political boundaries to the plot using our plot projection
        ax.add_feature(
            provinces, edgecolor="black", linewidth=0.4, transform=self.crs
        )
        ax.add_feature(
            countries, edgecolor="black", linewidth=0.4, transform=self.crs
        )

    def _add_colorbar(
        self, mesh: QuadMesh, fig: plt.Figure, label: str  # type: ignore
    ) -> None:

        # Add the colorbar
        caxes = fig.add_axes(self.cbar_box)
        cb = plt.colorbar(
            mesh,
            ticks=self.enhancement.cticks,
            orientation="horizontal",
            extend="both",
            cax=caxes,
        )

        # Create a minor tick locator for the colorbar
        minor_locator = MultipleLocator(5)

        # Set the colorbar tick characteristics
        cb.ax.xaxis.set_minor_locator(minor_locator)
        cb.ax.tick_params(
            labelsize=4.0,
            labelcolor="black",
            width=0.5,
            length=2.0,
            direction="out",
            pad=1.0,
        )
        cb.ax.tick_params(axis="both", which="minor", length=1.5, width=0.4)

        # Set the colorbar caption
        cb.set_label(label=label, size=5.0, color="black", weight="normal")

        # Set the colorbar characteristics
        cb.outline.set_linewidth(0.4)  # type: ignore[operator]

    def _add_crosshair(self, ax: plt.Axes) -> None:  # type: ignore
        # Add a crosshair to the plot
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        center_x = (xlims[1] + xlims[0]) / 2
        center_y = (ylims[1] + ylims[0]) / 2
        ax.axvline(x=center_x, color="red", linewidth=0.4)
        ax.axhline(y=center_y, color="red", linewidth=0.4)

    def _add_grid(
        self,
        ax: plt.Axes,  # type: ignore
        region: GeodeticRegion,
        labels: tuple[str, str],
    ) -> None:

        # Set gridline ticks characteristics
        ax.tick_params(
            left=True,
            right=True,
            bottom=True,
            top=True,
            labelleft=True,
            labelright=False,
            labelbottom=True,
            labeltop=False,
            length=0.0,
            width=0.05,
            labelsize=4.5,
            labelcolor="black",
        )

        # Add gridlines, set the longitude and latitude grid major tick locations,
        # set gridline characteristics, and configure gridline labels
        ax.gridlines(
            xlocs=region.xticks,
            ylocs=region.yticks,
            linewidth=0,
            linestyle="--",
            color="none",
            alpha=0.6,
            draw_labels=False,
        )

        # Set X-axis label characteristics
        ax.set_xticks(region.xticks, crs=self.crs)
        ax.xaxis.set_major_formatter(
            LongitudeFormatter(dateline_direction_label=True)
        )
        ax.set_xlabel(
            labels[0],
            color="black",
            fontsize=6.0,
            labelpad=3.0,
        )

        # Set Y-axis label characteristics
        ax.set_yticks(region.yticks, crs=self.crs)
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.set_ylabel(
            labels[1],
            color="black",
            fontsize=6.0,
            labelpad=3.0,
        )

        # Set the map limits, in `target_crs` (see definition above)
        ax.set_extent(region.extent, crs=self.crs)

    def _add_title(
        self,
        ax: plt.Axes,  # type: ignore
        title: str | tuple[str, str],
    ) -> None:
        # Set the title
        if isinstance(title, str):
            plt.title(title, fontsize=7.0)
            return

        ax.set_title(title[0], fontsize=6.5, loc="left")
        ax.set_title(title[1], fontsize=6.5, loc="right")

    def _add_watermark(self, fig: plt.Figure) -> None:  # type: ignore
        fig.text(
            self.watermark_loc[0],
            self.watermark_loc[1],
            "",
            horizontalalignment="right",
            verticalalignment="bottom",
            fontsize=4.0,
            color="gray",
            alpha=0.5,
            zorder=1000,
            transform=fig.transFigure,
        )

            
    def _plot_data(
        self,
        ax: plt.Axes,  # type: ignore
        data: SatImageData,
    ) -> QuadMesh:
        import numpy as np
        from scipy.interpolate import griddata
        from matplotlib.collections import QuadMesh
        from typing import cast
        from scipy.ndimage import gaussian_filter

    
        # Original coordinates y datos
        lon = data.grid.lon
        lat = data.grid.lat
        img = data.image
    
        # Crear un grid más fino para mejorar la resolución
        factor_resolucion = 5  # aumenta la resolución por este factor
        lon_fine = np.linspace(lon.min(), lon.max(), lon.shape[1] * factor_resolucion)
        lat_fine = np.linspace(lat.min(), lat.max(), lat.shape[0] * factor_resolucion)
        lon_fine_grid, lat_fine_grid = np.meshgrid(lon_fine, lat_fine)
    
        # Interpolación cúbica para suavizar la imagen
        img_fine = griddata(
            (lon.ravel(), lat.ravel()),
            img.ravel(),
            (lon_fine_grid, lat_fine_grid),
            method='cubic'  # también puedes usar 'linear' o 'nearest'
        )
        
        sigma_smooth = 1.0  # puedes ajustar entre 0.5 y 1.5
        img_smooth = gaussian_filter(img_fine, sigma=sigma_smooth)

    
        # Graficar con pcolormesh en el grid fino
        mesh = ax.pcolormesh(
            lon_fine_grid,
            lat_fine_grid,
            img_smooth,
            shading="auto",
            cmap=self.enhancement.cmap,
            norm=self.enhancement.cnorm,
            transform=data.grid.crs,
        )
    
        return cast(QuadMesh, mesh)


    
    def _setup_boxes(self, param: GSPlotParameter) -> None:
        self.fig_dpi = param.fig_dpi
        self.fig_size = (
            param.fig_width_px / self.fig_dpi,
            param.fig_height_px / self.fig_dpi,
        )
        self.axes_box = (
            param.left_margin_px / param.fig_width_px,
            param.bottom_margin_px / param.fig_height_px,
            (param.fig_width_px - param.left_margin_px - param.right_margin_px)
            / param.fig_width_px,
            (
                param.fig_height_px
                - param.top_margin_px
                - param.bottom_margin_px
            )
            / param.fig_height_px,
        )
        max_cbar_frac = 0.8  # máximo ancho como fracción del ancho total
        
        # ancho actual según márgenes
        cbar_width = (param.fig_width_px - param.left_margin_px - param.right_margin_px) / param.fig_width_px
        
        # limitar ancho al máximo
        cbar_width = min(cbar_width, max_cbar_frac)
        
        # centrar la barra
        cbar_x0 = (1 - cbar_width) / 2
        
        self.cbar_box = (
            cbar_x0,                                # x0 centrado
            param.cbar_bottom_px / param.fig_height_px,  # y0
            cbar_width,                             # width limitado
            param.cbar_height_px / param.fig_height_px   # height
        )


        self.plot_edges = (
            param.left_margin_px / param.fig_width_px,
            param.bottom_margin_px / param.fig_height_px,
            (param.fig_width_px - param.right_margin_px) / param.fig_width_px,
            (param.fig_height_px - param.top_margin_px) / param.fig_height_px,
        )
        self.watermark_loc = (
            (param.fig_width_px - param.watermark_bottom_px)
            / param.fig_width_px,
            param.watermark_bottom_px / param.fig_height_px,
        )
