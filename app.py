from viktor.core import ViktorController
from viktor.views import PlotlyView, PlotlyResult, GeometryResult, GeometryView, PDFView, PDFResult
from viktor.geometry import Polygon, RectangularExtrusion, BidirectionalPattern,LinearPattern
from parametrization import Parametrization
from pathlib import Path
from data import *
from plotly.subplots import make_subplots
import plotly.graph_objects as go


class Controller(ViktorController):
    label = "Viktor-OpenseesPy"
    parametrization = Parametrization(width=20)
    viktor_enforce_field_constraints = True

    @GeometryView("Modelo 3D", duration_guess=10, default_shadow=True, update_label="Generar Gráfico")
    def data_view(self, params, **kwargs):
        s = params.s1
        CONCRETE_SLAB = Material("Concrete", color=Color(220, 220, 220), threejs_roughness=1, threejs_opacity=1)
        CONCRETE_BEAM = Material("Concrete", color=Color(128, 128, 118), threejs_roughness=1, threejs_opacity=1)

        lx = s.ex * s.nx
        ly = s.ey * s.ny
        
        c_start = Point(0,0,0)
        c_end = Point(0,0,s.ez*s.nz)
        c_line = Line(c_start, c_end)
        col = RectangularExtrusion(s.cb, s.ch, c_line, material=CONCRETE_BEAM)

        cols = BidirectionalPattern(col,(1,0,0),(0,1,0), s.nx+1, s.ny+1, s.ex, s.ey)

        vx_start = Point(0,0,s.ez - s.vh/2)
        vx_end = Point(s.ex*s.nx,0,s.ez - s.vh/2)
        vx_line = Line(vx_start, vx_end)
        viga_x = RectangularExtrusion(s.vb, s.vh, vx_line, material=CONCRETE_BEAM, profile_rotation=90)
        vigasx = BidirectionalPattern(viga_x,(0,1,0),(0,0,1), s.ny+1, s.nz, s.ey, s.ez)

        vy_start = Point(0,0,s.ez - s.vh/2)
        vy_end = Point(0,s.ey*s.ny,s.ez - s.vh/2)
        vy_line = Line(vy_start, vy_end)
        viga_y = RectangularExtrusion(s.vb, s.vh, vy_line, material=CONCRETE_BEAM)
        vigasy = BidirectionalPattern(viga_y,(1,0,0),(0,0,1), s.nx+1, s.nz, s.ex, s.ez)

        TOLERANCE = 0.0005

        l_start = Point(s.nx * s.ex/2, s.ny * s.ey/2, s.ez+ TOLERANCE)
        l_end = Point(s.nx * s.ex/2, s.ny * s.ey/2, s.ez- s.le)
        l_ext = Line(l_start, l_end)
        slabz = RectangularExtrusion(
                width=lx + s.cb + TOLERANCE,
                height=ly + s.ch + TOLERANCE,
                line=l_ext,
                material=CONCRETE_SLAB)

        slabs = LinearPattern(slabz,(0,0,1), s.nz, s.ez)

        SOIL = Material("Soil", color=Color(171, 139, 98))
        p_1 = Point(-lx, -ly, 0)
        p_2 = Point(2 * lx, -ly, 0)
        p_3 = Point(2 * lx, 2 * ly, 0)
        p_4 = Point(-lx, 2 * ly, 0)
        suelo = [Polygon([p_1, p_2, p_3, p_4], material=SOIL)]

        return GeometryResult([cols,vigasx,vigasy,slabs]+suelo)

    @PDFView("Acerda de", duration_guess=1)
    def about(self, params, **kwargs):
        file_path = Path(__file__).parent / 'parametros.pdf'
        return PDFResult.from_path(file_path)

    @GeometryView("Modelo 3D", duration_guess=10)
    def three_d_model(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, u = analysis_static(Nodes, Elems, Diap, params)
        nodes = Nodes[:, 1:4]
        elems = Elems[:, 1:3]
        deform_scale = displacement(nodes, u, scale=params.s2.sec1.sc)
        group = plot_structure_geometry_view(nodes, elems, deform_scale=deform_scale)
        return GeometryResult(group)

    @PlotlyView("Derivas", duration_guess=10, update_label="Generar Resultados")
    def data_static2(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, u = analysis_static(Nodes, Elems, Diap, params)
        fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.05, subplot_titles=("Drift X (‰)", "Drift Y (‰)"))

        fig.add_scatter(x=[0] + df["DriftX(‰)"].to_list(), y=[0] + df["H(m)"].to_list(),
                        line=dict(color='black', width=1), row=1, col=1)
        fig.add_scatter(x=[7] * (params.s1.nz + 1), y=[0] + df["H(m)"].to_list(),
                        line=dict(color='red', width=3, dash='dash'), mode='lines', row=1, col=1)
        fig.add_scatter(x=[0] + df["DriftY(‰)"].to_list(), y=[0] + df["H(m)"].to_list(),
                        line=dict(color='black', width=1), row=1, col=2)
        fig.add_scatter(x=[7] * (params.s1.nz + 1), y=[0] + df["H(m)"].to_list(),
                        line=dict(color='red', width=3, dash='dash'), mode='lines', row=1, col=2)

        fig.update_layout(
            autosize=False,
            showlegend=False,
            margin=dict(l=300, r=300, t=50, b=50), )
        fig.update_layout(
            yaxis=dict(
                title="Altura (m)",
                tickmode='linear',
                tick0=0,
                dtick=s.ez),
            yaxis2=dict(
                tickmode='linear',
                tick0=0,
                dtick=s.ez),
        )

        return PlotlyResult(fig.to_json())

    @PlotlyView("Resultados", duration_guess=10, update_label="Generar Resultados")
    def data_static3(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, u = analysis_static(Nodes, Elems, Diap, params)

        fig = go.Figure(data=[go.Table(
            columnwidth=[50, 50, 50, 50, 50, 50],
            header=dict(values=list(df.columns),
                        line_color='darkslategray',
                        fill_color='royalblue',
                        align='center',
                        font=dict(color='white', size=12), ),
            cells=dict(values=[df['Nivel'], df['H(m)'], df['V(kN)'], df['Ux(cm)'], df['Uy(cm)'], df['0.75R*Ux(cm)'],
                               df['0.75R*Uy(cm)'], df['DriftX(‰)'], df['DriftY(‰)']],
                       fill_color='white',
                       align='center',
                       line_color='darkslategray', ))
        ])

        fig.update_layout(title_text='Tabla de Resultados del Análisis Lineal Estático', )

        return PlotlyResult(fig.to_json())
    
    @PDFView("Acerda de", duration_guess=1)
    def data_static(self, params, **kwargs):
        file_path = Path(__file__).parent / 'estatic_modal.pdf'
        return PDFResult.from_path(file_path)

    @PlotlyView("Modos", duration_guess=10, update_label="Generar Resultados")
    def data_dynamic(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, df2 = analysis_modal(Nodes, Elems, Diap, params)

        fig = go.Figure(data=[go.Table(
            columnwidth=[50, 50, 50, 50, 50, 50],
            header=dict(values=list(df.columns),
                        line_color='darkslategray',
                        fill_color='royalblue',
                        align='center',
                        font=dict(color='white', size=12), ),
            cells=dict(values=[df['Modo'], df['T(s)'], df['SumUx'], df['SumUy'], df['SumRz']],
                       fill_color='white',
                       align='center',
                       line_color='darkslategray', ))
        ])
        fig.update_layout(title_text='Tabla de Modos de Vibración', )

        return PlotlyResult(fig.to_json())

    @PlotlyView("Derivas", duration_guess=10, update_label="Generar Resultados")
    def data_dynamic2(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, df2 = analysis_modal(Nodes, Elems, Diap, params)

        fig = go.Figure()
        fig.add_scatter(x=[0] + df2["Δx(‰)"].to_list(), y=[0] + df2["H(m)"].to_list(), line=dict(color='red', width=1),
                        name="Drift X (‰)")
        fig.add_scatter(x=[0] + df2["Δy(‰)"].to_list(), y=[0] + df2["H(m)"].to_list(), line=dict(color='blue', width=1),
                        name="Drift Y (‰)")
        fig.add_scatter(x=[7] * (params.s1.nz + 1), y=[0] + df2["H(m)"].to_list(),
                        line=dict(color='black', width=3, dash='dash'), mode='lines', name="Limite E030(‰)")

        fig.update_layout(
            autosize=False,
            margin=dict(l=300, r=300, t=50, b=50), )
        fig.update_layout(xaxis=dict(title='Drift (‰)'),
                          yaxis=dict(
                              title='Altura (m)',
                              tickmode='linear',
                              tick0=0,
                              dtick=s.ez),
                          )

        return PlotlyResult(fig.to_json())

    @PlotlyView("Resultados", duration_guess=10, update_label="Generar Resultados")
    def data_dynamic3(self, params, **kwargs):
        s = params.s1
        Nodes, Elems, Diap = GeoModel(s.ex, s.ey, s.ez, s.nx, s.ny, s.nz)
        df, df2 = analysis_modal(Nodes, Elems, Diap, params)

        fig = go.Figure(data=[go.Table(
            columnwidth=[50, 50, 50, 50, 50, 50],
            header=dict(values=list(df2.columns),
                        line_color='darkslategray',
                        fill_color='royalblue',
                        align='center',
                        font=dict(color='white', size=12), ),
            cells=dict(values=[df2['Nivel'], df2['H(m)'], df2['Vx(kN)'], df2['Vy(kN)'], df2['Ux(cm)'], df2['Uy(cm)'],
                               df2['Δx(‰)'], df2['Δy(‰)']],
                       fill_color='white',
                       align='center',
                       line_color='darkslategray', ))
        ])

        fig.update_layout(title_text='Tabla de Resultados del Análisis Modal Espectral', )

        return PlotlyResult(fig.to_json())