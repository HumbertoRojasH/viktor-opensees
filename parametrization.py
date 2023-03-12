from viktor.parametrization import NumberField
from viktor.parametrization import OptionField
from viktor.parametrization import Text

from viktor.parametrization import ViktorParametrization
from viktor.parametrization import Step
from viktor.parametrization import Section

class Parametrization(ViktorParametrization):

    s1 = Step('Modelo Paramétrico',views=['data_view','about'])
    s1.t1 = Text('## Geometría del edificio')
    s1.nx = NumberField('Nx', default=3, flex=50)
    s1.ex = NumberField('Lx $(m)$', default=4, flex=50)
    s1.ny = NumberField('Ny', default=5, flex=50)
    s1.ey = NumberField('Ly $(m)$', default=4, flex=50)
    s1.nz = NumberField('Nz', default=8, flex=50)
    s1.ez = NumberField('Lz $(m)$', default=3, flex=50)
    s1.t2 = Text('## Elementos estructurales')
    s1.cb = NumberField('Cb $(m)$', default=0.6, flex=50)
    s1.ch = NumberField('Ch $(m)$', default=0.6, flex=50)
    s1.vb = NumberField('Vb $(m)$', default=0.3, flex=50)
    s1.vh = NumberField('Vh $(m)$', default=0.5, flex=50)
    s1.le = NumberField('Lh $(m)$', default=0.2, flex=50)
    s1.t3 = Text('## Propiedades del material')
    s1.fc = NumberField("f'c $(kg/cm^2)$", default=210, flex=100)
    s1.d = NumberField('ρ $(kg/m^3)$', default=2400, flex=100)

    s2 = Step('Análisis Lineal Estático',views=['data_static','three_d_model','data_static2','data_static3'])
    s2.sec1 = Section("Dirección de Análisis")
    s2.sec1.dir = OptionField('', options=["Sismo X-X","Sismo Y-Y"], variant='radio-inline', flex=100, default="Sismo X-X")
    s2.sec1.sc = NumberField('Factor de escala', default=200, flex=60)

    s2.sec2 = Section('Cargas de Gravedad')
    s2.sec2.t2 = Text('Peso Sísmico: 100%CM + 25%CV')
    s2.sec2.cv = NumberField('Carga viva $(kg/m^2)$', default=250, flex=100)
    s2.sec2.cl = NumberField('Carga por losa $(kg/m^2)$', default=300, flex=100)
    s2.sec2.ca = NumberField('Carga por acabados $(kg/m^2)$', default=100, flex=100)
    s2.sec2.ct = NumberField('Carga por tabiqueria $(kg/m^2)$', default=100, flex=100)

    s2.sec3 = Section('Parámetros de Sitio')
    s2.sec3.z = NumberField("Z", default=0.45, flex=50)
    s2.sec3.u = NumberField("U", default=1, flex=50)
    s2.sec3.s = NumberField("S", default=1, flex=50)
    s2.sec3.tp = NumberField("Tp", default=0.4, flex=50)
    s2.sec3.tl = NumberField("Tl", default=2.5, flex=50)
    s2.sec3.r = NumberField("R", default=8, flex=50)

    s3 = Step('Análisis Modal Espectral',views=['data_dynamic','data_dynamic2','data_dynamic3'])