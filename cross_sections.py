import numpy as np
import pandas as pd
from rich.pretty import pprint
from concreteproperties.material import Concrete, SteelBar
from concreteproperties.stress_strain_profile import (
    EurocodeNonLinear,
    RectangularStressBlock,
    SteelElasticPlastic,
)
from sectionproperties.pre.library.concrete_sections import concrete_circular_section
from concreteproperties.concrete_section import ConcreteSection
from matplotlib.figure import Figure

def createPileSection(section: pd.Series, materials: pd.DataFrame) -> ConcreteSection:
    '''
    Create Concrete Section Pile for the given concrete and rebar material
    '''
    concrete = Concrete(
        name=materials["name"]["Concrete"],
        density=2.4e-6,
        stress_strain_profile=EurocodeNonLinear(
            elastic_modulus=materials["E"]["Concrete"],
            ultimate_strain=0.0035,
            compressive_strength=materials["fck/fy"]["Concrete"],
            compressive_strain=0.0023,
            tensile_strength=3.8,
            tension_softening_stiffness=10e3,
        ),
        ultimate_stress_strain_profile=RectangularStressBlock(
            compressive_strength=materials["fck/fy"]["Concrete"]/materials["gamma"]["Concrete"],
            alpha=0.79,
            gamma=0.87,
            ultimate_strain=0.003,
        ),
        flexural_tensile_strength=3.8,
        colour="lightgrey",
    )

    steel = SteelBar(
        name=materials["name"]["Rebar"],
        density=7.85e-6,
        stress_strain_profile=SteelElasticPlastic(
            yield_strength=materials["fck/fy"]["Rebar"]/materials["gamma"]["Rebar"],
            elastic_modulus=materials["E"]["Rebar"],
            fracture_strain=0.05,
        ),
        colour="grey",
    )

    geom = concrete_circular_section(
        d=section["D"],
        n=64,
        dia=section["phi"],
        n_bar=int(section["nBars"]),
        n_circle=8,
        area_conc=np.pi * section["D"]**2 / 4,
        area_bar=np.pi * section["phi"]**2 / 4,
        cover=section["cover"],
        conc_mat=concrete,
        steel_mat=steel,
    )

    conc_sec = ConcreteSection(geom)
    # conc_sec.plot_section()

    return conc_sec

