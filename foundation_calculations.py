import numpy as np
from dataclasses import dataclass
from handcalcs.decorator import handcalc
from PyNite import FEModel3D
import pandas as pd
import math
import streamlit as st

@dataclass
class Material:
    '''
    Define general material
    '''
    name: str
    E: float
    G: float
    nu: float
    specificWeight: float

@dataclass
class Concrete(Material):
    '''
    Define concrete material
    '''
    fc: float

class Rebar(Material):
    '''
    Define reinforcing steel material
    '''
    fy: float

@dataclass
class Soil:
    '''
    Define soil object
    density: specific weight of soil
    phi: Internal friction angle of soil
    c: Cohesion of soil
    '''
    density: float
    phi: float
    c: float
    k: float
    thickness: float

    def CalculateSubgradeModulus(self,k):
        self.k = k
        return self.k


@dataclass
class Reaction():
    '''
    N: Axial force in kN
    Mx: Longitudinal bending moment in kNm
    My: Transversal bending moment in kNm
    T: Torsor moment in kNm
    Vx: Transversal shear force in kN
    Vy: Longitudinal shear force in kN
    '''
    combination: str = ""
    N: float=0.0
    Mx: float= 0.0
    My: float= 0.0
    T: float= 0.0
    Vx: float=0.0
    Vy: float=0.0
    isCalculated: bool=False
    typeOfReaction: str=""

@dataclass
class Pile:
    '''
    Definition of a single pile
    '''
    name: str
    coordinates: tuple
    diameter: float
    length: float
    pierReactions: dict[Reaction] = None
    area: float = 0.0
    pileReactions: dict[Reaction] = None
    concrete: Concrete = Concrete("C25",30000,10000,0.2,25,25)
    rebars: list[Rebar] = None

    def CalculatePileArea(self):
        self.area = np.pi*self.diameter**2/4



def CalculatePiles(Nx: float, Ny: float, Sx: float, Sy: float, diameter: float, length: float) -> dict[Pile]:
    zx, zy = (Nx-1)*Sx, (Ny-1)*Sy
    pileCounter = 1
    output = {}
    for i in range(Nx):
        for j in range(Ny):
            x = zx/2 - Sx*i
            y = zy/2 - Sy*j
            pile = Pile(
                name = f"P{pileCounter}",
                coordinates = (x,y),
                diameter = diameter,
                length = length
            )
            output[pile.name] = pile
            pileCounter += 1
    return output

def CalculateSinglePileReactions(piles: dict[Pile], reactions: dict[Reaction]):
    '''
    Calculate the single reaction of a pile given the reaction of the pier and piles geometry
    '''
    areaTotal, denominator2x, denominator2y, denominator3 = 0, 0, 0, 0
    for pile in piles.values():
        pile.CalculatePileArea()
        areaTotal += pile.area
        denominator2x += pile.area*pile.coordinates[0]**2
        denominator2y += pile.area*pile.coordinates[1]**2
        denominator3 += pile.area*(pile.coordinates[0]**2 + pile.coordinates[1]**2)
    
    for pile in piles.values():
        A = pile.area
        xi, yi = pile.coordinates[0], pile.coordinates[1]
        
        combosPerPile = {}
        for comboName, reaction in reactions.items():
            # Pier reaction per single combination
            N = reaction.N
            Mx, My, T = reaction.Mx, reaction.My, reaction.T
            Vx, Vy = reaction.Vx, reaction.Vy
            # Reactions of single pile
            Ni = A/areaTotal*N - A*yi/denominator2y*Mx + A*xi/denominator2x*My
            Vxi = A/areaTotal*Vx - yi*A**2/denominator3*T
            Vyi = A/areaTotal*Vy + xi*A**2/denominator3*T
            combosPerPile[comboName] = Reaction(combination=comboName, N = Ni, Vx = Vxi, Vy = Vyi)
    
        pile.pileReactions = combosPerPile

def CalculatePileFEModel(pile: Pile, soil: pd.DataFrame, stepSprings: float)->FEModel3D:
    beam = FEModel3D()
    num_nodes = round(pile.length/stepSprings)+1
    spacing = stepSprings
    # Create column "zStart" to make easier the calculus of Kh at each soil layer
    soil["zStart"] = 0
    for idx, row in soil.iterrows():
        if idx > 0: soil["zStart"][idx] = soil["zEnd"][idx-1] 

    # Generate nodes
    for i in range(num_nodes):
        y = -i*spacing
        beam.add_node("N"+str(i+1),X=0,Y=y,Z=0)
        # Add supports to the nodes
        if i==0: # top of the pile 
            beam.def_support("N"+str(i+1),support_DX=False,support_DY=False,support_DZ=True,support_RX=False,support_RY=True,support_RZ=True)
        elif i==num_nodes-1: # bottom of the pile
            beam.def_support("N"+str(i+1),support_DX=True,support_DY=True,support_DZ=True,support_RX=False,support_RY=False,support_RZ=False)
        else:
            pd.to_numeric(soil["kh"],errors="coerce")
            maskZ = (soil["zEnd"] <= y) & (soil["zStart"] > y)
            k = soil["kh"][maskZ].iloc[0]
            kCalc = k*pile.diameter*stepSprings
            beam.def_support_spring("N"+str(i+1),"DX",kCalc)

    # Create pile concrete material, member object and mechanical properties
    beam.add_material(pile.concrete.name,1000*pile.concrete.E,1000*pile.concrete.G,pile.concrete.nu,pile.concrete.specificWeight)
    I = math.pi*(pile.diameter/2)**4/4
    J = math.pi*(pile.diameter/2)**4/2
    beam.add_member(pile.name,"N1","N"+str(num_nodes),pile.concrete.name,I,I,J,pile.area)

    # Generate load combinations
    firstNodeTag = list(beam.Nodes.keys())[0] # Top pile node (first node created in the FEM model)
    counter = 1
    for combo, reaction in pile.pileReactions.items():
        if reaction.isCalculated == True:
            beam.add_load_combo(combo,{f"N_{counter}": 1.0, f"V_{counter}": 1.0})
            beam.add_node_load(firstNodeTag,"FY",-reaction.N,case=f"N_{counter}") # Positive pile effort from "pileReactions" object means compression
            V = (reaction.Vx**2 + reaction.Vy**2)**0.5
            beam.add_node_load(firstNodeTag,"FX",V,case=f"V_{counter}")
            counter +=1

            # print(f"Combo: {combo}, N = {reaction.N} kN, V = {V} kN")

    # Analyze beam
    beam.analyze()

    return beam

def ReadExcelPierReactions(fileName: str, directions:dict)->pd.DataFrame:
    '''
    Read Excel file with pier reactions are return dataframe filtered with relevant information
    '''
    df = pd.read_excel(fileName,header=1,skiprows=[2])
    effortsDirections = list(directions.values())
    return df












