import numpy as np
from dataclasses import dataclass
from handcalcs.decorator import handcalc
from PyNite import FEModel3D

@dataclass
class Material:
    '''
    Define general material
    '''
    name: str
    E: float
    G: float
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
    N: float
    Mx: float= 0.0
    My: float= 0.0
    T: float= 0.0
    Vx: float=0.0
    Vy: float=0.0

@dataclass
class Pile:
    '''
    Definition of a single pile
    '''
    name: str
    coordinates: tuple
    diameter: float
    length: float
    pierReaction: Reaction
    area: float = 0.0
    pileReaction: Reaction = None
    concrete: Concrete = Concrete("C25",30000,10000,25,25)
    rebars: list[Rebar] = None
    cover: float = 0.05

    def CalculatePileArea(self):
        self.area = np.pi*self.diameter**2/4



def CalculatePiles(Nx: float, Ny: float, Sx: float, Sy: float, diameter: float, length: float, pierReactions: Reaction) -> list[Pile]:
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
                length = length,
                pierReaction = pierReactions
            )
            output[pile.name] = pile
            pileCounter += 1
    return output

def CalculateSinglePileReactions(piles: dict):
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
        N = pile.pierReaction.N
        Mx, My, T = pile.pierReaction.Mx, pile.pierReaction.My, pile.pierReaction.T
        Vx, Vy = pile.pierReaction.Vx, pile.pierReaction.Vy

        # Reactions of single pile
        Ni = A/areaTotal*N - A*yi/denominator2y*Mx + A*xi/denominator2x*My
        Vxi = A/areaTotal*Vx - yi*A**2/denominator3*T
        Vyi = A/areaTotal*Vy + xi*A**2/denominator3*T
        pile.pileReaction = Reaction(N = Ni, Vx = Vxi, Vy = Vyi)

def CalculatePileFEModel(pile: Pile, k: float, stepSprings: float)->FEModel3D:
    beam = FEModel3D()
    num_nodes = round(pile.length/stepSprings)+1
    spacing = pile.length/round(pile.length/stepSprings)
    kCalc = k*pile.diameter*stepSprings

    # Generate nodes
    for i in range(num_nodes):
        beam.add_node("N"+str(i+1),0,-i*spacing,0)
        # Add supports to the nodes
        if i==0: # top of the pile 
            beam.def_support("N"+str(i+1),False,False,True,False,True,True)
        elif i==num_nodes-1: # bottom of the pile
            beam.def_support("N"+str(i+1),True,True,True,False, False, False)
        else:
            beam.def_support_spring("N"+str(i+1),"DX",kCalc)
 
    beam.add_material(pile.concrete.name,pile.concrete.E,pile.concrete.G,0.3,pile.concrete.specificWeight)
    I = (pile.diameter/2)**4/4
    J = 1
    beam.add_member(pile.name,"N1","N"+str(num_nodes),pile.concrete.name,I,I,J,pile.area)
    beam.add_node_load("N1","FY",-pile.pileReaction.N)
    V = (pile.pileReaction.Vx**2 + pile.pileReaction.Vy**2)**0.5
    beam.add_node_load("N1","FX",V)
    beam.analyze()
    return beam





