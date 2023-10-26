import streamlit as st
from handcalcs.decorator import handcalc
import foundation_calculations as fc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from PyNite import FEModel3D
import numpy as np

st.set_page_config(page_title="Pile designer - beta 0.1", layout="wide")
st.markdown("# Pile Designer - beta version 0.1")


colSideBar1, colSideBar2 =st.sidebar.columns(2)
colLeft, colRight = st.columns(2, gap="large")

# Pilecap geometry information
with colSideBar1:
    st.header("PILECAP GEOMETRY")
    cx = st.number_input("Lateral X cover", value= 0.5)
    cy = st.number_input("Lateral Y cover", value= 0.5)
    h = st.number_input("Depth of pilecap", value= 1.0)
    pileDiameter = st.number_input("Pile diameter",value = 1.0)
    pileLength = st.number_input("Pile length",value = 5)
    rows = st.number_input("Number of rows", value =2)
    columns = st.number_input("Number of columns", value =2)
    Sx = st.number_input("Transversal pile spacing", value =1)
    Sy = st.number_input("Longitudinal pile spacing", value =1)
    if min([cx, cy, h, pileDiameter, rows, columns, Sx, Sy]) < 0: 
        st.error(f"Geometry input values must be positive")
    elif min(Sx,Sy) <= pileDiameter:
        st.error(f"Piles are touching. Invalid geometry")

# Pier reactions information
with colSideBar2:
    st.header("PIER REACTIONS")
    N = st.number_input("Axial force [kN]", value = 100)
    Mlong = st.number_input("Longitudinal moment [kNm]", value = 100)
    Mtrans = st.number_input("Transversal moment [kNm]", value = 100)
    Vlong = st.number_input("Longitudinal shear [kN]", value = 100)
    Vtrans = st.number_input("Transversal shear [kN]", value = 100)
    T = st.number_input("Torsor moment [kNm]", value = 100)
    st.divider()
    st.header("BALLAST MODULUS")
    k = st.number_input("Ballast modulus [kN/m3]", value = 1000)
    stepSprings = st.number_input("Distance between springs [m]", value = 0.1)
    pierReaction = fc.Reaction(N=N, Mx=Mlong,My=Mtrans,T=T,Vx=Vtrans,Vy=Vlong)

st.sidebar.divider()

# Pile centers calculation
piles = fc.CalculatePiles(rows,columns,Sx,Sy,pileDiameter,pileLength,pierReaction)
fc.CalculateSinglePileReactions(piles)
x, y , pileNames = [], [], []
for pile in piles.values():
    x.append(pile.coordinates[0])
    y.append(pile.coordinates[1])
    pileNames.append(pile.name)

xPileMax, xPileMin = max(x) + pileDiameter/2 + cx, min(x) - pileDiameter/2 - cx
yPileMax, yPileMin = max(y) + pileDiameter/2 + cy, min(y) - pileDiameter/2 - cy

xCapPile = [xPileMin, xPileMax, xPileMax, xPileMin, xPileMin]
yCapPile = [yPileMin, yPileMin, yPileMax, yPileMax, yPileMin]

# Plot 
fig = go.Figure()
 # Plot pile center coordinates
fig.add_trace(
    go.Scatter(
        x = x,
        y = y,
        mode = "text",
        name = "Piles",
        text= pileNames,
        textfont=dict(color="red"),
        textposition = "middle center"
    )
) 
fig.add_trace(
    go.Scatter(
        x = xCapPile,
        y = yCapPile,
        line = {"color": "black"},
        mode = "lines",
        name="Pilecap",
        #fill = "toself", fillcolor="lightgrey",
        line_width = 1,
        showlegend=False
    )
)
# Plot pile circles and coordinates
for pile in piles.values():
    # Plot pile circles
    fig.add_shape(showlegend=False,
        type="circle",
        xref="x", yref="y",
        x0=pile.coordinates[0]-pileDiameter/2, y0=pile.coordinates[1]-pileDiameter/2, x1=pile.coordinates[0]+pileDiameter/2, y1=pile.coordinates[1]+pileDiameter/2,
        name = pile.name,
        line_width = 1,
        line_color="Black",
        #fillcolor="Grey",
        label = dict(text=f"{pile.name}")
    )

   
fig.layout.xaxis.scaleanchor = "y"
fig.layout.xaxis.scaleratio = 1
fig.layout.xaxis.title = "Transversal X [m]"
fig.layout.yaxis.title = "Longitudinal Y [m]"

with colLeft:
    st.header("Geometry")
    pileSelectedName = st.selectbox("Select pile", [pile.name for pile in piles.values()])
    #st.plotly_chart(fig)
    selectedPile = piles[pileSelectedName]
    fig.add_shape(showlegend=True,
        type="circle",
        xref="x", yref="y",
        x0=selectedPile.coordinates[0]-selectedPile.diameter/2, y0=selectedPile.coordinates[1]-selectedPile.diameter/2, x1=selectedPile.coordinates[0]+selectedPile.diameter/2, y1=selectedPile.coordinates[1]+selectedPile.diameter/2,
        name = selectedPile.name,
        line_width = 2,
        line_color="Red",
        #fillcolor="Grey",
        label = dict(text=f"{selectedPile.name}")
    )
    st.plotly_chart(fig)
    pileReactionsToDisplay = { "N [kN]": [round(pile.pileReaction.N,2) for pile in piles.values()],
                              "Vx [kN]": [round(pile.pileReaction.Vx,2) for pile in piles.values()],
                              "Vy [kN]": [round(pile.pileReaction.Vy,2) for pile in piles.values()],
                              "Vxy [kN]": [round((pile.pileReaction.Vx**2 + pile.pileReaction.Vy**2)**0.5,2)  for pile in piles.values()]}
    index = [pile.name for pile in piles.values()]
    df = pd.DataFrame(pileReactionsToDisplay,index)
    st.dataframe(df,use_container_width=True)

    with colRight:
        nPoints = 1000
        pileModel = fc.CalculatePileFEModel(selectedPile,k=k, stepSprings = stepSprings)
 
        momentPile = pileModel.Members[selectedPile.name].moment_array(Direction="Mz",n_points=nPoints)
        shearPile = pileModel.Members[selectedPile.name].shear_array(Direction="Fy",n_points=nPoints)
        displacementPile = pileModel.Members[selectedPile.name].deflection_array(Direction="dy", n_points=nPoints)
    
        figEfforts = make_subplots(rows=1, cols=3,subplot_titles=("MOMENT", "SHEAR", "DISPLACEMENT"), print_grid=True)
        figPileDisplacement = go.Figure()
        figEfforts.add_trace(
            go.Scatter(
                x = momentPile[1],
                y = -momentPile[0],
                line = {"color": "black"},
                mode = "lines",
                name="Moment",
                line_width = 2,
                showlegend=True
            ), row=1,col=1
        ) 
        figEfforts.add_trace(
            go.Scatter(
                x = shearPile[1],
                y = -shearPile[0],
                line = {"color": "red"},
                mode = "lines",
                name="Shear",
                line_width = 2,
                showlegend=True
            ), row=1,col=2
        ) 
        figEfforts.add_trace(
            go.Scatter(
                x = 1000*displacementPile[1],
                y = -displacementPile[0],
                line = {"color": "blue"},
                mode = "lines",
                name="Displacement",
                line_width = 2,
                showlegend=True
            ), row=1,col=3
        ) 
        figEfforts.layout.width = 700
        figEfforts.layout.height = 700
        figEfforts.layout.yaxis.title = "Pile height [m]"
        amplification = 1.2
        figEfforts.update_xaxes(title_text="Moment [kNm]", showgrid =True, range = [amplification*min(momentPile[1]),amplification*max(momentPile[1])], row=1, col=1)
        figEfforts.update_xaxes(title_text="Shear [kN]", showgrid =True, range = [amplification*min(shearPile[1]),amplification*max(shearPile[1])], row=1, col=2)
        figEfforts.update_xaxes(title_text="Displacement [mm]", showgrid =True, range = [1000*amplification*min(displacementPile[1]),1000*amplification*max(displacementPile[1])], row=1, col=3)
        st.plotly_chart(figEfforts)






