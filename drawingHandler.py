import plotly.graph_objects as go
from plotly.subplots import make_subplots
import foundation_calculations as fc
import pandas as pd
import numpy as np
from PyNite import FEModel3D
import streamlit as st

def plotCapPile(piles: dict[fc.Pile], cx: float, cy: float)-> go.Figure:
    pileDiameter = piles["P1"].diameter
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

    return fig

def plotPileEffortsSolution(pileModel: FEModel3D, selectedPile: fc.Pile, edited_dfSoil: pd.DataFrame, mrd: dict, edited_dfVerifications: pd.DataFrame, nPoints: int)-> go.Figure:
    
    # Add soil layers to the effort plots
    figEfforts = make_subplots(rows=1, cols=3,subplot_titles=("MOMENT", "SHEAR", "DISPLACEMENT"), print_grid=True)
    maxLimit = 1e10
    minLimit = -maxLimit
    colorSoil= ["cornsilk","goldenrod","darkgoldenrod","plum","violet","purple","darkmagenta"]

    for numberOfPlot in range(3):
        # Soils plot
        for idx in range(len(edited_dfSoil)):
            if idx>0:
                zSoilTop = edited_dfSoil.loc[idx-1,"zEnd"]
            else:
                zSoilTop = 0
            zSoilBot = edited_dfSoil.loc[idx,"zEnd"]
            soilName = edited_dfSoil.loc[idx,"name"]

            figEfforts.add_shape(type="rect",
                x0=minLimit, y0=zSoilBot, x1=maxLimit, y1=zSoilTop,
                line=dict(
                    color="black",
                    width=1,
                ),
                fillcolor=colorSoil[idx],
                opacity=0.6,
                layer="below",
                name=soilName,
                showlegend=True,
                row=1, col=numberOfPlot+1
            )

    # Add moment, shear and displacement solution plot
    translateToComboType = {}
    for k, v in selectedPile.pileReactions.items():
        translateToComboType[k] = v.typeOfReaction

    momentLimits = []
    shearLimits = []
    displacementLimits = []
    memberName = list(pileModel.Members)[0]
    for combo in pileModel.LoadCombos.keys():
        simplifiedComboLabel = translateToComboType[combo]
        momentPile = pileModel.Members[memberName].moment_array(Direction="Mz",combo_name=combo, n_points=nPoints)
        shearPile = pileModel.Members[memberName].shear_array(Direction="Fy",combo_name=combo,n_points=nPoints)
        displacementPile = pileModel.Members[memberName].deflection_array(Direction="dy", combo_name=combo, n_points=nPoints)
        # Calculate Max and Min limits for plotting
        momentLimits.append(max(momentPile[1]))
        momentLimits.append(min(momentPile[1]))
        shearLimits.append(max(shearPile[1]))
        shearLimits.append(min(shearPile[1]))
        displacementLimits.append(max(displacementPile[1]))
        displacementLimits.append(min(displacementPile[1]))

        figEfforts.add_trace(
            go.Scatter(
                x = momentPile[1],
                y = -momentPile[0],
                line = {"color": "black"},
                mode = "lines",
                name=simplifiedComboLabel,
                line_width = 2,
                showlegend=True,
                fill="tonextx"
            ), row=1,col=1
        ) 
        figEfforts.add_trace(
            go.Scatter(
                x = shearPile[1],
                y = -shearPile[0],
                line = {"color": "red"},
                mode = "lines",
                name=simplifiedComboLabel,
                line_width = 2,
                showlegend=True,
                fill="tonextx"
            ), row=1,col=2
        ) 
        figEfforts.add_trace(
            go.Scatter(
                x = 1000*displacementPile[1],
                y = -displacementPile[0],
                line = {"color": "blue"},
                mode = "lines",
                name=simplifiedComboLabel,
                line_width = 2,
                showlegend=True
            ), row=1,col=3
        ) 
    
    # Add resistance moments into the plot
    xMrdPos, xMrdNeg = [], []
    yMrd = []
    for idx, row in edited_dfVerifications.iterrows():
        sectionName = row["section"]
        mRdNmin, mRdNmax = mrd[sectionName]["Nmin"], mrd[sectionName]["Nmax"]
        mrdPlot = min(mRdNmin,mRdNmax)
        xMrdPos.append(mrdPlot)
        xMrdNeg.append(-mrdPlot)
        if idx == 0:
            y = 0
        else:
            y = -yMrd[idx-1] + edited_dfVerifications["L"][idx-1]
        
        yMrd.append(-y)

    xMrdPos.append(mrdPlot)
    xMrdNeg.append(-mrdPlot)
    yMrd.append(minLimit)
    mRdLabel = ["Mrd+", "Mrd-"]
    for idx, xMrd in enumerate([xMrdPos, xMrdNeg]):
        figEfforts.add_trace(
            go.Scatter(
                x = xMrd,
                y = yMrd,
                line = {"color": "red", "dash": "dash"},
                mode = "lines",
                name=mRdLabel[idx],
                line_width = 2,
                line_shape = "vh",
                showlegend=True
            ), row=1,col=1
        ) 

    figEfforts.layout.width = 800
    figEfforts.layout.height = 700
    figEfforts.layout.yaxis.title = "Pile height [m]"
    amplification = 1.3
    figEfforts.update_xaxes(title_text="Moment [kNm]", showgrid =True, range = [amplification*min(min(momentLimits),min(xMrdNeg)),amplification*max(max(momentLimits),max(xMrdPos))], row=1, col=1)
    figEfforts.update_xaxes(title_text="Shear [kN]", showgrid =True, range = [amplification*min(shearLimits),amplification*max(shearLimits)], row=1, col=2)
    figEfforts.update_xaxes(title_text="Displacement [mm]", showgrid =True, range = [1000*amplification*min(displacementLimits),1000*amplification*max(displacementLimits)], row=1, col=3)
    pileLength = -max(momentPile[0])
    for i in range(1,4,1):
        figEfforts.update_yaxes(range=[pileLength,0],row=1,col=i)

    return figEfforts

def addSelectedPiletoExistingPlot(fig: go.Figure, pile: fc.Pile):
    fig.add_shape(showlegend=True,
        type="circle",
        xref="x", yref="y",
        x0=pile.coordinates[0]-pile.diameter/2, y0=pile.coordinates[1]-pile.diameter/2, x1=pile.coordinates[0]+pile.diameter/2, y1=pile.coordinates[1]+pile.diameter/2,
        name = pile.name,
        line_width = 2,
        line_color="Red",
        #fillcolor="Grey",
        label = dict(text=f"{pile.name}")
        )