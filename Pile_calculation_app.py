import streamlit as st
from handcalcs.decorator import handcalc
import foundation_calculations as fc
from drawingHandler import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from PyNite import FEModel3D
import numpy as np
from concreteproperties.concrete_section import ConcreteSection
import cross_sections as cs
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

# Set main page parameters
st.set_page_config(page_title="Pile designer - beta 0.1", layout="wide")
st.markdown("# Pile Designer - beta version 0.1")

# Main page columns definition
colLeft, colRight = st.columns(2, gap="large")

# Initialize variable in Session_State
if "selectBoxPile" not in st.session_state:
    st.session_state["selectBoxPile"] = "P1"
if "selectBoxSections" not in st.session_state:
    st.session_state["selectBoxSections"] = "Section 1"

# Pilecap geometry information
with st.sidebar:
    st.header("PILECAP GEOMETRY")
    cx = st.number_input("Lateral X cover", value= 0.5)
    cy = st.number_input("Lateral Y cover", value= 0.5)
    h = st.number_input("Depth of pilecap", value= 1.0)
    pileDiameter = st.number_input("Pile diameter",value = 1.0)
    pileLength = st.number_input("Pile length",value = 9)
    rows = st.number_input("Number of rows", value =2)
    columns = st.number_input("Number of columns", value =2)
    Sx = st.number_input("Transversal pile spacing", value =2)
    Sy = st.number_input("Longitudinal pile spacing", value =2)
    if min([cx, cy, h, pileDiameter, rows, columns, Sx, Sy]) < 0: 
        st.error(f"Geometry input values must be positive")
    elif min(Sx,Sy) <= pileDiameter:
        st.error(f"Piles are touching. Invalid geometry")
    st.divider()
    st.header("BALLAST MODULUS")
    #k = st.number_input("Ballast modulus [kN/m3]", value = 1000)
    stepSprings = st.number_input("Distance between springs [m]", value = 0.5)


piles = fc.CalculatePiles(rows,columns,Sx,Sy,pileDiameter,pileLength)
fig = plotCapPile(piles,cx,cy)

#LEFT COLUMN PAGE
with colLeft:
    tabGeometry,tabSoil, tabPierReactions, tabPileEfforts, tabSection, tabVerification = st.tabs(["Geometry","Soil", "Pier reactions","Pile efforts","Cross-section design", "Verification"])
    # 1. Geometry tab
    with tabGeometry:
        st.plotly_chart(fig)
        
    # 2. Soil tab
    with tabSoil:
        columnsSoil = ["name","zEnd","kh","gamma","phi","c"]
        initialSoilValues = [["Layer 1", -5,10000, 0, 0, 0.0],["Layer 2", -100,20000, 0, 0, 0.0]]
        dfSoil = pd.DataFrame(data=initialSoilValues,columns=columnsSoil)
        
        configSoil = {
            "name": st.column_config.TextColumn("Layer name",width="medium", required=True),
            "zEnd": st.column_config.NumberColumn("End of layer (m)",min_value=None,max_value=0),
            "kh": st.column_config.NumberColumn("Kh (kN/m3)",min_value=0,max_value=None),
            "gamma": st.column_config.NumberColumn("γ (kN/m3)",min_value=0,max_value=None),
            "phi": st.column_config.NumberColumn("Φ (degrees)",min_value=0,max_value=None),
            "c": st.column_config.NumberColumn("c' (kN/m2)",min_value=0,max_value=None)
            #"color": st.column_config.SelectboxColumn("Color", options=colorSoil)
        }
        edited_dfSoil = st.data_editor(dfSoil,column_config=configSoil, num_rows="dynamic", disabled=["gamma","phi","c"])

    # 3. Pier reactions tab (from Excel)
    
    with tabPierReactions:
        uxCol, uyCol, uzCol = st.columns(3,gap="medium")
        directions = ["F1", "F2", "F3"]
        with uxCol:
            ux = st.selectbox("Longitudinal direction axis", options=directions, index = 0)
        with uyCol:
            uy = st.selectbox("Trasversal direction axis", options=directions, index = 1)
        with uzCol:
            uz = st.selectbox("Vertical direction axis", options=directions, index = 2)
        globalDirections= {"Ux": ux, "Uy": uy, "Uz": uz, "Rx": ux.replace("F","M"), "Ry": uy.replace("F","M"), "Rz": uz.replace("F","M")}
        updatedFileName = st.file_uploader("Choose Excel file",type="xlsx")
        if updatedFileName is None:
            updatedFileName = "Pier reactions Initial.xlsx" # Excel inicial para que no de error

        if updatedFileName is not None:
            dfPierReactions = fc.ReadExcelPierReactions(updatedFileName, globalDirections)
            pierReactions = {}
            for idx, row in dfPierReactions.iterrows():
                N = row[globalDirections["Uz"]]
                Vlong = row[globalDirections["Ux"]]
                Vtrans = row[globalDirections["Uy"]]
                Mlong = row[globalDirections["Rx"]]
                Mtrans = row[globalDirections["Ry"]]
                T = row[globalDirections["Rz"]]
                if row["StepType"] == "":
                    comboName = row["OutputCase"]
                else:
                    comboName = row["OutputCase"]+" ("+row["StepType"]+")"
                pierReaction = fc.Reaction(comboName, N=N, Mx=Mlong,My=Mtrans,T=T,Vx=Vtrans,Vy=Vlong)
                pierReactions[comboName] = pierReaction

            fc.CalculateSinglePileReactions(piles,pierReactions)
            st.dataframe(dfPierReactions,use_container_width=True, hide_index=True)
        
    # 4. Pile efforts for the different load combinations coming from Excel
    with tabPileEfforts:
        df = pd.DataFrame()
        for pile in piles.values():
            pileReactionsToDisplay = {"Pile": [pile.name for combo, reactions in pile.pileReactions.items()],
                                      "Combination": [combo for combo in pile.pileReactions.keys()],
                                    "N [kN]": [round(reactions.N,2) for combo, reactions in pile.pileReactions.items()],
                                    "Vx [kN]": [round(reactions.Vx,2) for combo, reactions in pile.pileReactions.items()],
                                    "Vy [kN]": [round(reactions.Vy,2) for combo, reactions in pile.pileReactions.items()],
                                    "Vxy [kN]": [round((reactions.Vx**2 + reactions.Vy**2)**0.5,2)  for combo, reactions in pile.pileReactions.items()]}
            index = [i+len(df) for i in range(len(pile.pileReactions))]
            dfAux = pd.DataFrame(pileReactionsToDisplay,index)
            df = pd.concat([df,dfAux],axis=0)
        # Group by pier and pile (only one pier in beta version 0.1 is considered) to get the maximum and minimum N and Vxy per pile
        group = df.groupby(["Pile"], as_index=False)
        #maxPileEfforts = group[["N [kN]", "Vxy [kN]"]].max()
        maxPileEffortsIdx = group[["N [kN]", "Vxy [kN]"]].idxmax()
        #minPileEfforts = group[["N [kN]", "Vxy [kN]"]].min()
        minPileEffortsIdx = group[["N [kN]", "Vxy [kN]"]].idxmin()

        # DataFrames with maximum and minimum N and Vxy and their respective concomitants
        dfMaxN = df.loc[maxPileEffortsIdx["N [kN]"]]
        dfMaxN["Max/Min"] = "Max N"
        dfMaxVxy = df.loc[maxPileEffortsIdx["Vxy [kN]"]]
        dfMaxVxy["Max/Min"] = "Max Vxy"
        dfMinN = df.loc[minPileEffortsIdx["N [kN]"]]
        dfMinN["Max/Min"] = "Min N"
        dfMinVxy = df.loc[minPileEffortsIdx["Vxy [kN]"]]
        dfMinVxy["Max/Min"] = "Min Vxy"

        # write DataFrames
        dfFinal = pd.concat([dfMaxN,dfMinN,dfMaxVxy],axis=0) # I don't concatenate dfMinVxy because makes no sense since Vxy > 0
        st.subheader("Summary of maximum and minimum piles reactions")
        st.dataframe(dfFinal.sort_values(["Pile"],ascending=True),use_container_width=True)
        st.subheader("Total piles reactions")
        st.dataframe(df,use_container_width=True)

        # Update piles dictionary with the combinations that should be calculated (max/min)
        for pile in piles.values():
            maskPileName = dfFinal["Pile"] == pile.name
            dfFiltered = dfFinal.loc[maskPileName]
            for idx, row in dfFiltered.iterrows():
                combo = row["Combination"]
                reaction = pile.pileReactions[combo]
                reaction.isCalculated=True
                reaction.typeOfReaction=row["Max/Min"]

    # 5. Cross-section design
    with tabSection:
        st.subheader("Materials")
        # Define concrete material
        columnsMaterials = ["name","fck/fy","E", "gamma"]
        fck = 25
        Ec = 22000*((fck+8)/10)**0.3
        initialMaterialValues = [["Concrete C30",fck,round(Ec,0),1.5], ["Rebar B500S",500,210000,1.15]]
        dfMaterials = pd.DataFrame(data=initialMaterialValues, index=["Concrete", "Rebar"], columns=columnsMaterials)
        
        configMaterials = {
            "name": st.column_config.TextColumn("Name"),
            "fck/fy": st.column_config.NumberColumn("fck or fy (MPa)",min_value=0,max_value=None),
            "E": st.column_config.NumberColumn("E (MPa)",min_value=0,max_value=None),
            "gamma": st.column_config.NumberColumn("γULS",min_value=0,max_value=None)
        }
            
        edited_dfMaterials = st.data_editor(dfMaterials,column_config=configMaterials, num_rows="static", key="dfMaterials")
        
        # Define pile section parameters
        st.subheader("Pile section definition")
        columnsSections = ["name","D","nBars", "phi","cover"]
        initialSectionValues = [["Section 1",1000*pileDiameter,10,16,50], ["Section 2",1000*pileDiameter,10,20,50]]
        dfSections = pd.DataFrame(data=initialSectionValues, columns=columnsSections)
        
        configSections = {
            "name": st.column_config.TextColumn("Name",required=True),
            "D": st.column_config.NumberColumn("Diameter (mm)",min_value=0,max_value=None,required=True),
            "nBars": st.column_config.NumberColumn("Nº of rebars",min_value=0,max_value=None, required=True),
            "phi": st.column_config.SelectboxColumn("Φ (mm)",options=[10,12,16,20,25,32],required=True),
            "cover": st.column_config.NumberColumn("Cover (mm)",min_value=0,max_value=None,required=True)
        }
        edited_dfSections = st.data_editor(dfSections,column_config=configSections, num_rows="dynamic")

        # Create pile section objects and plots
        pileSections = {}
        figSections, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,5))
        for idx, row in edited_dfSections.iterrows():
            sec = cs.createPileSection(row,edited_dfMaterials)
            pileSections[row["name"]] = sec
        sectionOptions = [sectionName for sectionName in pileSections.keys()]
        selectedPileSectionName = st.selectbox("Select pile section", options=sectionOptions, key="selectBoxSections")
        axSection = pileSections[selectedPileSectionName].plot_section(title=f"Geometry ({selectedPileSectionName})",ax=axes[0])
        axSection.legend().set_visible(False)

        momentInteractionDiagram = pileSections[selectedPileSectionName].moment_interaction_diagram(progress_bar=False)
        axInteraction = momentInteractionDiagram.plot_diagram(ax=axes[1])

        pileSelectedName = st.session_state["selectBoxPile"]
        #pileSelectedName = "P1" # Used only to debug
        Nmin = dfFinal["N [kN]"][(dfFinal["Pile"] == pileSelectedName) & (dfFinal["Max/Min"] == "Min N")].iloc[0]
        Nmax = dfFinal["N [kN]"][(dfFinal["Pile"] == pileSelectedName) & (dfFinal["Max/Min"] == "Max N")].iloc[0]

        N = [Nmin, Nmax] # Nmin, Nmax of selected pile
        Mrd = {} # dictionary of bending moment resistance for Nmin and Nmax
        for section in sectionOptions:
            for n in N:
                ultimateBendingCapacityNmin = pileSections[section].ultimate_bending_capacity(n=N[0]*1e3) # Nmin
                ultimateBendingCapacityNax = pileSections[section].ultimate_bending_capacity(n=N[1]*1e3) # Nmax
                Mrd[section] = {"Nmin": ultimateBendingCapacityNmin.m_xy/1e6, "Nmax": ultimateBendingCapacityNax.m_xy/1e6}
        
        colMrdNmax, colMrdNmin = st.columns(2)
        with colMrdNmin:
            MrdNmin = Mrd[selectedPileSectionName]["Nmin"]
            st.metric(label=f"Bending resistance moment  for Nmin = {N[0]} kN", value=f"{round(MrdNmin,2)} kNm")
            #st.pyplot(axSection.figure)
        with colMrdNmax:
            MrdNmax = Mrd[selectedPileSectionName]["Nmax"]
            st.metric(label=f"Bending resistance moment  for Nmax = {N[1]} kN", value=f"{round(MrdNmax,2)} kNm")
            #st.pyplot(axInteraction.figure)
        st.pyplot(figSections)

    with tabVerification:
        columnsVerification = ["name","L", "section"]
        initialVerificationValues = [["Reinforcement 1",pileLength/3, sectionOptions[0]], ["Reinforcement 2",pileLength/3, sectionOptions[1]],["Reinforcement 3",pileLength/3, sectionOptions[0]]]
        dfVerification = pd.DataFrame(data=initialVerificationValues, columns=columnsVerification)
        configVerification = {
        "name": st.column_config.TextColumn("Name",required=True),
        "L": st.column_config.NumberColumn("Length (m)",min_value=0,max_value=None,required=True),
        "section": st.column_config.SelectboxColumn("Pile Section",options=sectionOptions,required=True)
                    }
        edited_dfVerification = st.data_editor(dfVerification,column_config=configVerification, num_rows="dynamic",use_container_width=True)


# Create pile model with PyNite and plot moment, shear and displacement of selected pile
pileSelectedName = st.session_state["selectBoxPile"]
#pileSelectedName = "P1" # Only for debug
for pile in piles.values():
    E = edited_dfMaterials["E"]["Concrete"]
    nu = 0.2
    G = E/(2*(1+nu))
    fck = edited_dfMaterials["fck/fy"]["Concrete"]
    concreteName = edited_dfMaterials["name"]["Concrete"]
    pile.concrete = fc.Concrete(concreteName,E,G,nu,25,fck)

selectedPile = piles[pileSelectedName]
pileModel = fc.CalculatePileFEModel(selectedPile,soil=edited_dfSoil, stepSprings = stepSprings)
nPoints = 1000
figEfforts = plotPileEffortsSolution(pileModel,selectedPile,edited_dfSoil,Mrd,edited_dfVerification,nPoints)

        
# RIGHT COLUMN PAGE (PLOTS)
with colRight:
    momentCol, shearCol, dispCol = st.columns(3)
    pileSelectedName = st.selectbox("Select pile", options=[pile.name for pile in piles.values()], key="selectBoxPile")
    selectedPile = piles[pileSelectedName]

    #addSelectedPiletoExistingPlot(fig,selectedPile) #  No dibuja el pilote seleccionado y no sé porqué

    #pileModel = fc.CalculatePileFEModel(selectedPile,soil=edited_dfSoil, stepSprings = stepSprings)
    #figEfforts = plotPileEffortsSolution(pileModel,edited_dfSoil,nPoints)

    st.plotly_chart(figEfforts)

    # Metrics of maximum and minimum moment, shear and displacement of all 4 load combinations.
    memberName = list(pileModel.Members)[0]
    momentAux = pileModel.Members[memberName].moment_array(Direction="Mz",combo_name=list(pileModel.LoadCombos)[0], n_points=nPoints)
    dimension = len(momentAux[0])
    momentsPile, shearsPile, displacementsPile = np.zeros([2,dimension]), np.zeros([2,dimension]), np.zeros([2,dimension])
    for combo in list(pileModel.LoadCombos):
        moment = pileModel.Members[memberName].moment_array(Direction="Mz",combo_name=combo, n_points=nPoints)
        shear = pileModel.Members[memberName].shear_array(Direction="Fy",combo_name=combo,n_points=nPoints)
        displacement= pileModel.Members[memberName].deflection_array(Direction="dy", combo_name=combo, n_points=nPoints)
        momentsPile = np.concatenate([momentsPile,moment], axis=1)
        shearsPile = np.concatenate([shearsPile,shear],axis=1)
        displacementsPile = np.concatenate([displacementsPile,displacement],axis=1)
    with momentCol:
        st.metric(label="Moment Max", value=f"{round(max(momentsPile[1]),2)} kNm")
        st.metric(label="Moment Min", value=f"{round(min(momentsPile[1]),2)} kNm")
    with shearCol:
        st.metric(label="Shear Max", value=f"{round(max(shearsPile[1]),2)} kN")
        st.metric(label="Shear Min", value=f"{round(min(shearsPile[1]),2)} kN")
    with dispCol:
        st.metric(label="Displacement Max", value=f"{round(1000*max(displacementsPile[1]),2)} mm")
        st.metric(label="Displacement Min", value=f"{round(1000*min(displacementsPile[1]),2)} mm")

    









