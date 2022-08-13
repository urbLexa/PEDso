# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:52:37 2021

@author: abruck
"""

from pyomo.environ import *
import pandas as pd
import numpy as np
import helpers as hp
import evClassKm



def build_model(in_tech = "Input_year_Technology.xlsx", in_loc = "Input_year_Location.xlsx", in_timeseries = "Input_year_Timeseries.xlsx",\
                in_emisFules = "Input_year_Emissions_Fuels.xlsx", in_tariffs = "Input_year_Tariffs.xlsx", in_ev = "Input_year_EV.xlsx",\
                    in_other = "Input_year_Other.xlsx", projectlife = 20, dt = 60, consHeat = "No", consCool = "No", consEV = "No", kW2X = 1):
    """
    

    Parameters
    ----------
    in_tech : TYPE, optional
        DESCRIPTION. The default is "Input_year_Technology.xlsx".
    in_loc : TYPE, optional
        DESCRIPTION. The default is "Input_year_Location.xlsx".
    in_timeseries : TYPE, optional
        DESCRIPTION. The default is "Input_year_Timeseries.xlsx".
    in_emisFules : TYPE, optional
        DESCRIPTION. The default is "Input_year_Emissions_Fuels.xlsx".
    in_tariffs : TYPE, optional
        DESCRIPTION. The default is "Input_year_Tariffs.xlsx".
    in_ev : TYPE, optional
        DESCRIPTION. The default is "Input_year_EV.xlsx".
    in_other : TYPE, optional
        DESCRIPTION. The default is "Input_year_Other.xlsx".
    projectlife : TYPE, optional
        DESCRIPTION. The default is 20.
    dt : TYPE, optional
        DESCRIPTION. The default is 60.
    consHeat : TYPE, optional
        DESCRIPTION. The default is "No".
    consCool : TYPE, optional
        DESCRIPTION. The default is "No".
    consEV : TYPE, optional
        DESCRIPTION. The default is "No".
    kW2X : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    dth = dt/60 #timestep in hours
    

    timesteps = 365 * 24 * 60 / dt
    date_rng = pd.date_range(start='2019-01-01 00:00', end='2019-12-31 23:59', freq = str(dt)+('min'))
    date_rng_hourly = pd.date_range(start='2019-01-01 00:30', end='2019-12-31 23:59', freq = "H")

    ts_i = np.arange(timesteps, dtype = int)
    
    
    timesteps_wholeyear = 365 * 24 * 60 / dt
    

    """
    INPUT DATA SECTION
    """
    
    
    """
    LOAD(S)
    """
    loadData = pd.read_excel(in_timeseries, sheet_name = "LoadData")
    loadChangeData = pd.read_excel(in_timeseries, sheet_name = "Load_Change", index_col="Variable")
    load_e_change = loadChangeData.at["load_e_change","Value"]
    load_h_change = loadChangeData.at["load_h_change","Value"]
    load_c_change = loadChangeData.at["load_c_change","Value"]
    load_e = (loadData.loc[:,"Load_e [kW]"].values)/kW2X
    if consHeat == "Yes":
        load_h = (loadData.loc[:,"Load_h [kW]"].values)/kW2X
    if consCool == "Yes":
        load_c = (loadData.loc[:,"Load_c [kW]"].values)/kW2X


    
    """
    TARIFFS
    """
    tariff_changes = pd.read_excel(in_tariffs, sheet_name = "Tariff_el_change", index_col="Variable")
    kWdia_change_rate = tariff_changes.at["kWdia_change_percent","Value"]
    tar_feedin_change_rate = tariff_changes.at["tar_feedin_change_percent","Value"]
    tariff_change_rate = tariff_changes.at["tariff_change_percent","Value"]
    

    tariffs = pd.read_excel(in_tariffs, sheet_name = "Tariff_el")
    if tariffs.at[0,"Tariff [€/kW]"] != None:
        tar_kWdia = tariffs.at[0,"Tariff [€/kW]"] * kW2X
    else:
        tar_kWdia = 0        
    tariffs = tariffs.drop(["Tariff [€/kW]"], axis =1)
    
    if tariffs["Tariff_feed [€/kWh]"].isnull().values.any():
        tar_feedin = tariffs.at[0,"Tariff_feed [€/kWh]"] * kW2X
        tariffs = tariffs.drop(["Tariff_feed [€/kWh]"], axis =1)
    tariffs_dt = hp.tariffsizer(tariffs, timesteps, dt)
    if "Tariff_feed [€/kWh]" in tariffs_dt.columns:
        tar_feedin = (tariffs_dt.loc[:,"Tariff_feed [€/kWh]"].values) * kW2X
    tariff = (tariffs_dt.loc[:,"Tariff [€/kWh]"].values) * kW2X
    if all(tariff == tariff[0]) == True:
        tariff = tariff[0]
    
    
    """
    TECHNOLOGY SELECTION
    """
    ders_total = pd.read_excel(in_tech, sheet_name = "Selection&Cost")
    ders_df = ders_total[ders_total["Consider"] == 1]
    ders_df = ders_df.set_index(['DER'], drop = False)
    
    #All distributed energy sources that can output electricity
    d_elOut = ders_df[ders_df['el_out'] == 1].index.values
    
    #All considered tech (including grid) that can output and take in electrictiy
    elOut = np.append(d_elOut, 'grid')
    elIn = ders_df[ders_df['el_in'] == 1].index.values
    elIn = np.append(elIn, ['grid', 'load_e'])
    
    #All considered tech that can output and take in heat
    hOut = ders_df[ders_df['heat_out'] == 1].index.values
    hIn = ders_df[ders_df['heat_in'] == 1].index.values
    hIn = np.append(hIn, ['load_h'])
    
    #All considered tech that can output and take in cooling
    cOut = ders_df[ders_df['cool_out'] == 1].index.values
    cIn = ders_df[ders_df['cool_in'] == 1].index.values
    cIn = np.append(cIn, ['load_c'])

    # All considered tech that generates electricity
    d_elGen = d_elOut[d_elOut!='battery']
    
    #fuel and emission technology frames
    fueltec = ders_df[ders_df['fuels_in'] == 1].index.values
    bfueltec = ders_df[ders_df['fuel_type'] == "bf"].index.values
    ffueltec = ders_df[ders_df['fuel_type'] == "ff"].index.values
    emistec = np.append(fueltec, 'grid')
    
    
    """
    OTHER TECH-DATA
    """
    #PV
    pv_facts = pd.read_excel(in_tech, sheet_name = "pv", index_col="Variable")
    pv_eff = pv_facts.at["efficiency","Value"] #Panel efficiency
    pv_pr = pv_facts.at["pr","Value"] #Performance ratio 
    pv_gcr_ground = pv_facts.at["gcr_ground","Value"] # ground coverage ratio flat roof/ground
    pv_gcr_roof = pv_facts.at["gcr_roof","Value"] # ground coverage ratio roof
    pv_gcr_wall = pv_facts.at["gcr_Wall","Value"] # ground coverage ratio wall
    pv_curt = pv_facts.at["curtailment","Value"] # curtailment allowed or not
    if pv_curt == 1:
        elIn = np.append(elIn, 'curtail')
        
    #Battery
    battery_facts = pd.read_excel(in_tech, sheet_name = "battery", index_col="Variable")
    batt_capa2pow = battery_facts.at["cap2power","Value"] # Relationship of power to capacity --> power = capa * capa2pow
    batt_socIni = battery_facts.at["SOC_ini","Value"] # Starting State of charge of battery in kWh
    batt_eff = battery_facts.at["efficiency","Value"] # Round trip efficiency of the battery
    batt_grid2batt = battery_facts.at["grid2batt","Value"] # if 0 no grid to batt allowed
    
    #Electric Boiler
    elBoiler_facts = pd.read_excel(in_tech, sheet_name = "boiler_electric", index_col="Variable")
    elBoiler_eff = elBoiler_facts.at["efficiency","Value"] #Boiler efficiency
    
    #Gas Boiler
    gasBoiler_facts = pd.read_excel(in_tech, sheet_name = "boiler_gas", index_col="Variable")
    gasBoiler_eff = gasBoiler_facts.at["efficiency","Value"] #Boiler efficiency
    
    #Thermal Storage System
    tss_facts = pd.read_excel(in_tech, sheet_name = "tss", index_col="Variable")
    tss_eff = tss_facts.at["efficiency","Value"] #TSS roundtrip efficiency
    tss_socIni = tss_facts.at["SOC_ini","Value"] # initial state of charge
    tss_capa2pow = tss_facts.at["cap2power","Value"] # Relationship of power to capacity
    tss_self = tss_facts.at["selfdis","Value"] # selfdischarge per hour
    
    #Micro CHP Plant
    mchp_facts = pd.read_excel(in_tech, sheet_name = "mCHP", index_col="Variable")
    mchp_eff_e = mchp_facts.at["efficiency_e","Value"] #electric efficiency 
    mchp_eff_h = mchp_facts.at["efficiency_h","Value"] #heating efficiency 
    mchp_p2h = mchp_eff_h / mchp_eff_e
    mchp_mincap = mchp_facts.at["minCapa","Value"] #minimal Capacity to be used 
    
    # Wind turbine 
    wt_facts = pd.read_excel(in_tech, sheet_name = "wind", index_col="Variable")
    wt_cutInWS = wt_facts.at["cutInWS","Value"] #cut in wind speed
    wt_ratedWS = wt_facts.at["ratedWS","Value"] #rated wind speed
    wt_cutOutWS = wt_facts.at["cutOutWS","Value"] #cut out wind speed
    wt_ratedP = wt_facts.at["ratedPower","Value"] #rated power
    wt_linear = hp.lineequ(wt_cutInWS, 0, wt_ratedWS, wt_ratedP)
    wt_m = wt_linear[0]
    wt_t = wt_linear[1]
    
    #Solar Thermal
    st_facts = pd.read_excel(in_tech, sheet_name = "solarThermal", index_col="Variable")
    st_eff_0 = st_facts.at["eff_0","Value"] # Optical efficiency
    st_a1 = st_facts.at["a1","Value"] # Thermal loss factor 1
    st_a2 = st_facts.at["a2","Value"] # Thermal loss factor 2
    st_Tf = st_facts.at["Tf","Value"] # mean heating fluid temperature
    
    #Absorption Chiller
    absCh_facts = pd.read_excel(in_tech, sheet_name = "abs_chiller", index_col="Variable")
    absCOP = absCh_facts.at["COP","Value"] # COP
    
    #Ground-Source Heat Pump
    gshp_facts = pd.read_excel(in_tech, sheet_name = "gshp", index_col="Variable")
    gshpCOP_h = gshp_facts.at["COPh","Value"] # COP heating mode constant
    gshpCOP_c = gshp_facts.at["COPc","Value"] # COP cooling mode constant
    
    """
    OTHER INPUT DATA
    """
    
    other = pd.read_excel(in_other, sheet_name = "Other", index_col="Variable")
    
    i_rate = other.at["i_rate", "Value"] #interest rate of investment
    selfConsumption = other.at["min_selfCons", "Value"] #minimum value for self consumption
    pedBalance = other.at["pedBalance", "Value"]        #None, ElectricitPEF, PEF_dynamic; What type of balancing is expected
    primaryEn_fac_el = other.at["primaryEn_fac_el", "Value"]    # 2 https://ec.europa.eu/energy/sites/ener/files/documents/final_report_pef_eed.pdf; 
    #standard value for primary energy factor for electricity
    fuel_export_in =  other.at["fuel_export", "Value"] # Annual fuel export in kWh of higher heating value
    fuel_export_change_ann =  other.at["fuel_export_change_ann", "Value"] # Annual fuel export change 
    fuel_export = {}
       
    
    
    lhv_factor = other.at["lhv_factor", "Value"] #https://eurogas.org/website/wp-content/uploads/2018/03/Statistics_2011_091211.pdf --> 0.9 conversion factor HHV/LHV

    grid_limit = other.at["gridlim", "Value"]
    grid_limitValue = other.at["gridLimitValue", "Value"] / kW2X

    pe_factor_df = pd.read_excel(in_other, sheet_name = "pe_factor", index_col="Resource")
    #https://ec.europa.eu/energy/sites/ener/files/documents/final_report_pef_eed.pdf
    pe_factor = pe_factor_df.T.to_dict("records")[0]
    #primary energy factor electricity according to generation fuel/technology

    merit_order_df = pd.read_excel(in_other, sheet_name = "Merit_Order", index_col="Order")
    merit_order = merit_order_df.T.to_dict("records")[0]
    #merit order of electricity generation technology
    
    loc_df = pd.read_excel(in_loc, sheet_name = "Location", index_col="Variable")
    loc = loc_df.T.to_dict("records")[0] #location dict with lat, long, alt, temp and alb

    
    pv_area_g = pd.read_excel(in_loc, sheet_name = "pv_ground_flatroof", index_col="Variable")
    pv_groundareaAvail = pv_area_g.at["pv_groundareaAvail","Value"] #free area ground or flat roofs in m2
    pv_roofwallareaAvail = pd.read_excel(in_loc, sheet_name="pv_tiltroof_wall",\
                                         index_col="Az") #free area tilt roof or wall  in m2
    
    pv_ground_exist = pd.read_excel(in_loc, sheet_name = "pv_exist_flat",\
                                    index_col = "Az")
    pv_ground_exist = pv_ground_exist.squeeze().to_dict()
    pv_ground_exist = {x:y for x,y in pv_ground_exist.items() if y!=0}
    if bool(pv_ground_exist) == False:
        pv_ground_exist[0] = 0
    
    pv_roof_exist = pd.read_excel(in_loc, sheet_name = "pv_exist_tilt",\
                                  index_col = "Az")
    # pv_roof_exist = pv_roof_exist[(pv_roof_exist.T != 0).any()]
    # pv_roof_exist = pv_roof_exist.loc[:,(pv_roof_exist != 0).any(axis = 0)]
    # if pv_roof_exist.empty == True:
    #     pv_roof_exist = pd.DataFrame(index = [0], columns = [15])
    #     pv_roof_exist = pv_roof_exist.fillna(0)
    
    pv_wall_exist = pd.read_excel(in_loc, sheet_name = "pv_exist_wall",\
                                    index_col = "Az")
    pv_wall_exist = pv_wall_exist.squeeze().to_dict()
    pv_wall_exist = {x:y for x,y in pv_wall_exist.items() if y!=0}
    if bool(pv_wall_exist) == False:
        pv_wall_exist[0] = 0
     
    #Fuel Selection   
    fuel_df = pd.read_excel(in_emisFules, sheet_name = "Fuels", index_col="Type")
    ff = fuel_df.at["ff", "Fuel"] #Fossil Fuel
    bf = fuel_df.at["bf", "Fuel"] #Bio Fuel
    
    #Higher heating value
    hhv_df = pd.read_excel(in_emisFules, sheet_name = "hhv", index_col="Fuel")
    hhv = hhv_df.T.to_dict("records")[0]
    # kWh/m3 hhv --> depends on gas quality --> https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:52018XC0614(02)&rid=4
    #https://ec.europa.eu/energy/intelligent/projects/sites/iee-projects/files/projects/documents/redubar_a_register_of_all_gas_regulations.pdf
    #kWh/m3 hhv --> depends on origin and if scrubbed

    #fuel costs (fixed and variable)
    fuelcost_df = pd.read_excel(in_emisFules, sheet_name = "CostFuels", index_col="Fuel")
    cost_fuels = fuelcost_df.T.to_dict("records")[0]
    
    
    #Emission Factors
    emissions_df = pd.read_excel(in_emisFules, sheet_name = "CO2_factors", index_col="Fuel")
    gridEmissRed = emissions_df.at["grid", "Reduction"]
    emissions_df = emissions_df.drop(["Reduction"], axis = 1)
    emissions_df = emissions_df* kW2X
    emiss_fac = emissions_df.T.to_dict("records")[0]  #kg/kWh
    
    #Emission Decisions
    emissDecisions_df = pd.read_excel(in_emisFules, sheet_name = "CO2", index_col="Variable")
    CO2_cap = emissDecisions_df.at["CO2_cap","Value"] # if true cap needs to be specified --> PED 0!
    CO2_cap_value = emissDecisions_df.at["CO2_cap_value","Value"]
    cap_change_ann = emissDecisions_df.at["cap_change_ann","Value"]
    CO2_cost = emissDecisions_df.at["CO2_cost","Value"] # decides if CO2 cost should be considered or not
    CO2_price = emissDecisions_df.at["CO2_price","Value"] # €/ton of CO2
    CO2_price_increase_ann = emissDecisions_df.at["price_increase_ann","Value"] # annual increase of price in %
    
    #Grid Mix data
    gen = pd.read_excel(in_timeseries, sheet_name = "GridMixData")
    g_tot = gen.total
    g_share = gen.drop("total", axis = 1).div(g_tot, axis = 0)
    for i in g_share.columns:
        g_share[i] *=pe_factor[i]
    pef_av = g_share.sum(axis = 1)
    gen["pef_av"] = pef_av
    
    
    
    
    if consHeat == "Yes":
        h = 1
    else:
        h = 0
    if consCool == "Yes":
        c = 1
    else:
        c = 0
        
    """
    METEO DATA
    """
    #Area for PV available and existing PV installations
    r_az = np.array([0,45,90,135,180,225,270,315])
    r_tilt = np.array([30])
    w_az = np.array([0,45,90,135,180,225,270,315])
    
    irrad_tilt_tilt_df = pd.DataFrame(columns = r_tilt, index = r_az, dtype = "O")
    irrad_tilt_wall_df = pd.DataFrame(columns = ["wall"], index = r_az, dtype = "O")
    
    tmy = pd.read_excel(in_timeseries, sheet_name = "MeteoData")
    irrad = tmy.drop(["temp", "windspeed"], axis = 1)
    
    irrad_tilt_flat = hp.irrad_tilt_v3(date_rng_hourly,irrad, loc, 180, loc["lat"], dt) / kW2X
    
    for bet in irrad_tilt_tilt_df.columns:
        for i in irrad_tilt_tilt_df.index:
            irrad_tilt_tilt_df.at[i,bet] = hp.irrad_tilt_v3(date_rng_hourly,irrad, loc, i, bet, dt) / kW2X
            
    if "pv_wall" in ders_df.index.values:
        for bet in irrad_tilt_wall_df.columns:
            for i in irrad_tilt_wall_df.index:
                irrad_tilt_wall_df.at[i,bet] = hp.irrad_tilt_v3(date_rng_hourly,irrad, loc, i, 90, dt) / kW2X
    
    irrad_tilt_flat_exist = {}
    for az in pv_ground_exist:
        irrad_tilt_flat_exist[az] = hp.irrad_tilt_v3(date_rng_hourly,irrad, loc, az, loc["lat"], dt) / kW2X
    
    tmy_dt = tmy.loc[tmy.index.repeat(60/dt)].reset_index(drop = True)
    temp = tmy_dt.loc[:,"temp"].values
    windspeed = tmy_dt.loc[:,"windspeed"].values

    COPc = (0.064 * 27 + 2.2726)-1
    
    
    """
    Electric Vehicle prep and chargers
    """
    
    if consEV == "Yes":
        evInput = pd.read_excel(in_ev, sheet_name = "EV_Input", index_col="Variable")
        v2g = evInput.at["v2g", "Value"] #Can EV battery used to cover demand?
        btm = evInput.at["btm", "Value"] #Grid interaction allowed?
        evSocIni = evInput.at["socIni", "Value"] #initial SOC of evs 
        chTarf = pd.read_excel(in_ev, sheet_name = "EV_Tariffs")
        chFeedIn = chTarf.loc[:,"Feed_in"].values
        chPublic = chTarf.loc[:,"publicCh"].values
        #creation of EVs
        evData = pd.read_excel(in_ev, sheet_name = "EV_creation", index_col="Name")
        elvehicles = {}
        for ind in evData.index:
            elvehicles[ind] = evClassKm.electricVehicle(evData.at[ind,"Capacity"],\
                evData.at[ind,"kW"],evData.at[ind,"avail"], evData.at[ind,"chAway"],\
                    evData.at[ind,"avkm"], evData.at[ind,"eff"], evData.at[ind,"Visiting"])
        
        EVS = []
        for ev in elvehicles:
            EVS.append(ev)
            
        avails = {} #availability of capacity
        ts_end = {} #ts before vehicle leaves
        ts_arrive = {} #ts of arrival
        ts_avail = {} #ts when vehicle is available for charging
        
        ev_socArr = {} #soc at time of arrival of ev
        ev_socMin = {} #soc that is required before leaving
        ev_pMax = {} #Max charging power
        ev_cons = {} #battery consumption when gone
        ev_chAway = {} #amount charged when gone
        for ev in EVS:
            end = []
            arr = []
            av = []
            avails[ev] = elvehicles[ev].avail()
            
            for ts in ts_i:
                if avails[ev][ts] != 0:
                    av.append(ts)
                if ts != len(ts_i)-1:
                    if avails[ev][ts] == 0 and avails[ev][ts+1] != 0:
                        arr.append(ts)
                    elif avails[ev][ts] !=0 and avails[ev][ts+1] == 0:
                        end.append(ts)
            ts_avail[ev] = av
            ts_arrive[ev] = arr
            ts_end[ev] = end
        
            ev_socMin[ev] = elvehicles[ev].socMincap()
            ev_pMax[ev] = elvehicles[ev].p_chMax
            ev_cons[ev] = elvehicles[ev].cons()
            ev_chAway[ev] = elvehicles[ev].chAway
            
            elOut = np.append(elOut, ev)
            elIn = np.append(elIn, ev)
        

        chData = pd.read_excel(in_ev, sheet_name = "Chargers", index_col="Variable")
        
        p_cs = {} #avail charging powers
        
        for i in chData.index:
            p_cs[i] = chData.at[i, "Value"]
        
        chData = chData.drop(["Value"], axis=1)
        chData['DER'] = chData.index
        ders_df = ders_df.append(chData)
        
    
    ders_df['CAPEX'] = ders_df['CAPEX'] * kW2X
    ders_df['OM_fix'] = ders_df['OM_fix'] * kW2X
    ders_df['OM_var'] = ders_df['OM_var'] * kW2X
    
    ders_df["an"] = (1 - (1/(1+i_rate)**ders_df["lifetime"]))/i_rate #Annuity factor
    ders_c = ders_df[ders_df['continously'] == 1].index.values
    ders_d = ders_df[ders_df['continously'] == 0].index.values
    
    #Creating optimization model
    m = ConcreteModel()
    '''
    Sets and Parameters
    '''
    
    '''
    Model Variable defintion and initialization with input data:
        
        c_c: Capacity of continuous technology
        c_d: Capacity of discrete technology
        power: 3 dimensional matrix with all power flows at each timestep
            elOut: Anything that can supply electricity
            elIn: Anything that can consume electricity
            ts_i: Timestep
        gridImport: Power imported from the grid at each ts
        gridExport: Power exported to the grid at each ts
        cost_I: Investement cost --> Overall CAPEX in year 0
        cost_OM: Operation and maintenance cost each year
        revenue: Revenue from electricty sales
    '''
    #Capacity of continous technology
    m.c_c = Var(ders_c, within = NonNegativeReals)
    #Capacity of discrete technology --> only Wind Turbines or chargers
    m.c_d = Var(ders_d, within = NonNegativeIntegers)
    
    #Power matrix --> From where to where does electricty flow at each ts
    m.power = Var(elOut, elIn, ts_i, within = NonNegativeReals, initialize = 0)
    #Equivalent heat matrix
    m.heat = Var(hOut, hIn, ts_i, within = NonNegativeReals, initialize = 0)
    #Equivalent cooling matrix
    m.cool = Var(cOut, cIn, ts_i, within = NonNegativeReals, initialize = 0)
    
    #Annual electrical energy output by each generation/conversion device
    m.elEner = Var(d_elOut, within = NonNegativeReals)
    #Equivalent heat output annually
    m.hEner = Var(hOut, within = NonNegativeReals)
    #Equivalent cooling output annually
    # m.cEner = Var(cOut, within = NonNegativeReals)    
    
    
    #Electricity grid import
    m.gridImport = Var(ts_i, within = NonNegativeReals)
    m.gridImportPE = Var(ts_i, within = NonNegativeReals) #PE for Primary Energy converted
    #Electricity grid export
    m.gridExport = Var(ts_i, within = NonNegativeReals)
    m.gridExportPE = Var(ts_i, within = NonNegativeReals) #PE for Primary Energy converted
    #Electricity self consumption rate
    m.selfCons = Var(ts_i, within = NonNegativeReals)
    m.selfCons_annual = Var(within = NonNegativeReals)
    # highest value of gridImport vector
    m.maxgridImport = Var(within = NonNegativeReals)
    
    # CO2 emissions of all emitting tec plus grid import
    m.emiss_CO2 = Var(emistec, ts_i, within = NonNegativeReals, initialize = 0)
    # CO2 emissions anually
    m.emiss_CO2_annual = Var(emistec, within = NonNegativeReals, initialize = 0)
    
    #if fueltec.size != 0:
    # Fuel Consumption in m3 of chosen fuel
    m.cons_fuel_vol = Var(fueltec, ts_i, within = NonNegativeReals, initialize = 0)
    # Fuel Consumption in kWh of chosen fuel
    m.cons_fuel_e = Var(fueltec,  ts_i, within = NonNegativeReals, initialize = 0)
    # Local particulate matter emissions
    m.emiss_PM = Var(fueltec,  ts_i, within = NonNegativeReals, initialize = 0)

    m.cost_I = Var(within = NonNegativeReals, initialize = 0)
    m.cost_OM_fix = Var(within = NonNegativeReals, initialize = 0)
    m.cost_OM_var = Var(within = NonNegativeReals, initialize = 0)
    m.cost_fuel = Var(within = NonNegativeReals, initialize = 0)
    m.cost_CO2 = Var(within = NonNegativeReals, initialize = 0)
    m.cost_Ext = Var(within = NonNegativeReals, initialize = 0)
    m.revenue = Var(within = NonNegativeReals, initialize = 0)
    
    '''
    Objective function:
        
        Annualized cost minimisation
    '''
    
    def obj_rule(m):
        return m.cost_I + m.cost_OM_fix + m.cost_OM_var + m.cost_fuel +\
            m.cost_CO2 + m.cost_Ext - m.revenue
    m.obj = Objective(rule = obj_rule)
    
    
    '''
    Economic Constraints:
        
        Investment: sum of the capacity of all considered technologies times 
            their CAPEX; differenciating among continuous and descrete technology
            divided by annuation factor
        
        Reveue: Indexed over year; sum of grid infeed for each der that has
            electricity as an output at each timestep; If the optimization is
            not accounting for the whole year, the value will be multiplied by
            the respective term in the end
            
        O&M fix costs: Indexed over years; sum of all operation and maintenance related
            annual fixed costs per kW, kWh or unit of installed capacity as well as
            contracted power cost per kW dia
        
        O&M variable costs: Indexed over year; Non-fuel cost per kWh of energy
            generation plus electricity imports according to tariff
        
        Fuel costs: Fuel consumption in kWh of HHV times cost per kWh plus monthly
            fix cost
        
        CO2 costs: Cost for each kWh of CO2 produced locally excluding grid
        
        External costs: TBD
    '''
    
    def investment_rule(m):
        return m.cost_I == sum(m.c_c[d_c] * ders_df.at[d_c, 'CAPEX'] /ders_df.at[d_c,"an"] for 
                               d_c in ders_c) + sum(m.c_d[d_d] * ders_df.at[d_d, 'CAPEX']\
                                / ders_df.at[d_d,"an"] for d_d in ders_d)
    m.const_inv = Constraint(rule = investment_rule)

        
    if consEV == "Yes":
        
        if type(tar_feedin) == np.ndarray:
            def revenue_rule(m):
                return m.revenue == sum(m.gridExport[ts] * dth * tar_feedin[ts] \
                    + sum(m.power[derOut,ev,ts] * elvehicles[ev].visiting * chFeedIn[ts] \
                          for derOut in elOut for ev in EVS) for ts in ts_i)
            m.const_revenue = Constraint(rule = revenue_rule)
        else:
            def revenue_rule(m):
                return m.revenue == sum(m.gridExport[ts] * dth * tar_feedin \
                    + sum(m.power[derOut,ev,ts] * elvehicles[ev].visiting * chFeedIn[ts] \
                          for derOut in elOut for ev in EVS) for ts in ts_i) 
            m.const_revenue = Constraint(rule = revenue_rule)
            
        if type(tariff) == np.ndarray:
            def omVar_rule(m):
                return m.cost_OM_var == sum(m.elEner[d_e] * ders_df.at[d_e,"OM_var"] for d_e in d_elOut) + \
                    sum(m.hEner[d_h] * ders_df.at[d_h,"OM_var"] for d_h in hOut) + \
                    sum(sum(m.power[ev,derIn,ts] * elvehicles[ev].visiting * tariff[ts] \
                          for derIn in elIn for ev in EVS) + m.gridImport[ts] * tariff[ts]\
                        for ts in ts_i) * dth 
            m.const_costOMVar = Constraint(rule = omVar_rule)
        else:
            def omVar_rule(m):
                return m.cost_OM_var == sum(m.elEner[d_e] * ders_df.at[d_e,"OM_var"] for d_e in d_elOut) + \
                    sum(m.hEner[d_h] * ders_df.at[d_h,"OM_var"] for d_h in hOut) + \
                    sum(sum(m.power[ev,derIn,ts] * elvehicles[ev].visiting * tariff \
                          for derIn in elIn for ev in EVS) + m.gridImport[ts] * tariff \
                        for ts in ts_i) * dth 
            m.const_costOMVar = Constraint(rule = omVar_rule)   

    else:
        if type(tar_feedin) == np.ndarray:
            def revenue_rule(m):
                return m.revenue == sum(m.gridExport[ts] * dth * tar_feedin[ts] \
                    for ts in ts_i)
            m.const_revenue = Constraint(rule = revenue_rule)
        else:
            def revenue_rule(m):
                return m.revenue == sum(m.gridExport[ts] * dth * tar_feedin \
                    for ts in ts_i) 
            m.const_revenue = Constraint(rule = revenue_rule)

        if type(tariff) == np.ndarray:
            def omVar_rule(m):
                return m.cost_OM_var == sum(m.elEner[d_e] * ders_df.at[d_e,"OM_var"] for d_e in d_elOut) + \
                    sum(m.hEner[d_h] * ders_df.at[d_h,"OM_var"] for d_h in hOut) + \
                    sum(m.gridImport[ts] * tariff[ts]\
                        for ts in ts_i) * dth 
            m.const_costOMVar = Constraint(rule = omVar_rule)
        else:
            def omVar_rule(m):
                return m.cost_OM_var == sum(m.elEner[d_e] * ders_df.at[d_e,"OM_var"] for d_e in d_elOut) + \
                    sum(m.hEner[d_h] * ders_df.at[d_h,"OM_var"] for d_h in hOut) + \
                    sum(m.gridImport[ts] * tariff \
                        for ts in ts_i) * dth 
            m.const_costOMVar = Constraint(rule = omVar_rule)   

    def omFix_rule(m):
        return m.cost_OM_fix == sum(m.c_c[d] * ders_df.at[d,"OM_fix"] for d in ders_c) + \
            sum(m.c_d[dd] * ders_df.at[dd,"OM_fix"] for dd in ders_d) + \
            m.maxgridImport * 1.2 * tar_kWdia * 365
    m.const_costOMFix = Constraint(rule = omFix_rule)     
    
    if fueltec.size != 0:
        def fuelCost_rule(m):
            return m.cost_fuel == (sum(m.cons_fuel_e[tec, ts] * cost_fuels["NG"]\
                for tec in ffueltec for ts in ts_i) + sum(m.cons_fuel_e[tec, ts] *\
                cost_fuels["BG"] for tec in bfueltec for ts in ts_i))
        m.const_costFuel = Constraint(rule = fuelCost_rule)                                                           
    
    #here emistec if grid should be accounted for or ffueltec/fueltec if grid should be excluded
    if CO2_cost == True:
        def CO2Cost_rule(m):
            return m.cost_CO2 == sum(m.emiss_CO2_annual[tec] for tec in emistec)\
                *0.001 * CO2_price 
        m.const_emissCost = Constraint(rule = CO2Cost_rule)

    '''
    PED Constraints:
        
        PED Balance Electricity: Only consideres electricity for the PED balance;
            heating and cooling if present will not be considered, but still need to be covered.
            Grid Exports - Grid Imports >= 0
            
        PED Balance Primary Energy (PEF): Takes into account all energy carriers;
            A constant primary energy factor for electricity is taken. 
            (Grid Export PE + fuel export) - (Grid Import PE + Fuel consumption) >= 0
            
        PED Balance Primary Energy Dynamic: Takes into account all energy carriers;
            A dynamic primary energy factor is used to account for the higher value of 
            renewably generated energy.
            (Grid Export + fuel export) - (Grid Import PEdyn + fuel cons) >= 0
            
        Self Consumption:
            Self consumption of each timestep determined by the power used by everything 
            that accepts electricity (Including the grid) and generated by generators
            (no storage) minus the generated power going to the grid
            
        Self Consumption annual:
            Self for annual measure
            
        Self Consumption minimum:
            Annual self consumption needs to be bigger than electric energy generated over
            year times the input variable (0-1), while 0 is no self consumption requirement
            and 1 means that all the energy generated needs to be consumed directly
            
        CO2 Annual rule:
            Sums up emissions of all emitting technology at each timestep
            
        CO2_cap rule:
            Emissions need to be smaller than cap annually
    '''
    if pedBalance == "Electricity":
        def pedBalance_ruleEL(m):
            return sum(m.gridExport[ts] for ts in ts_i) -\
                sum(m.gridImport[ts] for ts in ts_i) >= 0
        m.const_pedBalanceEL = Constraint(rule = pedBalance_ruleEL)

    if pedBalance == "PEF": #change export primary en fac! should not be the same
        def pedBalance_rulePEF(m):
            return (sum(m.gridExport[ts] for ts in ts_i) * primaryEn_fac_el * dth + fuel_export -\
                sum(m.gridImport[ts] for ts in ts_i) * primaryEn_fac_el * dth -\
                sum(m.cons_fuel_e[tec, ts] for tec in ffueltec for ts in ts_i))  >= 0
        m.const_pedBalancePEF = Constraint(rule = pedBalance_rulePEF)
        
    if pedBalance == "PEF_dynamic":
        def pedBalance_rulePEF_dyn(m):
            return (sum(m.gridExportPE[ts] for ts in ts_i) *dth + fuel_export -\
                sum(m.gridImportPE[ts] for ts in ts_i) *dth -\
                sum(m.cons_fuel_e[tec, ts] for tec in ffueltec for ts in ts_i))  >= 0
        m.const_pedBalancePEF_dyn = Constraint(rule = pedBalance_rulePEF_dyn)
    
    
    # annual CO2 emission in kg per technology
    if emistec.size != 0:
        def CO2annual_rule(m, tec):
            return m.emiss_CO2_annual[tec] == sum(m.emiss_CO2[tec, ts] \
                    for ts in ts_i)
        m.const_CO2annual = Constraint(emistec, rule = CO2annual_rule)
    
    if CO2_cap == True and fueltec.size != 0:
        def CO2cap_rule(m):
            return sum(m.emiss_CO2_annual[tec] for tec in fueltec) <=\
                CO2_cap_value 
        m.const_CO2cap = Constraint(rule = CO2cap_rule) 
        
        
    '''
    Grid Constraints:
        grid2grid: at each timestep electricity cannot flow from the grid to 
            the grid
            
        import_rule: 
            at each timestep the sum over all electricty that is supplied
            by the grid to consumers/storage is the grid import
        
        importPE_rule: 
            at each ts the sum over all electricity supplied by the grid
            times the pef of the current grid mix
            
        export_rule: 
            at each timestep the sum of electricty supplied to the grid
            by generators/storage is the export
            
        exportPE_rule:
            at each timestep the sum of electricity export is counted to be worth
            the pef of the most expenisve tech at the merit order curve
        
        gridemiss_rule:
            Sum of grid related CO2 emissions over timesteps   
            
        maxgridImport_rule:
            Finds maximum grid power imported
    '''
    
    def grid2grid_rule(m, ts):
        return m.power['grid','grid', ts] == 0
    m.const_grid2grid = Constraint(ts_i, rule = grid2grid_rule)
    
    def import_rule(m, ts): #maybe OUT
        return m.gridImport[ts] == sum(m.power['grid',derIn,ts] for derIn in elIn)
    m.const_import = Constraint(ts_i, rule = import_rule)    
    
    if pedBalance == "PEF_dynamic":
        def importPE_rule(m, ts):
            if gen.at[ts,merit_order[max(merit_order)]] != 0:
                return m.gridImportPE[ts] == sum(m.power['grid',derIn,ts] for derIn in elIn[elIn != "load_e"])\
                    * pe_factor[merit_order[max(merit_order)]] + m.power['grid','load_e',ts]\
                        * gen.at[ts, "pef_av"]
            elif gen.at[ts,merit_order[max(merit_order)-1]] != 0:
                return m.gridImportPE[ts] == sum(m.power['grid',derIn,ts] for derIn in elIn[elIn != "load_e"])\
                    * pe_factor[merit_order[max(merit_order)-1]] + m.power['grid','load_e',ts]\
                        * gen.at[ts, "pef_av"]
            else:
                return m.gridImportPE[ts] == sum(m.power['grid',derIn,ts] for derIn in elIn[elIn != "load_e"])\
                    * pe_factor[merit_order[max(merit_order)-2]] + m.power['grid','load_e',ts]\
                        * gen.at[ts, "pef_av"]
        m.const_importPE = Constraint(ts_i, rule = importPE_rule)
    
    def export_rule(m, ts): #maybe OUT
        return m.gridExport[ts] == sum(m.power[derOut,'grid',ts] for derOut in elOut) 
    m.const_export = Constraint(ts_i, rule = export_rule)
    
    if pedBalance == "PEF_dynamic":
        def exportPE_rule(m, ts):
            if gen.at[ts,merit_order[max(merit_order)]] != 0:
                return m.gridExportPE[ts] == m.gridExport[ts] * pe_factor[merit_order[max(merit_order)]]
            elif gen.at[ts,merit_order[max(merit_order)-1]] != 0:
                return m.gridExportPE[ts] == m.gridExport[ts] * pe_factor[merit_order[max(merit_order)-1]]
            else:
                return m.gridExportPE[ts] == m.gridExport[ts] * pe_factor[merit_order[max(merit_order)-2]]
            
        m.const_exportPE = Constraint(ts_i, rule = exportPE_rule)
    
    def gridemiss_rule(m, ts):
        return m.emiss_CO2["grid", ts] == m.gridImport[ts] * (emiss_fac["grid"])
    m.const_gridEmiss = Constraint(ts_i, rule = gridemiss_rule)
    
    def maxgridImport_rule(m, ts):
        return m.maxgridImport >= m.gridImport[ts]
    m.const_maxgridImport = Constraint(ts_i, rule = maxgridImport_rule)
    
    # if grid_limit == True:
    #     def gridImportLimit_rule(m,ts):
    #         return m.gridImport[ts] <= grid_limitValue
    #     m.const_gridImportLimit = Constraint(ts_i, rule = gridImportLimit_rule)
        
    #     def gridExportLimit_rule(m,ts):
    #         return m.gridExport[ts] <= grid_limitValue
    #     m.const_gridExportLimit = Constraint(ts_i, rule = gridExportLimit_rule)
        
    if consEV == "Yes":
        if btm == True:
            # Behind The Meter --> No grid interaction with PV or Bat or evs
            def btm_rule(m, ts):
                return m.power["pv","grid", ts] + m.power["bat","grid", ts] + \
                    m.power["grid", "bat", ts] + sum( m.power["grid", evs, ts] +\
                    m.power[evs, "grid", ts] for evs in EVS) == 0
            m.const_btm = Constraint(ts_i, rule = btm_rule)
        
    '''
    Annual Energy Constraints:
        
        elEner_rule:
            Sum of all power output over each timestep for any device that can
            output electricity
        
        hEner_rule:
            Equivalent for heat
            
        cEner_rule:
            Sames
    '''
    def elEner_rule(m, elSup):
        return m.elEner[elSup] == sum(m.power[elSup, elCons, ts] \
            for elCons in elIn for ts in ts_i) * dth
    m.const_annualelEner = Constraint(d_elOut, rule = elEner_rule)
    
    if h == 1:
        def hEner_rule(m, hSup):
            return m.hEner[hSup] == sum(m.heat[hSup, hCons, ts] \
                for hCons in hIn for ts in ts_i) * dth
        m.const_annualhEner = Constraint(hOut, rule = hEner_rule)

    # if c == 1:
    #     def cEner_rule(m, cSup):
    #         return m.cEner[cSup] == sum(m.cool[cSup, cCons, ts] \
    #             for cCons in cIn for ts in ts_i) * dth
    #     m.const_annualcEner = Constraint(cOut, rule = cEner_rule)
    
    '''
    Load Constraints:
        Electricity: Sum of power by each available electricity DER/grid at each
            timestep towards the load needs to meet the load at each timestep
        
        Heat: Sum of heat by each available heat DER at each
            timestep towards the load needs to meet the load at each timestep
        
        Cooling: Sames
    '''
    
    def eloadBalance_rule(m, ts):
        return sum(m.power[d,'load_e',ts] for d in elOut) == load_e[ts]
    m.const_eloadBalance = Constraint(ts_i, rule = eloadBalance_rule)

    if h == 1:
        def hloadBalance_rule(m, ts):
            return sum(m.heat[d,'load_h',ts] for d in hOut) == load_h[ts]
        m.const_hloadBalance = Constraint(ts_i, rule = hloadBalance_rule)
        
    if c == 1:
        def cloadBalance_rule(m, ts):
            return sum(m.cool[d,'load_c',ts] for d in cOut) == load_c[ts] 
        m.const_cloadBalance = Constraint(ts_i, rule = cloadBalance_rule)
    
    '''
    PV Constraints:
        
        pv_ground constraints encompass installations on the ground or flat 
        roofs
        
        pv_roof constraints encompass installations on all titled roofs
        
        pv_wall constraints encompass installations on walls with 90ª tilt
    '''
    
    #PV covered roof area at different angles
    m.pv_roofPanelArea = Var(r_az, r_tilt, within = NonNegativeReals, initialize = 0)
    m.pv_wallPanelArea = Var(w_az, within = NonNegativeReals, initialize = 0)
    m.pv_groundPanelArea = Var(within = NonNegativeReals, initialize = 0)

    
    if "pv_ground" in ders_df.index.values:
        
        # Relation between installed capacity and panel area        
        def pv_groundPanelAreaRule(m):
            return m.c_c["pv_ground"] == m.pv_groundPanelArea * pv_eff
        m.const_pvGroundArea = Constraint(rule = pv_groundPanelAreaRule) 
        
        # Total amount of PV power generated by ground generation at each timestep --> needs to be fastened!!!!
        def pv_groundGenRule(m, ts):
            return sum(m.power["pv_ground", elCons, ts] for elCons in elIn) == (m.pv_groundPanelArea * irrad_tilt_flat[ts])* pv_eff * pv_pr +\
                sum(pv_ground_exist[az] * irrad_tilt_flat_exist[az][ts]\
                    for az in pv_ground_exist) * pv_eff * pv_pr
        m.const_pvGroundGen = Constraint(ts_i, rule = pv_groundGenRule)
        
        
        #Limits the area available for ground and flat roof pv
        
        if "solar_thermal_ground" not in ders_df.index.values: 
            def panel_groundAreaLimit(m):
                return m.pv_groundPanelArea <= pv_groundareaAvail * pv_gcr_ground
            m.const_panelGroundAreaLimit = Constraint(rule = panel_groundAreaLimit)
        
    if "pv_roof" in ders_df.index.values:
        
        # Relation between installed capacity and panel area
        def pv_roofPanelAreaRule(m):
            return m.c_c["pv_roof"] == sum(m.pv_roofPanelArea[az, tilt] for az in r_az \
                for tilt in r_tilt) * pv_eff
        m.const_pvRoofArea = Constraint(rule = pv_roofPanelAreaRule)
        
        # Total amount of tilted roof pv generation at each timestep --> needs to be fastened!!!
        def pv_roofGenRule(m, ts):
            return sum(m.power["pv_roof", elCons, ts] for elCons in elIn) == sum((m.pv_roofPanelArea[az, tilt] + \
                pv_roof_exist.at[az,tilt]) * irrad_tilt_tilt_df.at[az,tilt][ts]\
                for az in r_az for tilt in r_tilt) * pv_eff * pv_pr
        m.const_pvRoofGen = Constraint(ts_i, rule = pv_roofGenRule)
    
        
        #Limits the area available for ground and flat roof pv
        if "solar_thermal_roof" not in ders_df.index.values:
            def panel_roofAreaLimit(m, az, tilt):
                return m.pv_roofPanelArea[az, tilt] <= pv_roofwallareaAvail.at[az, tilt] * pv_gcr_roof
            m.const_panelRoofAreaLimit = Constraint(r_az, r_tilt, rule = panel_roofAreaLimit)
        
    if "pv_wall" in ders_df.index.values: #Needs to be fixed to fit to new input
        
        # Relation between installed capacity and panel area
        def pv_wallPanelAreaRule(m):
            return m.c_c["pv_wall"] == sum(m.pv_wallPanelArea[az] for az in w_az)\
                * pv_eff
        m.const_pvWallArea = Constraint(rule = pv_wallPanelAreaRule)
        
        # Total amount of 90º wall pv generation at each timestep ---> needs to be fastened!!!
        def pv_wallGenRule(m, ts):
            return m.pv_wallGen[ts] == (sum(m.pv_wallPanelArea[az] * \
                hp.irrad_tilt(date_rng[ts], irrad.iloc[ts,:], loc, az, 90)\
                for az in w_az) + sum(pv_wall_exist[az] * \
                hp.irrad_tilt(date_rng[ts], irrad.iloc[ts,:], loc, az, 90)\
                for az in pv_wall_exist)) * pv_eff * pv_pr
        m.const_pvWallGen = Constraint(ts_i, rule = pv_wallGenRule)
    
        # Total wall generated PV power needs to be equal as the sum of 
        # wall generated PV power to all consumers
        def pv_wallGen2X(m, ts):
            return m.pv_wallGen[ts] == \
                sum(m.power["pv_wall", elCons, ts] for elCons in elIn)
        m.const_pvWallGen2X = Constraint(ts_i, rule = pv_wallGen2X)
        
        #Limits the area available for wall mounted pv
        def pv_wallAreaLimit(m, az):
            return m.pv_wallPanelArea[az] <= pv_roofwallareaAvail.at[az, 90] * pv_gcr_wall
        m.const_pvWallAreaLimit = Constraint(w_az, rule = pv_wallAreaLimit)
    
    
    '''
    Battery Constraints
    '''
    if "battery" in ders_df.index.values:
        
        #battery state of charge in kWh
        m.batt_soc = Var(ts_i, within = NonNegativeReals, initialize = 0)
        
        def batt2batt_rule(m, ts):
            return m.power["battery", "battery", ts] == 0
        m.const_batt2batt = Constraint(ts_i, rule = batt2batt_rule)
    
          
        # SOC maximal the installed capacity of batteries
        def battSOCmax_rule(m, ts):
            return m.batt_soc[ts] <= m.c_c["battery"]
        m.const_battSOCmax = Constraint(ts_i, rule = battSOCmax_rule)
        
        
        # SOC update each ts
        def battSOCupdate_rule(m, ts):
            if ts == 0:
                return m.batt_soc[ts] == batt_socIni + \
                    (sum(m.power[elSup, "battery", ts] for elSup in elOut) * \
                     sqrt(batt_eff) - sum(m.power["battery", elCons, ts] for elCons in elIn) /\
                        sqrt(batt_eff)) * dth
            else:
                return m.batt_soc[ts] == m.batt_soc[ts-1] +\
                    (sum(m.power[elSup, "battery", ts] for elSup in elOut) * \
                     sqrt(batt_eff) - sum(m.power["battery", elCons, ts] for elCons in elIn) /\
                        sqrt(batt_eff)) * dth
        m.const_battSOCupdate = Constraint(ts_i, rule = battSOCupdate_rule)
        
        # Discharge power limited at each timestep
        def battpowerDisc_rule(m, ts):
            return sum(m.power["battery", elCons, ts] for elCons in elIn) / sqrt(batt_eff) <= \
                m.c_c["battery"] * batt_capa2pow
        m.const_battPowDisc = Constraint(ts_i, rule = battpowerDisc_rule)
        
        #Charging power limted at each timestep
        def battpowerCh_rule(m, ts):
            return sum(m.power[elSup, "battery", ts] for elSup in elOut) * sqrt(batt_eff) <= \
                m.c_c["battery"] * batt_capa2pow
        m.const_battPowCh = Constraint(ts_i, rule = battpowerCh_rule)
        
        if batt_grid2batt == 0:
            def grid2batt_rule(m, ts):
                return m.power["grid", "battery", ts] == 0
            m.const_grid2batt = Constraint(ts_i, rule = grid2batt_rule)
     
    """
    Hot Water Storage
    """
    if "thermal_storage" in ders_df.index.values:
        
        #storage state of charge in kWh
        m.tss_soc = Var(ts_i, within = NonNegativeReals, initialize = 0)
        
        def tss2tss_rule(m, ts):
            return m.heat["thermal_storage", "thermal_storage", ts] == 0
        m.const_tss2tss = Constraint(ts_i, rule = tss2tss_rule)
    
          
        # SOC maximal the installed capacity of tss
        def tssSOCmax_rule(m, ts):
            return m.tss_soc[ts] <= m.c_c["thermal_storage"]
        m.const_tssSOCmax = Constraint(ts_i, rule = tssSOCmax_rule)
        
        
        # SOC update each ts
        def tssSOCupdate_rule(m, ts):
            if ts == 0:
                return m.tss_soc[ts] == tss_socIni * tss_self + \
                    (sum(m.heat[hSup, "thermal_storage", ts] for hSup in hOut) * \
                     sqrt(tss_eff) - sum(m.heat["thermal_storage", hCons, ts] for hCons in hIn) /\
                        sqrt(tss_eff)) * dth
            else:
                return m.tss_soc[ts] == m.tss_soc[ts-1] * tss_self +\
                    (sum(m.heat[hSup, "thermal_storage", ts] for hSup in hOut) * \
                     sqrt(tss_eff) - sum(m.heat["thermal_storage", hCons, ts] for hCons in hIn) /\
                        sqrt(tss_eff)) * dth
        m.const_tssSOCupdate = Constraint(ts_i, rule = tssSOCupdate_rule)
        
        # Discharge power limited at each timestep
        def tsspowerDisc_rule(m, ts):
            return sum(m.heat["thermal_storage", hCons, ts] for hCons in hIn) / sqrt(tss_eff) <= \
                m.c_c["thermal_storage"] * tss_capa2pow
        m.const_tssPowDisc = Constraint(ts_i, rule = tsspowerDisc_rule)
        
        #Charging power limted at each timestep
        def tsspowerCh_rule(m, ts):
            return sum(m.heat[hSup, "thermal_storage", ts] for hSup in hOut) * sqrt(tss_eff) <= \
                m.c_c["thermal_storage"] * tss_capa2pow
        m.const_tssPowCh = Constraint(ts_i, rule = tsspowerCh_rule)
        
    """
    EVs
    """
    if consEV == "Yes":
        
        m.soc_ev = Var(EVS, ts_i, within = NonNegativeReals, initialize = 0) #SOC of evs
        
        #EVs cannot send electricity to itself or among each other
        def ev2ev_rule(m, ts):
            return sum(m.power[evsout, evsin, ts] for evsout in EVS for evsin in EVS) == 0
        m.const_ev2ev = Constraint(ts_i, rule = ev2ev_rule)
        
        #if V2G is true the chargers work bidirectionally
        if v2g == False:
            def v2g_rule(m,ts):
                return sum(m.power[evs, elCons,ts] for evs in EVS for elCons in elIn) == 0
            m.const_v2g = Constraint(ts_i, rule = v2g_rule)
        
        # SOC needs to be smaller than the battery capa of each car
        def evSOCmax_rule(m, ev, ts):
            return m.soc_ev[ev, ts] <= elvehicles[ev].capa
        m.const_evSOXmax = Constraint(EVS, ts_i, rule = evSOCmax_rule)
                
        # SOC update each ts
        def evSOCupdate_rule(m, ev, ts):
            if ts == 0:
                return m.soc_ev[ev, ts] == evSocIni*elvehicles[ev].capa + \
                    (sum(m.power[elSup, ev, ts] for elSup in elOut) * \
                     sqrt(batt_eff) - sum(m.power[ev, elCons, ts] for elCons in elIn) /\
                        sqrt(batt_eff)) * dth
            elif ts in ts_arrive[ev]:
                return m.soc_ev[ev, ts] == m.soc_ev[ev,ts-1] - ev_cons[ev] + ev_chAway[ev]
            elif ts in ts_avail[ev]:
                return m.soc_ev[ev, ts] == m.soc_ev[ev, ts-1] +\
                    (sum(m.power[elSup, ev, ts] for elSup in elOut) * \
                     sqrt(batt_eff) - sum(m.power[ev, elCons, ts] for elCons in elIn) /\
                        sqrt(batt_eff)) * dth
            else:
                return m.soc_ev[ev,ts] == m.soc_ev[ev,ts-1]
            
        m.const_evSOCupdate = Constraint(EVS, ts_i, rule = evSOCupdate_rule)
    
        # minimum SOC at end of stationary period
        
        def evSOCmin_rule(m, ev, ts):
            if ts in ts_end[ev]:
                return m.soc_ev[ev, ts] >= ev_socMin[ev]
            else:
                return m.soc_ev[ev, ts] >= 0
        m.const_evSOCmin = Constraint(EVS, ts_i, rule = evSOCmin_rule)
        
        
        # Discharge power limited at each timestep
        def evpowerDisc_rule(m, ev, ts):
            if ts in ts_avail[ev]:
                return sum(m.power[ev, elCons, ts] for elCons in elIn) / sqrt(batt_eff) <= ev_pMax[ev]
            else:
                return sum(m.power[ev, elCons, ts] for elCons in elIn) == 0
        m.const_evPowDisc = Constraint(EVS, ts_i, rule = evpowerDisc_rule)
        
        #Charging power limted at each timestep
        def evpowerCh_rule(m, ev, ts):
            if ts in ts_avail[ev]:
                return sum(m.power[elSup, ev, ts] for elSup in elOut) * sqrt(batt_eff) <= ev_pMax[ev]
            else:
                return sum(m.power[elSup, ev, ts] for elSup in elOut) == 0
        m.const_evPowCh = Constraint(EVS, ts_i, rule = evpowerCh_rule)
        
        
        """
        Chargers
        """
        
        m.b_cs1 = Var(EVS, within = Binary, initialize = 0) #binary variable for cs1 for each ev
        m.b_cs2 = Var(EVS, within = Binary, initialize = 0) #binary variable for cs2 for each ev
        m.b_cs3 = Var(EVS, within = Binary, initialize = 0) #binary variable for cs3 for each ev
        
        #Defining the charging power allowed by the charger
        def chargingPower_rule(m, ev, ts):
            return sum(m.power[elSup, ev, ts] for elSup in elOut) <= p_cs["cs1"] * m.b_cs1[ev] +\
                p_cs["cs2"] * m.b_cs2[ev] + p_cs["cs3"] * m.b_cs3[ev]
        m.const_chargingPower = Constraint(EVS, ts_i, rule = chargingPower_rule)
        
        #only one type of charger can be used for one vehicle
        def chargerBin_rule(m,ev):
            return m.b_cs1[ev] + m.b_cs2[ev] + m.b_cs3[ev] <=1
        m.const_bincharg = Constraint(EVS,rule = chargerBin_rule)
        
        #Amount of CS1 chargers
        def nrCS1_rule(m):
            return m.c_d["cs1"] == sum(m.b_cs1[ev] for ev in EVS)
        m.const_nrCS1 = Constraint(rule = nrCS1_rule)
        
        #Amount of CS2 chargers
        def nrCS2_rule(m):
            return m.c_d["cs2"] == sum(m.b_cs2[ev] for ev in EVS)
        m.const_nrCS2 = Constraint(rule = nrCS2_rule)
        
        #Amount of CS3 chargers
        def nrCS3_rule(m):
            return m.c_d["cs3"] == sum(m.b_cs3[ev] for ev in EVS)
        m.const_nrCS3 = Constraint(rule = nrCS3_rule)

    
    '''
    Air sourced Heat Pumps
    '''
    
    m.b_asgshp = Var(within = Binary) #binary variable of GSHP (0) vs ASHP (1)

    
    if "heat_pump_air" in ders_df.index.values:
        
        if c==0:
        
            def ashpHeat_rule(m,ts):
                if temp[ts] <= -20:
                    return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) == \
                        sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut)
                elif temp[ts] >= 35:
                    return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) == \
                        sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut) * 4.5
                else:
                    return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) == \
                        sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut)*\
                            (0.064 * temp[ts] + 2.2726)
            m.const_ashpHeat = Constraint(ts_i, rule = ashpHeat_rule)
            
            def ashpCapa_rule(m,ts):
                return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) <=\
                    m.c_c["heat_pump_air"]
            m.const_ashpCapa = Constraint(ts_i, rule = ashpCapa_rule)
            
            def asvsgs_rule(m):
                return m.c_c["heat_pump_air"] <= 9999 * m.b_asgshp
            m.constASvsGS = Constraint(rule = asvsgs_rule)

        elif c==1:
            #Binary variables for heating or cooling mode of heat pump 1 for heating 0 for cooling
            m.b_hph = Var(ts_i, within=Binary)
            
            
            def ashpcons_rule(m,ts):
                if temp[ts] <= -20:
                    return sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut) == \
                        sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) +\
                        sum(m.cool["heat_pump_air", cCons, ts] for cCons in cIn)/COPc
                elif temp[ts] >= 35:
                    return sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut) == \
                        sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn)/4.5 +\
                        sum(m.cool["heat_pump_air", cCons, ts] for cCons in cIn)/COPc
                else:
                    return sum(m.power[elSup, "heat_pump_air",ts] for elSup in elOut) == \
                        sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn)/\
                        (0.064 * temp[ts] + 2.2726) + sum(m.cool["heat_pump_air", cCons, ts]\
                        for cCons in cIn)/COPc
            m.const_ashpcons = Constraint(ts_i, rule = ashpcons_rule)
            

            def ashpheatdecision_rule(m,ts):
                return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) <=\
                    9999 * m.b_hph[ts] + 0 * (1-m.b_hph[ts])
            m.const_ashpheatdecision = Constraint(ts_i, rule = ashpheatdecision_rule)
            
            def ashpcooldecision_rule(m,ts):
                return sum(m.cool["heat_pump_air", cCons, ts] for cCons in cIn) <=\
                    0 * m.b_hph[ts] + 9999 * (1-m.b_hph[ts])
            m.const_ashpcooldecision = Constraint(ts_i, rule = ashpcooldecision_rule)

            def ashpCapah_rule(m,ts):
                return sum(m.heat["heat_pump_air", hCons, ts] for hCons in hIn) <=\
                    m.c_c["heat_pump_air"] 
            m.const_ashphCapa = Constraint(ts_i, rule = ashpCapah_rule)

            def ashpCapac_rule(m,ts):
                return sum(m.cool["heat_pump_air", cCons, ts] for cCons in cIn) <=\
                    m.c_c["heat_pump_air"] * (COPc/COPc+1) 
            m.const_ashpcCapa = Constraint(ts_i, rule = ashpCapac_rule)
            
            def asvsgs_rule(m):
                return m.c_c["heat_pump_air"] <= 9999 * m.b_asgshp
            m.constASvsGS = Constraint(rule = asvsgs_rule)
            
            
    '''
    Ground Sourced Heat Pump 
    '''
    if "heat_pump_ground" in ders_df.index.values:
        
        
        if c==0:
        
            def gshpHeat_rule(m,ts):
                return sum(m.heat["heat_pump_ground", hCons, ts] for hCons in hIn) == \
                    sum(m.power[elSup, "heat_pump_ground",ts] for elSup in elOut) * gshpCOP_h
            m.const_gshpHeat = Constraint(ts_i, rule = gshpHeat_rule)
            
            def gshpCapa_rule(m,ts):
                return sum(m.heat["heat_pump_ground", hCons, ts] for hCons in hIn) <=\
                    m.c_c["heat_pump_ground"]
            m.const_gshpCapa = Constraint(ts_i, rule = gshpCapa_rule)

            def asvsgs_rule2(m):
                return m.c_c["heat_pump_ground"] <= 9999 * (1-m.b_asgshp)
            m.constASvsGS2 = Constraint(rule = asvsgs_rule2)

        elif c==1:
            #Binary variables for heating or cooling mode of heat pump 1 for heating 0 for cooling
            m.b_gshph = Var(ts_i, within=Binary)
            
            
            def gshpcons_rule(m,ts):
                return sum(m.power[elSup, "heat_pump_ground",ts] for elSup in elOut) == \
                    sum(m.heat["heat_pump_ground", hCons, ts] for hCons in hIn)/gshpCOP_h +\
                    sum(m.cool["heat_pump_ground", cCons, ts] for cCons in cIn)/gshpCOP_c
            m.const_gshpcons = Constraint(ts_i, rule = gshpcons_rule)
            

            def gshpheatdecision_rule(m,ts):
                return sum(m.heat["heat_pump_ground", hCons, ts] for hCons in hIn) <=\
                    9999 * m.b_gshph[ts] + 0 * (1-m.b_gshph[ts])
            m.const_gshpheatdecision = Constraint(ts_i, rule = gshpheatdecision_rule)
            
            def gshpcooldecision_rule(m,ts):
                return sum(m.cool["heat_pump_ground", cCons, ts] for cCons in cIn) <=\
                    0 * m.b_gshph[ts] + 9999 * (1-m.b_gshph[ts])
            m.const_gshpcooldecision = Constraint(ts_i, rule = gshpcooldecision_rule)

            def gshpCapah_rule(m,ts):
                return sum(m.heat["heat_pump_ground", hCons, ts] for hCons in hIn) <=\
                    m.c_c["heat_pump_ground"]
            m.const_gshphCapa = Constraint(ts_i, rule = gshpCapah_rule)

            def gshpCapac_rule(m,ts):
                return sum(m.cool["heat_pump_ground", cCons, ts] for cCons in cIn) <=\
                    m.c_c["heat_pump_ground"] * (gshpCOP_c/gshpCOP_c+1)
            m.const_gshpcCapa = Constraint(ts_i, rule = gshpCapac_rule)
            
            def asvsgs_rule2(m):
                return m.c_c["heat_pump_ground"] <= 9999 * (1-m.b_asgshp)
            m.constASvsGS2 = Constraint(rule = asvsgs_rule2)

        
    '''
    Electric Boiler 
    '''
    if "boiler_electric" in ders_df.index.values:
        
        # Heat generation by electric boiler from electricity
        def electricBoilerHeat_rule(m, ts):
            return sum(m.heat["boiler_electric", hCons, ts] for hCons in hIn) ==\
               sum(m.power[elSup, "boiler_electric", ts] for elSup in elOut) *\
                   elBoiler_eff
        m.const_electricBoilerHeat = Constraint(ts_i, rule = electricBoilerHeat_rule)
        
        # Max heat generated can only be as high as rated capacity
        def electricBoilerCapa_rule(m, ts):
            return sum(m.heat["boiler_electric", hCons, ts] for hCons in hIn) <=\
                m.c_c["boiler_electric"]
        m.const_electricBoilerCapa = Constraint(ts_i, rule = electricBoilerCapa_rule)            
            
            
    '''
    Micro Combined Heat and Power 
    '''
    if "chp_micro" in ders_df.index.values:
        # Power Ouput can be maximal the rated electric capacity
        def capa_rule(m, ts):
            return sum(m.power["chp_micro", eCons, ts] for eCons in elIn) <= \
                m.c_c["chp_micro"]
        m.const_Capa_mCHP = Constraint(ts_i, rule = capa_rule)
        # Heat Output per power output --> constant assumption
        def p2h_rule(m, ts):
            return sum(m.heat["chp_micro", hCons, ts] for hCons in hIn) ==\
                sum(m.power["chp_micro", eCons, ts] for eCons in elIn) * \
                    mchp_p2h
        m.const_p2h_mCHP = Constraint(ts_i, rule = p2h_rule)
        # No generation below 50% capacity allowed --> above constant efficiency assumed
        def minProduction_rule(m, ts):
            return sum(m.power["chp_micro", eCons, ts] for eCons in elIn) >= \
                0.5 * m.c_c["chp_micro"]
        m.const_minProd_mCHP = Constraint(ts_i, rule = minProduction_rule)
        
        # Primary energy consumption of fuel in kWh HHV
        def EConsmCHP_rule(m, ts):
            return m.cons_fuel_e["chp_micro", ts] == (sum(m.power["chp_micro", eCons, ts] for eCons in elIn) /\
                mchp_eff_e)
        m.const_EConsmCHP = Constraint(ts_i, rule = EConsmCHP_rule)

        # Fuel consumption in m3 of NG
        def VolConsmCHP_rule(m, ts):
            return m.cons_fuel_vol["chp_micro", ts] == m.cons_fuel_e["chp_micro", ts] / hhv["NG"]
        m.const_VolConsmCHP = Constraint(ts_i, rule = VolConsmCHP_rule)

        def chpEmiss_rule(m, ts):
            return m.emiss_CO2["chp_micro", ts] == m.cons_fuel_e["chp_micro", ts] * emiss_fac["NG"]
        m.const_chpEmiss = Constraint(ts_i, rule = chpEmiss_rule)

    '''
    Micro Combined Heat and Power (Bio)
    '''
    if "chp_micro_b" in ders_df.index.values:
        # Power Ouput can be maximal the rated electric capacity
        def capa_ruleB(m, ts):
            return sum(m.power["chp_micro_b", eCons, ts] for eCons in elIn) <= \
                m.c_c["chp_micro_b"]
        m.const_Capa_mCHPb = Constraint(ts_i, rule = capa_ruleB)
        # Heat Output per power output --> constant assumption
        def p2h_ruleB(m, ts):
            return sum(m.heat["chp_micro_b", hCons, ts] for hCons in hIn) ==\
                sum(m.power["chp_micro_b", eCons, ts] for eCons in elIn) * \
                    mchp_p2h
        m.const_p2h_mCHPb = Constraint(ts_i, rule = p2h_ruleB)
        # No generation below 50% capacity allowed --> above constant efficiency assumed
        def minProduction_ruleB(m, ts):
            return sum(m.power["chp_micro_b", eCons, ts] for eCons in elIn) >= \
                0.5 * m.c_c["chp_micro_b"]
        m.const_minProd_mCHPb = Constraint(ts_i, rule = minProduction_ruleB)
        # Primary energy consumption of fuel in kWh HHV
        def EConsmCHP_ruleB(m, ts):
            return m.cons_fuel_e["chp_micro_b", ts] == (sum(m.power["chp_micro_b", eCons, ts] for eCons in elIn) /\
                mchp_eff_e)
        m.const_EConsmCHPb = Constraint(ts_i, rule = EConsmCHP_ruleB)

        # Fuel consumption in m3 of BG        
        def VolConsmCHP_ruleB(m, ts):
            return m.cons_fuel_vol["chp_micro_b", ts] == m.cons_fuel_e["chp_micro_b", ts] / hhv["BG"]
        m.const_VolConsmCHPb = Constraint(ts_i, rule = VolConsmCHP_ruleB)

        def chpbEmiss_rule(m, ts):
            return m.emiss_CO2["chp_micro_b", ts] == m.cons_fuel_e["chp_micro_b", ts] * emiss_fac["BG"]
        m.const_chpbEmiss = Constraint(ts_i, rule = chpbEmiss_rule)
        
    

    '''
    Wind Turbine:
        
        windPower_rule:
            power output of WT at each ts depending on windspeed needs to be
            equal to the sum of power going to consumers
            
        Area restrictions: TBD
        
        Include m.c_d["wind"] --> how many wind turbines are there --> TBD
    '''
    
    if "wind" in ders_df.index.values:
        
        def windPower_rule(m, ts):
            if windspeed[ts] < wt_cutInWS:
                return sum(m.power["wind", elCons, ts] for elCons in elIn) == 0
            elif wt_cutInWS <= windspeed[ts] < wt_ratedWS:
                return sum(m.power["wind", elCons, ts] for elCons in elIn) ==\
                    windspeed[ts] * wt_m + wt_t 
            elif wt_ratedWS <= windspeed[ts] < wt_cutOutWS:
                return sum(m.power["wind", elCons, ts] for elCons in elIn) ==\
                    wt_ratedP
            else:
                return sum(m.power["wind", elCons, ts] for elCons in elIn) == 0
            
        m.const_windPower = Constraint(ts_i, rule = windPower_rule)

    '''
    Gas Boiler 
    '''
    
    if "boiler" in ders_df.index.values:
        
        def BoilerHeat_rule(m, ts):
            return sum(m.heat["boiler", hCons, ts] for hCons in hIn) <= \
                m.c_c["boiler"]
        m.const_BoilerCapa = Constraint(ts_i, rule = BoilerHeat_rule)
        
        def EConsBoiler_rule(m, ts):
            return m.cons_fuel_e["boiler", ts] == (sum(m.heat["boiler", hCons, ts] for hCons in hIn) /\
                gasBoiler_eff)
        m.const_EConsBoiler = Constraint(ts_i, rule = EConsBoiler_rule)
            
        def VolConsBoiler_rule(m, ts):
            return m.cons_fuel_vol["boiler", ts] == m.cons_fuel_e["boiler", ts] / hhv["NG"]
        m.const_VolConsBoiler = Constraint(ts_i, rule = VolConsBoiler_rule)   
        
        def boilerEmiss_rule(m, ts):
            return m.emiss_CO2["boiler", ts] == m.cons_fuel_e["boiler", ts] * emiss_fac["NG"]
        m.const_boilerEmiss = Constraint(ts_i, rule = boilerEmiss_rule)

        
    '''
    Gas Boiler (Biofuel)
    '''
    
    if "boiler_b" in ders_df.index.values:
        
        def BoilerBHeat_rule(m, ts):
            return sum(m.heat["boiler_b", hCons, ts] for hCons in hIn) <= \
                m.c_c["boiler_b"]
        m.const_BoilerBCapa = Constraint(ts_i, rule = BoilerBHeat_rule)
        
        def EConsBoilerB_rule(m, ts):
            return m.cons_fuel_e["boiler_b", ts] == (sum(m.heat["boiler_b", hCons, ts] for hCons in hIn) /\
                gasBoiler_eff)
        m.const_EConsBoilerB = Constraint(ts_i, rule = EConsBoilerB_rule) 

        def VolConsBoilerB_rule(m, ts):
            return m.cons_fuel_vol["boiler_b", ts] == m.cons_fuel_e["boiler_b", ts] / hhv["BG"]
        m.const_VolConsBoilerB = Constraint(ts_i, rule = VolConsBoilerB_rule)
        
        def boilerBEmiss_rule(m, ts):
            return m.emiss_CO2["boiler_b", ts] == m.cons_fuel_e["boiler_b", ts] * emiss_fac["BG"]
        m.const_boilerBEmiss = Constraint(ts_i, rule = boilerBEmiss_rule)
    
    '''
    Solar Thermal collectors
    '''
    m.st_groundPanelArea = Var(within = NonNegativeReals, initialize = 0)
    m.st_roofPanelArea = Var(r_az, r_tilt, within = NonNegativeReals, initialize = 0)
    #Solar Thermal covered roof area at different angles

    
    # if "solar_thermal" in ders_df.index.values:
        
    #     def stUsefulHeat_rule(m,ts):          
    #         return sum(m.heat["solar_thermal", hCons, ts] for hCons in hIn) == \
    #             (st_eff_0-((st_a1*(st_Tf-temp[ts]))/(irrad_tilt_flat[ts] *1000)) - \
    #             ((st_a2*(st_Tf-temp[ts])**2)/(irrad_tilt_flat[ts] * 1000))) * \
    #             m.st_groundPanelArea * irrad_tilt_flat[ts] + \
    #             sum((st_eff_0-((st_a1*(st_Tf-temp[ts]))/(irrad_tilt_tilt_df.at[az,tilt][ts] *1000)) - \
    #             ((st_a2*(st_Tf-temp[ts])**2)/(irrad_tilt_tilt_df.at[az,tilt][ts] * 1000))) *\
    #             m.st_roofPanelArea[az, tilt] * irrad_tilt_tilt_df.at[az,tilt][ts]\
    #             for az in r_az for tilt in r_tilt)
    #     m.const_stUsefulHeat = Constraint(ts_i,rule = stUsefulHeat_rule)
            
    #     def stpanelArea_rule(m):
    #         return m.c_c["solar_thermal"] == (m.st_groundPanelArea + \
    #         sum(m.st_roofPanelArea[az, tilt] for az in r_az for tilt in r_tilt)) * 0.717 # official conversion factor 
    #     m.const_stpanelArea = Constraint(rule = stpanelArea_rule)
    
    if "solar_thermal_ground" in ders_df.index.values:
        
        def stUsefulHeatG_rule(m,ts):          
            return sum(m.heat["solar_thermal_ground", hCons, ts] for hCons in hIn) == \
                m.st_groundPanelArea*(st_eff_0*float(irrad_tilt_flat[ts]) *1000 - st_a1*(st_Tf-temp[ts]) - \
                st_a2*((st_Tf-temp[ts])**2))/1000
        m.const_stUsefulHeatG = Constraint(ts_i,rule = stUsefulHeatG_rule)
            
        def stpanelAreaG_rule(m):
            return m.c_c["solar_thermal_ground"] == m.st_groundPanelArea * 0.717 # official conversion factor 
        m.const_stpanelAreaG = Constraint(rule = stpanelAreaG_rule)
        
        if "pv_ground" in ders_df.index.values:
            def panel_groundAreaLimit(m):
                return m.pv_groundPanelArea/pv_gcr_ground + m.st_groundPanelArea*1.7 <= pv_groundareaAvail
            m.const_panelGroundAreaLimit = Constraint(rule = panel_groundAreaLimit)
        else:
            def panel_groundAreaLimit(m):
                return m.st_groundPanelArea*1.7 <= pv_groundareaAvail
            m.const_panelGroundAreaLimit = Constraint(rule = panel_groundAreaLimit)
        
    if "solar_thermal_roof" in ders_df.index.values:
        
        def stUsefulHeatR_rule(m,ts):          
            return sum(m.heat["solar_thermal_roof", hCons, ts] for hCons in hIn) == \
                sum((st_eff_0*irrad_tilt_tilt_df.at[az,tilt][ts] *1000-st_a1*(st_Tf-temp[ts]) - \
                st_a2*((st_Tf-temp[ts])**2)) * m.st_roofPanelArea[az, tilt]/1000 for az in r_az for tilt in r_tilt)
        m.const_stUsefulHeatR = Constraint(ts_i,rule = stUsefulHeatR_rule)
            
        def stpanelAreaR_rule(m):
            return m.c_c["solar_thermal_roof"] == sum(m.st_roofPanelArea[az, tilt] \
                    for az in r_az for tilt in r_tilt) * 0.717 # official conversion factor 
        m.const_stpanelAreaR = Constraint(rule = stpanelAreaR_rule)
        
        if "pv_roof" in ders_df.index.values:
            def panel_roofAreaLimit(m, az, tilt):
                return m.pv_roofPanelArea[az, tilt]/pv_gcr_roof + m.st_roofPanelArea[az, tilt] * \
                    1.7 <= pv_roofwallareaAvail.at[az, tilt]
            m.const_panelRoofAreaLimit = Constraint(r_az, r_tilt, rule = panel_roofAreaLimit)
        else:
            def panel_roofAreaLimit(m, az, tilt):
                return m.st_roofPanelArea[az, tilt] * 1.7 <= pv_roofwallareaAvail.at[az, tilt]
            m.const_panelRoofAreaLimit = Constraint(r_az, r_tilt, rule = panel_roofAreaLimit)
        
    '''
    Absorption Chiller
    '''
    
    if "abs_chiller" in ders_df.index.values:
        
        def absCold_rule(m,ts):
            return sum(m.cool["abs_chiller", cCons, ts] for cCons in cIn) == \
                sum(m.heat[elSup, "abs_chiller", ts] for hSup in hOut) *\
                    absCOP
        m.const_absCold = Constraint(ts_i,rule = absCold_rule)
         
        def absCapa_rule(m,ts): 
             return sum(m.cool["abs_chiller", cCons, ts] for cCons in cIn) <=\
                 m.c_c["abs_chiller"]
        m.const_absCapa = Constraint(ts_i,rule = absCapa_rule)
            
            
            
    
            

    return m