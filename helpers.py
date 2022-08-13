# -*- coding: utf-8 -*-
"""
Created on Wed May 20 12:03:42 2020

@author: Axel
"""

import pandas as pd
import numpy as np
from pvlib import solarposition
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import seaborn as sns

def ders_consider(ders):
    '''
    Deletes the rows with DERs not considered from ders panda

    Parameters
    ----------
    ders : panda Dataframe
        Dataframe with distributed generators considerations as discribed in 
        data_prep.

    Returns
    -------
    dersActive : panda Dataframe
        adapted ders only containing those DERs to be considered.

    '''
    
    dersActive = ders[ders["Consider"] == 1]
    
    return dersActive


def irrad_tilt(time, irrad, loc, pan_az, pan_tilt):
    '''
    Calculates the tilted irradiance at current timestep

    Parameters
    ----------
    time : Panda DatetimeIndex
        DatetimeIndex of current timestep.
    irrad : Panda Dataframe
        Dataframe with global horizontal, direct normal and diffuse horrizontal
        irradiance in W/m2
    loc : Dict
        Python dictionary with values for the follwing keys:
            latitude: "lat" --> float
            longitude: "long" --> float
            altitude: "alt" --> float
            av. temperature: "temp" --> float
            albedo: "alb". --> float
    pan_az : Int
        Panel/roof azimuth angle in °.
    pan_tilt : Int
        Panel/roof tilt angle inn °.

    Returns
    -------
    irrad_tilt : Float
        Tilted irradiance in kW/m2 for given surface and timestep.

    '''
    
    position_sun = solarposition.get_solarposition(time, loc['lat'], loc['long'], \
                                                   altitude = loc['alt'], pressure = None, \
                                                       temperature = loc['temp'])
    sun_zenith = position_sun.at[time, 'zenith']
    sun_azimuth = position_sun.at[time, 'azimuth']
    
    cos_tilt = np.cos(np.radians(sun_zenith)) * np.cos(np.radians(pan_tilt)) + \
        np.sin(np.radians(sun_zenith)) * np.sin(np.radians(pan_tilt)) * \
            np.cos(np.radians(sun_azimuth) - np.radians(pan_az))
    
    tilt_dir = irrad.at["direct_n"] * np.maximum(0 , cos_tilt / \
                        np.sin(np.radians(90-sun_zenith)))
    tilt_diff = irrad.at["diffuse_h"] * (1 + np.cos(np.radians(pan_tilt))) * 0.5
    tilt_refl = irrad.at["global_h"] * loc['alb'] * (1 - np.cos(np.radians(pan_tilt))) * 0.5
    
    irrad_tilt = (tilt_diff + tilt_dir + tilt_refl) / 1000
    
    return irrad_tilt

def irrad_tiltMW(time, irrad, loc, pan_az, pan_tilt):
    '''
    Calculates the tilted irradiance at current timestep

    Parameters
    ----------
    time : Panda DatetimeIndex
        DatetimeIndex of current timestep.
    irrad : Panda Dataframe
        Dataframe with global horizontal, direct normal and diffuse horrizontal
        irradiance in W/m2
    loc : Dict
        Python dictionary with values for the follwing keys:
            latitude: "lat" --> float
            longitude: "long" --> float
            altitude: "alt" --> float
            av. temperature: "temp" --> float
            albedo: "alb". --> float
    pan_az : Int
        Panel/roof azimuth angle in °.
    pan_tilt : Int
        Panel/roof tilt angle inn °.

    Returns
    -------
    irrad_tilt : Float
        Tilted irradiance in MW/m2 for given surface and timestep.

    '''
    
    position_sun = solarposition.get_solarposition(time, loc['lat'], loc['long'], \
                                                   altitude = loc['alt'], pressure = None, \
                                                       temperature = loc['temp'])
    sun_zenith = position_sun.at[time, 'zenith']
    sun_azimuth = position_sun.at[time, 'azimuth']
    
    cos_tilt = np.cos(np.radians(sun_zenith)) * np.cos(np.radians(pan_tilt)) + \
        np.sin(np.radians(sun_zenith)) * np.sin(np.radians(pan_tilt)) * \
            np.cos(np.radians(sun_azimuth) - np.radians(pan_az))
    
    tilt_dir = irrad.at["direct_n"] * np.maximum(0 , cos_tilt / \
                        np.sin(np.radians(90-sun_zenith)))
    tilt_diff = irrad.at["diffuse_h"] * (1 + np.cos(np.radians(pan_tilt))) * 0.5
    tilt_refl = irrad.at["global_h"] * loc['alb'] * (1 - np.cos(np.radians(pan_tilt))) * 0.5
    
    irrad_tilt = (tilt_diff + tilt_dir + tilt_refl) / 1000000
    
    return irrad_tilt

def i_PyomoValues(var):
    array = np.fromiter(var.get_values().values(), dtype = float)
    return array

def lineequ(aX, aY, bX, bY):
    m = (aY-bY)/(aX-bX)
    t = aY-(m*aX)
    return m, t

def tariffsizer(tariff_df, timesteps, dt):
    """
    

    Parameters
    ----------
    tariff_df : Pandas
        Dataframe with tariff (and feedin tariff).
    timesteps : int
        Nr of timesteps in model
    dt : int
        Resolution of model
    
    Returns
    -------
    Df with tariffs to fit to timesteps.

    """
    
    if tariff_df.index.size == 24:
        tariff_dt = tariff_df.loc[tariff_df.index.repeat(60/dt)].reset_index(drop = True)
        tariff = pd.concat([tariff_dt]*int(timesteps/tariff_dt.index.size),ignore_index = True)
    elif tariff_df.index.size == 48:
        tariff_dt = tariff_df.loc[tariff_df.index.repeat(30/dt)].reset_index(drop = True)
        tariff = pd.concat([tariff_dt]*int(timesteps/tariff_dt.index.size),ignore_index = True)
    elif tariff_df.index.size == 96:
        tariff = pd.concat([tariff_df]*int(timesteps/tariff_df.index.size),ignore_index = True)
    elif tariff_df.index.size == 8760 or 17520:
        tariff = pd.concat([tariff_df]*int(timesteps/tariff_df.index.size),ignore_index = True)
    
    return tariff


def irrad_tilt_v2(time, irrad, loc, pan_az, pan_tilt, resolution):
    '''
    Calculates the tilted irradiance at current timestep

    Parameters
    ----------
    time : Panda DatetimeIndex range
        Complete date range
    irrad : Panda Dataframe
        Dataframe with global horizontal, direct normal and diffuse horrizontal
        irradiance in W/m2
    loc : Dict
        Python dictionary with values for the follwing keys:
            latitude: "lat" --> float
            longitude: "long" --> float
            altitude: "alt" --> float
            av. temperature: "temp" --> float
            albedo: "alb". --> float
    pan_az : Int
        Panel/roof azimuth angle in °.
    pan_tilt : Int
        Panel/roof tilt angle inn °.

    Returns
    -------
    irrad_tilt : np.array
        Tilted irradiance in kW/m2 for given surface and timestep.

    '''
    irrad_tilt = np.empty(len(irrad.index))
    for i in range(len(irrad.index)):        
    
        position_sun = solarposition.get_solarposition(time[i], loc['lat'], loc['long'], \
                                                       altitude = loc['alt'], pressure = None, \
                                                           temperature = loc['temp'])
        sun_zenith = position_sun.at[time[i], 'zenith']
        sun_azimuth = position_sun.at[time[i], 'azimuth']
    
        cos_tilt = np.cos(np.radians(sun_zenith)) * np.cos(np.radians(pan_tilt)) + \
            np.sin(np.radians(sun_zenith)) * np.sin(np.radians(pan_tilt)) * \
                np.cos(np.radians(sun_azimuth) - np.radians(pan_az))
        iDir = irrad.at[i,"direct_n"]
        if iDir <= 1:
            tilt_dir = 0
        elif sun_zenith > 85 or sun_zenith < 0:
            tilt_dir = 0
        else:
            tilt_dir = iDir * np.maximum(0 , cos_tilt / \
                            np.sin(np.radians(90-sun_zenith)))
        iDiff = irrad.at[i,"diffuse_h"]
        if iDiff <= 1:
            tilt_diff = 0
        else:
            tilt_diff = iDiff * (1 + np.cos(np.radians(pan_tilt))) * 0.5
        iGlob = irrad.at[i,"global_h"]
        if iGlob <= 1:
            tilt_refl = 0
        else:
            tilt_refl = iGlob * loc['alb'] * (1 - np.cos(np.radians(pan_tilt))) * 0.5
    
        irrad_tilt[i] = (tilt_diff + tilt_dir + tilt_refl) / 1000
    
    
    return np.repeat(irrad_tilt,60/resolution)

def irrad_tilt_v3(time, irrad, loc, pan_az, pan_tilt, resolution):
    '''
    Calculates the tilted irradiance at current timestep

    Parameters
    ----------
    time : Panda DatetimeIndex range
        Complete date range
    irrad : Panda Dataframe
        Dataframe with global horizontal, direct normal and diffuse horrizontal
        irradiance in W/m2
    loc : Dict
        Python dictionary with values for the follwing keys:
            latitude: "lat" --> float
            longitude: "long" --> float
            altitude: "alt" --> float
            av. temperature: "temp" --> float
            albedo: "alb". --> float
    pan_az : Int
        Panel/roof azimuth angle in °.
    pan_tilt : Int
        Panel/roof tilt angle inn °.

    Returns
    -------
    irrad_tilt : np.array
        Tilted irradiance in kW/m2 for given surface and timestep.

    '''       
    
    position_sun = solarposition.get_solarposition(time, loc['lat'], loc['long'], \
                                                   altitude = loc['alt'], pressure = None, \
                                                       temperature = loc['temp'])
    sun_zenith = position_sun.zenith.values
    sun_azimuth = position_sun.azimuth.values

    cos_tilt = np.cos(np.radians(sun_zenith)) * np.cos(np.radians(pan_tilt)) + \
        np.sin(np.radians(sun_zenith)) * np.sin(np.radians(pan_tilt)) * \
            np.cos(np.radians(sun_azimuth) - np.radians(pan_az))

    
    tilt_dir = irrad.direct_n.values * (np.maximum(0 , cos_tilt / \
                        np.sin(np.radians(90-sun_zenith))))

    high_values_flags = sun_zenith > 85   
    tilt_dir[high_values_flags] = 0  # All low values set to 0
    low_values_flags = sun_zenith < 0   
    tilt_dir[low_values_flags] = 0  # All low values set to 0    
    
    tilt_diff = irrad.diffuse_h.values * (1 + np.cos(np.radians(pan_tilt))) * 0.5

    tilt_refl = irrad.global_h.values * loc['alb'] * (1 - np.cos(np.radians(pan_tilt))) * 0.5

    irrad_t = (tilt_diff + tilt_dir + tilt_refl) / 1000
    
    irrad_tilt = np.asarray(irrad_t)
    low_values_flags = irrad_tilt < 0.01  # Where values are low
    irrad_tilt[low_values_flags] = 0  # All low values set to 0
    
    return np.repeat(irrad_tilt,60/resolution)

def mVal(m):
    return m.obj.expr()

def varVal(Var):
    return Var.get_values()

def statusQ(m):
    d = {}
    d["obj"] = mVal(m)
    d["cost_I"] = m.cost_I.value
    x = varVal(m.cost_OM_fix) 
    d["cost_OMfix"] = np.array(list(x.values()))
    x = varVal(m.cost_OM_var)
    d["cost_OMvar"] = np.array(list(x.values()))
    x = varVal(m.revenue)
    d["rev"] = np.array(list(x.values()))
    x = varVal(m.gridImport)
    d["gridIn"] = np.array(list(x.values()))
    
    return d

def PED(m):
    d = {}
    d["obj"] = mVal(m)
    d["cost_I"] = m.cost_I.value
    x = varVal(m.cost_OM_fix) 
    d["cost_OMfix"] = np.array(list(x.values()))
    x = varVal(m.cost_OM_var)
    d["cost_OMvar"] = np.array(list(x.values()))
    x = varVal(m.cost_fuel) 
    d["cost_fuel"] = np.array(list(x.values()))
    x = varVal(m.cost_CO2)
    d["cost_CO2"] = np.array(list(x.values()))
    x = varVal(m.cost_Ext)
    d["cost_Ext"] = np.array(list(x.values()))
    x = varVal(m.revenue)
    d["rev"] = np.array(list(x.values()))
    d["capa_c"] = varVal(m.c_c)
    d["power"] = varVal(m.power)
    d["heat"] = varVal(m.heat)
    d["cool"] = varVal(m.cool)
    d["elEner"] = varVal(m.elEner)
    d["hEner"] = varVal(m.hEner)
    x = varVal(m.gridImport)
    d["gridIn"] = np.array(list(x.values()))
    x = varVal(m.gridImportPE)
    d["gridInPE"] = np.array(list(x.values()))
    x = varVal(m.gridExport)
    d["gridEx"] = np.array(list(x.values()))
    x = varVal(m.gridExportPE)
    d["gridExPE"] = np.array(list(x.values()))
    x = varVal(m.selfCons)
    d["selfCons"] = np.array(list(x.values()))
    x = varVal(m.selfCons_annual)
    d["selfCons_annual"] = np.array(list(x.values()))
    x = varVal(m.maxgridImport)
    d["gridInMax"] = np.array(list(x.values()))
    x = varVal(m.emiss_CO2)
    d["emiss_CO2"] = np.array(list(x.values())) 
    x = varVal(m.emiss_CO2_annual)
    d["emiss_CO2_annual"] = np.array(list(x.values()))
    x = varVal(m.cons_fuel_vol)
    d["cons_fuel_vol"] = np.array(list(x.values()))
    x = varVal(m.cons_fuel_e)
    d["cons_fuel_e"] = np.array(list(x.values()))
    x = varVal(m.emiss_PM)
    d["emiss_PM"] = np.array(list(x.values()))
    # x = varVal(m.pv_Gen)
    # d["pv_Gen"] = np.array(list(x.values()))
    # x = varVal(m.pv_groundGen)
    # d["pv_groundGen"] = np.array(list(x.values()))
    # x = varVal(m.pv_roofGen)
    # d["pv_roofGen"] = np.array(list(x.values()))
    # x = varVal(m.pv_wallGen)
    # d["pv_wallGen"] = np.array(list(x.values()))
    d["pv_groundPanelArea"] = m.pv_groundPanelArea.value
    # d["pv_roofGen"] = varVal(m.pv_roofGen)
    # d["pv_wallGen"] = varVal(m.pv_wallGen)
    x = varVal(m.batt_soc)
    d["batt_soc"] = np.array(list(x.values()))
    x = varVal(m.tss_soc)
    d["tss_soc"] = np.array(list(x.values()))
    return d

def PED_ev(m):
    d = {}
    d["obj"] = mVal(m)
    d["cost_I"] = m.cost_I.value
    x = varVal(m.cost_OM_fix) 
    d["cost_OMfix"] = np.array(list(x.values()))
    x = varVal(m.cost_OM_var)
    d["cost_OMvar"] = np.array(list(x.values()))
    x = varVal(m.cost_fuel) 
    d["cost_fuel"] = np.array(list(x.values()))
    x = varVal(m.cost_CO2)
    d["cost_CO2"] = np.array(list(x.values()))
    x = varVal(m.cost_Ext)
    d["cost_Ext"] = np.array(list(x.values()))
    x = varVal(m.revenue)
    d["rev"] = np.array(list(x.values()))
    d["capa_c"] = varVal(m.c_c)
    d["capa_d"] = varVal(m.c_d)
    d["power"] = varVal(m.power)
    d["heat"] = varVal(m.heat)
    d["cool"] = varVal(m.cool)
    d["elEner"] = varVal(m.elEner)
    d["hEner"] = varVal(m.hEner)
    x = varVal(m.gridImport)
    d["gridIn"] = np.array(list(x.values()))
    x = varVal(m.gridImportPE)
    d["gridInPE"] = np.array(list(x.values()))
    x = varVal(m.gridExport)
    d["gridEx"] = np.array(list(x.values()))
    x = varVal(m.gridExportPE)
    d["gridExPE"] = np.array(list(x.values()))
    x = varVal(m.selfCons)
    d["selfCons"] = np.array(list(x.values()))
    x = varVal(m.selfCons_annual)
    d["selfCons_annual"] = np.array(list(x.values()))
    x = varVal(m.maxgridImport)
    d["gridInMax"] = np.array(list(x.values()))
    x = varVal(m.emiss_CO2)
    d["emiss_CO2"] = np.array(list(x.values())) 
    x = varVal(m.emiss_CO2_annual)
    d["emiss_CO2_annual"] = np.array(list(x.values()))
    x = varVal(m.cons_fuel_vol)
    d["cons_fuel_vol"] = np.array(list(x.values()))
    x = varVal(m.cons_fuel_e)
    d["cons_fuel_e"] = np.array(list(x.values()))
    x = varVal(m.emiss_PM)
    d["emiss_PM"] = np.array(list(x.values()))
    x = varVal(m.soc_ev)
    d["soc_ev"] = np.array(list(x.values()))
    d["pv_groundPanelArea"] = m.pv_groundPanelArea.value
    x = varVal(m.batt_soc)
    d["batt_soc"] = np.array(list(x.values()))
    d["cs1"] = varVal(m.b_cs1)
    d["cs2"] = varVal(m.b_cs2)
    d["cs3"] = varVal(m.b_cs3)
    return d



def capaInstalled(d):
    for key in d["capa_c"]:
        print(key + ": " + str(d["capa_c"][key]))

def exImPower(d,dt, filename):
    dic = {}
    exp = d["gridEx"][:int((8760*60/dt))]
    im = d["gridIn"][:int((8760*60/dt))]
    dic["minEx"] = np.min(exp)
    dic["maxEx"] = np.max(exp)
    dic["share0Ex"] = (len(exp)-np.count_nonzero(exp))/len(exp)
    dic["minIm"] = np.min(im)
    dic["maxIm"] = np.max(im)
    dic["share0Im"] = (len(im)-np.count_nonzero(im))/len(im)
    # plt.hist(exp[exp!=0],20,label = "Export Power", color = "#20a9a6")
    # plt.ylabel("Appearance Annual")
    # plt.xlabel("Power [kW]")
    # plt.legend()
    # plt.savefig("exporthist.png",dpi = 300)
    # plt.hist(im[im!=0],20,label = "Import Power", color = "#ff6361")
    # plt.legend()
    # plt.savefig("hist_Comp-png",dpi=300)
    # plt.clf()
    # plt.hist(im[im!=0],30, label = "Import Power")
    # plt.ylabel("Appearance Annual")
    # plt.xlabel("Power [kW]")
    # plt.legend()
    # plt.savefig("importhist.png",dpi = 300)
    # plt.clf()
    # dens = sns.kdeplot(exp[exp!=0], label = "Export Power").get_figure()
    # plt.ylabel("Appearance Density")
    # plt.xlabel("Power [kW]")
    # plt.legend()
    # dens.savefig("exportdens.png",dpi = 300)
    # dens = sns.kdeplot(im[im!=0], label = "Import Power").get_figure()
    # dens.savefig("imExComparedens.png",dpi = 300)
    # plt.clf()
    # dens = sns.kdeplot(im[im!=0], label = "Import Power").get_figure()
    # plt.ylabel("Appearance Density")
    # plt.xlabel("Power [kW]")
    # plt.legend()
    # dens.savefig("importdens.png",dpi = 300)
    # plt.clf()
    
    fig, ax = plt.subplots()
    ax.grid()
    sns.distplot(im[im!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
             color="#ff6361", label = "Import Power")
    for p in ax.patches:  # turn the histogram upside down
        p.set_height(-p.get_height())
    for l in ax.lines:  # turn the kde curve upside down
        l.set_ydata(-l.get_ydata())
    sns.distplot(exp[exp!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
              color="#20a9a6", label = "Export Power")
    ax.set_xticks(np.arange(0, round((dic["maxEx"] + 10)/10)*10, 40))
    ax.set_yticks(np.arange(0.0, 0.08, 0.01))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # pos_ticks = np.array([t for t in ax.get_yticks() if t > 0])
    # ticks = np.concatenate([-pos_ticks[::-1], [0], pos_ticks])
    ticks = np.arange(-0.07,0.03,0.01)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{abs(t):.2f}' for t in ticks])
    ax.spines['bottom'].set_position('zero')
    ax.set_ylabel("Density",loc = "center")
    ax.set_xlabel("Grid Exchange Power [kW]",loc = "right")
    ax.set_title("Distribution of grid exchange power observations")
    ax.legend(loc = 'lower right')
    fig.savefig(filename + ".png",dpi = 300)
    return dic

def exImPowerComp(d,dd,dt, filename):
    dic = {}
    exp1 = d["gridEx"][:int((8760*60/dt))]
    im1 = d["gridIn"][:int((8760*60/dt))]
    maxEx1 = np.max(exp1)
    exp2 = dd["gridEx"][:int((8760*60/dt))]
    im2 = dd["gridIn"][:int((8760*60/dt))]
    maxEx2 = np.max(exp2)
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    sns.distplot(im1[im1!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
             color="#ff6361", label = "Import Power", ax = ax1)
    for p in ax1.patches:  # turn the histogram upside down
        p.set_height(-p.get_height())
    for l in ax1.lines:  # turn the kde curve upside down
        l.set_ydata(-l.get_ydata())
    sns.distplot(exp1[exp1!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
              color="#20a9a6", label = "Export Power", ax = ax1)
    ax1.set_xticks(np.arange(0, round((maxEx1 + 10)/10)*10, 40))
    ax1.set_yticks(np.arange(0.0, 0.08, 0.01))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # pos_ticks = np.array([t for t in ax.get_yticks() if t > 0])
    # ticks = np.concatenate([-pos_ticks[::-1], [0], pos_ticks])
    ticks = np.arange(-0.07,0.03,0.01)
    ax1.set_yticks(ticks)
    ax1.set_yticklabels([f'{abs(t):.2f}' for t in ticks])
    ax1.spines['bottom'].set_position('zero')
    ax1.set_ylabel("Density",loc = "center")
    ax1.set_xlabel("Grid Exchange Power [kW]",loc = "right")
    ax1.set_title("Distribution of grid exchange power observations")
    ax1.legend(loc = 'lower right')
    
    sns.distplot(im2[im2!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
             color="#ff6361", label = "Import Power", ax = ax2)
    for p in ax2.patches:  # turn the histogram upside down
        p.set_height(-p.get_height())
    for l in ax2.lines:  # turn the kde curve upside down
        l.set_ydata(-l.get_ydata())
    sns.distplot(exp2[exp2!=0], hist=True, kde=True, hist_kws={'edgecolor': 'black'}, kde_kws={'linewidth': 2}, bins=10,
              color="#20a9a6", label = "Export Power", ax = ax2)
    ax2.set_xticks(np.arange(0, round((maxEx2 + 10)/10)*10, 100))
    ax2.set_yticks(np.arange(0.0, 0.08, 0.01))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # pos_ticks = np.array([t for t in ax.get_yticks() if t > 0])
    # ticks = np.concatenate([-pos_ticks[::-1], [0], pos_ticks])
    ticks = np.arange(-0.07,0.03,0.01)
    ax2.set_yticks(ticks)
    ax2.set_yticklabels([f'{abs(t):.2f}' for t in ticks])
    ax2.spines['bottom'].set_position('zero')
    ax2.set_ylabel("Density",loc = "center")
    ax2.set_xlabel("Grid Exchange Power [kW]",loc = "right")
    ax2.set_title("Distribution of grid exchange power observations")
    ax2.legend(loc = 'lower right')
    fig.legend()
    fig.savefig(filename + ".png",dpi = 300)
    
    
    
    
def NPVCost(d):
    """
    

    Parameters
    ----------
    d : dict
        Output of PED/statusQ functions with all saved variables and results
        of optimization.

    Returns
    -------
    dic : dict
        dictonary with all relevant cost component of optimized solution.

    """
    dic = {}
    dic["fix"] = d["cost_OMfix"]
    dic["var"] = d["cost_OMvar"]
    #dic["fuel"] = d["cost_fuel"]
    #dic["CO2"] = d["cost_CO2"]    
    #dic["Ext"] = d["cost_Ext"]
    dic["rev"] = d["rev"]
    dicc = {}
    for key in dic:
        dicc[key + "_share"] = sum((dic[key][y])*(1/(1.05)**(y+1)) for y in range(len(dic[key])))
    dic = {**dic, **dicc}
    dic["NPV"] = d["obj"]
    dic["I"] = d["cost_I"]
    return dic
    
def plotNPV(d, names, filename):
    labels = names
    x = np.arange(len(labels))
    width = 0.35
    NPV = np.empty(0)
    inv = np.empty(0)
    fix = np.empty(0)
    var = np.empty(0)
    rev = np.empty(0)
    for i in d:
        NPV = np.append(NPV,i["NPV"]/1000000)
        inv = np.append(inv,-i["I"]/1000000)
        fix = np.append(fix,-i["fix_share"]/1000000)
        var = np.append(var,-i["var_share"]/1000000)
        rev = np.append(rev,i["rev_share"]/1000000)
    
    sns.set_theme(style='darkgrid')
    fig, ax = plt.subplots()
    ax.grid() 
    ax.bar(x - width/2,inv,width,label = "CAPEX", color="#bc5090")
    ax.bar(x - width/2,fix,width, bottom=inv, label = "Fix Costs", color="#58508d")
    ax.bar(x - width/2,var,width, bottom=(fix+inv), label = "Variable Costs",color="#ff6361")
    ax.bar(x - width/2,rev,width, bottom=0, label = "Revenue", color="#20a9a6")
    ax.bar(x + width/2,NPV,width,label = "NPV", color="#ffa600")
    
    ax.set_ylabel("NPV Composition [Mio €]")
    ax.set_title("NPV Composition by Scenario")
    for i in x:
        ax.text(i+width/2, NPV[i]-0.1, str(np.around(NPV[i],2)), color='black', ha = "center", size = "x-small")
    ax.set_xticks(x)
    ax.tick_params(axis='x', rotation=45)
    ax.set_xticklabels(labels)
    ax.grid(linestyle='-', linewidth='0.5', color='white')

    fontP = FontProperties()
    fontP.set_size('xx-small')
    ax.legend(prop=fontP)
    fig.subplots_adjust(bottom=0.30)
    fig.savefig(filename + ".png", dpi = 300)
    

def ex2imRatio(d):
    return sum(d["li"])/sum(d["gridIn"])

def ex2imRatioPE(d):
    return sum(d["gridExPE"])/sum(d["gridInPE"])

def CO2em(d, in_timeseries, CO2_factors):
    """

    Parameters
    ----------
    d : dict
        Result of model.
    in_timeseries : string
        link to timeseries file that includes generation mix.
    CO2_factors : dict
        dict with CO2 factors that have the same keys as gen mix file.

    Returns
    -------
    CO2 emissions over whole time

    """
    CO2_fac = []
    im = d["gridIn"]
    gen = pd.read_excel(in_timeseries, sheet_name = "GridMixData")
    for row in gen.index:
        dd = CO2_factors.copy()
        for key in CO2_factors:
            dd[key] = gen.at[row,key]/gen.at[row,"total"]*CO2_factors[key]
        CO2_fac.append(sum(dd.values()))
    years = int(len(im)/len(gen))
    CO2_fac = CO2_fac * years
    CO2_fac = np.array(CO2_fac)
    return sum(CO2_fac*im)

def plotLoadCov(d, in_timeseries, in_tech,filename1,filename2,filename3):
    plt.style.use('seaborn')
    ders_total = pd.read_excel(in_tech, sheet_name = "Selection&Cost")
    ders_df = ders_total[ders_total["Consider"] == 1]
    ders_df = ders_df.set_index(['DER'], drop = False)
    
    ders_c = ders_df[ders_df['continously'] == 1].index.values
    ders_d = ders_df[ders_df['continously'] == 0].index.values
    #All distributed energy sources that can output electricity
    d_elOut = ders_df[ders_df['el_out'] == 1].index.values
    
    #All considered tech (including grid) that can output and take in electrictiy
    elOut = np.append(d_elOut, 'grid')
    elIn = ders_df[ders_df['el_in'] == 1].index.values
    elIn = np.append(elIn, ['grid', 'load_e'])
    d_elGen = d_elOut[d_elOut!='battery']
    
    loadData = pd.read_excel(in_timeseries, sheet_name = "LoadData")
    load_e_in = loadData.loc[:,"Load_e [kW]"].values
    
    pv = d["pv_Gen"]
    soc = d["batt_soc"]
    power = d["power"]



    grid2bat = []
    bat2grid = []
    bat2X = []
    X2load = {}
    X2X = {}
    batP = {
      "pv2bat": [],
      "grid2bat": [],
      "bat2load": [],
      "bat2grid": []
    }
    for out in elOut:
        X2load[out] = []
        X2X[out] = []
    X2X["pv2grid"] = []
    X2X["pv2bat"] = []
    pv2bat = []
    pv2X = {}
    for elin in elIn:
        pv2X[elin] = []
    for i in range(8760):
        grid2bat.append(power["grid","battery",1,i])
        bat2grid.append(power["battery","grid",1,i])
        bat2X.append(sum(power["battery",ein,1,i] for ein in elIn))
        batP["pv2bat"].append(sum(power[gen,"battery",1,i] for gen in d_elGen))
        batP["grid2bat"].append(power["grid","battery",1,i])
        batP["bat2load"].append((power["battery","load_e",1,i])*(-1))
        batP["bat2grid"].append((power["battery","grid",1,i])*(-1))
        for out in elOut:
            X2load[out].append(power[out,"load_e",1,i])
            X2X[out].append(power[out,"load_e",1,i])
        X2X["pv2grid"].append(sum(power[gen,"grid",1,i] for gen in d_elGen))
        X2X["pv2bat"].append(sum(power[gen,"battery",1,i] for gen in d_elGen))
        pv2bat.append(sum(power[gen,"battery",1,i] for gen in d_elGen))
        for elin in elIn:
            pv2X[elin].append(sum(power[gen,elin,1,i] for gen in d_elGen))
    grid2bat = np.array(grid2bat)
    bat2grid = np.array(bat2grid)*(-1)
    bat2X = np.array(bat2X)*(-1)
    pv2bat = np.array(pv2bat)

    elInn = elIn[elIn != "grid"]
    
    x2load = pd.DataFrame(X2load)
    x2load["pv2load"] = x2load["pv_roof"] + x2load["pv_ground"]
    x2load = x2load.drop(['pv_roof', 'pv_ground'], axis=1)
    x2load = x2load[["pv2load", "battery", "grid"]]
    X2X = pd.DataFrame(X2X)
    X2X["pv2load"] = X2X["pv_roof"] + X2X["pv_ground"]
    x2x = X2X.drop(['pv_roof', 'pv_ground'], axis=1)
    cols = x2x.columns.tolist()
    colsnew = ['grid', 'battery', 'pv2load', 'pv2grid', 'pv2bat']
    x2x = x2x[colsnew]
    
    
    x2loadPlot = x2load[5640:5688]
    x2loadPlot = x2loadPlot.reset_index(drop = True)
    x2xPlot = x2x[5640:5688]
    x2xPlot = x2xPlot.reset_index(drop = True)
    
    PV2X = pd.DataFrame(pv2X)
    colsnew = ["load_e", "battery", "grid"]
    PV2X = PV2X[colsnew]
    PV2XPlot = PV2X[5640:5688]
    PV2XPlot = PV2XPlot.reset_index(drop = True)
    pv = pv[5640:5688]
    
    batP.keys()
    bat = pd.DataFrame(batP)
    batPlot = bat[5640:5688]
    batPlot = batPlot.reset_index(drop = True)


    fig, ax = plt.subplots()
    sns.set()
    ax.grid()
    line0 = sns.lineplot(ax = ax, data = pv, linewidth=0, color="#ffa600")
    areas = x2loadPlot.plot.area(ax = ax, linewidth=0, color={'grid':'#bc5090', 'battery':'#20a9a6', 'pv2load':'#ffa600'})
    line = sns.lineplot(ax = ax, data = load_e_in[5640:5688], linewidth=3, label = "load")
    line2 = sns.lineplot(ax = ax, data = pv, linewidth=1, color="#ffa600", label = "PV power")
    ax.set_xlabel('Hour of day', fontsize=15)
    ax.set_ylabel('Electricity Generation [kW]', fontsize=15)
    ax.set_title('Load coverage',fontsize=17)
    ax.set_xticks([0,4,8,12,16,20,24,28,32,36,40,44,48])
    
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [0,4,8,12,16,20,0,4,8,12,16,20,0]
    
    ax.set_xticklabels(labels)
    fig.savefig(filename1 + ".png", dpi = 300)
    
    
    fig, ax = plt.subplots()
    sns.set()
    areas = PV2XPlot.plot.area(ax = ax, linewidth=0 , color={'load_e':'#ffa600', 'battery':'#c3800e', 'grid':'#996512'})
    line = sns.lineplot(ax = ax, data = load_e_in[5640:5688], linewidth=3, label = "load")
    ax.set_xlabel('Hour of day', fontsize=15)
    ax.set_ylabel('PV Power Generation [kW]', fontsize=15)
    ax.set_title('PV Power Generation Distribution',fontsize=17)
    
    ax.set_xticks([0,4,8,12,16,20,24,28,32,36,40,44,48])
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [0,4,8,12,16,20,0,4,8,12,16,20,0]
    ax.set_xticklabels(labels)
    fig.savefig(filename2 + ".png", dpi = 300)
    
    fig, ax = plt.subplots() 
    sns.set()
    df_neg, df_pos = batPlot.clip(upper=0), batPlot.clip(lower=0)
    df_pos.plot.area(ax=ax, stacked=True, linewidth=0, color={'pv2bat':'#ffa600', 'grid2bat':'#20a9a6', 'bat2load':'#bc5090', 'bat2grid':'#dcdcdc'})
    ax.set_prop_cycle(None)
    df_neg.rename(columns=lambda x: '_' + x).plot.area(ax=ax, stacked=True, linewidth=0, color = '#bc5090')
    ax.set_ylim([df_neg.sum(axis=1).min(), df_pos.sum(axis=1).max()])
    ax2 = ax.twinx()
    sns.lineplot(data = soc[5640:5688], palette = ['#bc5090'], ax = ax2, linewidth=3)
    ax.set_title("Battery Charge and Discharge Power", fontsize=17)
    ax.set_xlabel ("Hour of day")
    ax.set_ylabel ("Battery Charge/Discharge [kW]")
    ax2.set_ylabel ("State of Charge [kWh]")
    ax.legend(loc='upper right')
    ax2.grid(False)
    ax2.legend(loc='upper left',labels = ["SOC"])
    
    ax.set_xticks([0,4,8,12,16,20,24,28,32,36,40,44,48])
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [0,4,8,12,16,20,0,4,8,12,16,20,0]
    ax.set_xticklabels(labels)
    fig.savefig(filename3 + ".png", dpi = 300)
    