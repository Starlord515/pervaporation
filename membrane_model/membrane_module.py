#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 12:54:19 2026

@author: ben

Models the output mass fractions, temperature, and pressure.
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint, solve_ivp
from thermo.nrtl import *
from membrane_utility import parse_templates, wt_to_x
from scipy.optimize import minimize
import plotly.graph_objects as go

def component_membrane_flux(wfm, s, c):
    flux = c.Q * (c.psat * wfm[0] * c.gamma - s.permeate_pressure * c.permeate_y)
    return flux


def boundary_layer_flux(wfm, s, c):
    #print('a = ', (c.fw - c.pw), '\nb = ', (wfm[0] - c.pw), '\n')

    flux = s.rho_tot * c.k * np.log((wfm[0] - c.pw)/(c.fw - c.pw))
    return flux


# Numerically solve for membrane flux at a point along the length of the module.
# Need to match yl (permeate molar fraction) and weight
def solve_component_flux(s, c):

    def objective_flux_fun(wfm, s, c):
        # Flux across membrane.
        f1 = component_membrane_flux(wfm, s, c)

        # Flux across boundary layer.
        f2 = boundary_layer_flux(wfm, s, c)

        return abs(f1 / f2 - 1)
    
    bounds = [(0, 1)]
    guess =  [0.9 * c.fw]
    wfm_solution = minimize(objective_flux_fun, guess, args=(s,c), method='Nelder-Mead', bounds=bounds)
    component_flux_solution = component_membrane_flux(wfm_solution.x, s, c)
    return component_flux_solution


# Solve 2nd order ODE for mass balance, given w at z=0, and dwdz at z=0.
def component_weight_fraction_derivative(s, c, dwdz, flux):
    # Differential for mass balance, assuming steady-state. Must equal to 0.
    diff = ( s.u * dwdz + flux / (s.rho_tot * s.membrane_height) ) / c.Dax
    return diff


def permiate_mass_differential(s, flux_array):
    dmdz = s.dh * sum(flux_array) #Todo. Check what this equation means.
    return dmdz


def odes(z, x, system_data, components, dflog):
    # Assign each weight fraction ODE to an element of x.
    w_ethanol = x[0]
    dwdt_ethanol = x[1]
    w_water = x[2]
    dwdt_water = x[3]
    permeate_mass_flow = x[4]

    print('z: ', z, ' : ', x, '\n')
    
    # Update system class object.
    system_data.update(new_mass_fractions = [w_water, w_ethanol],
                       new_mass_flow = system_data.total_mass_flow - permeate_mass_flow,
                       new_T = system_data.temp,
                       new_P = system_data.pressure
                       )

    # Update component class objects.
    for comp in components: comp.update(system_data=system_data)

    # Assign dwdz for x[i] at which i is a weight fraction.
    i = 0
    for j in range(system_data.N_components):
        x[i] = x[i+1]
        i+=2

    # Calculate the component SECOND derivatives.
    flux_array = []
    i = 1
    for comp in components:
        # Solve for component membrane flux.
        flux = solve_component_flux(s=system_data, c=comp)
        flux_array.append(flux)

        # Calculate 2nd derivative, assign to dwdz.
        x[i] = component_weight_fraction_derivative(s=system_data, c=comp, dwdz=x[i], flux=flux)
        i += 2

    # Calculate permeate derivative
    x[4] = permiate_mass_differential(s=system_data, flux_array=flux_array)

    return x


    
def main():
    # Extract data from templates. Create class objects for system data and component properties.
    system_data, components = parse_templates()

    # Initial weight fractions, and second derivatives are 0.
    x0 = []
    for wt in system_data.mass_fractions:
        x0.append(wt)
        x0.append(wt)

    # Append initial permiate mass flow (0)
    x0.append(0)

    # Space of z values, along length of reactor.
    z = np.linspace(0, system_data.l, 10)

    # Create dataframe to store values
    dflog = pd.DataFrame()

    # Solve ODEs.
    #sol = odeint(odes, y0=x0, t=z, args=(system_data, components, dflog))
    sol = solve_ivp(fun = odes, t_span=(0, system_data.l), y0=x0, method='RK45', args=(system_data, components, dflog))

    # Extract each component solution.
    df = pd.DataFrame()
    df['z'] = sol.t
    i = 0
    fig = go.Figure()
    for comp in components:
        df[comp.label] = sol.y[i]

        # Show as figure.
        fig.add_trace(go.Scatter(x=sol.t, y=sol.y[i], name=comp.label))

        i += 2

    fig.add_trace(go.Scatter(x=sol.t, y=sol.y[-1], name='Permeate Flowrate', yaxis='y2', mode='lines'))
    fig.show()
    print(df)
    
    return 0

main()