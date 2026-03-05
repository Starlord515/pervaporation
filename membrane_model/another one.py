#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 12:54:19 2026

@author: Ben Alexander

Models the output mass fractions.
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from membrane_utility import parse_templates, wt_to_x
from scipy.optimize import minimize
import plotly.graph_objects as go


def component_membrane_flux(wfm_x, s, c):
    # Calculate gammas for membrane surface.
    ge = s.ge
    ge.to_T_xs(s.temp, [wfm_x, 1 - wfm_x])
    gamma = ge.gammas()[c.n]

    flux = c.Q * (c.psat * wfm_x * gamma - s.permeate_pressure * c.permeate_y)
    return flux


def boundary_layer_flux(wfm, s, c):
    flux = s.rho_tot * c.k * np.log((wfm - c.pw) / (c.fw - c.pw))
    return flux


# Numerically solve for membrane flux at a point along the length of the module.
# Need to match yl (permeate molar fraction) and weight
def solve_component_flux(s, components):
    def obj_fun_wfm(wfm, s, components):
        boundary_fluxes = []
        for comp in components:
            boundary_fluxes.append(boundary_layer_flux(wfm[comp.n], s, comp))
        # print(wfm,' : ',sum(wfm), '=====',boundary_fluxes)
        flux_range = (np.max(boundary_fluxes) - np.min(boundary_fluxes)) ** 2 * 100000
        return flux_range

    # objective_flux_fun = lambda wfm:  - component_membrane_flux(wfm, s, c)
    bounds = [(0, 1) for i in range(len(components))]
    guess = [0.9 * c.fw for c in components]

    # Find solution to mass fraction at membrane surface.
    constraint = {'type': 'eq', 'fun': lambda wfm: sum(wfm) - 1}
    wfm_solution = minimize(obj_fun_wfm,
                            guess,
                            method='Nelder-Mead',
                            constraints=[constraint],
                            bounds=bounds,
                            args=(s, components))

    # Calculate flux across membrane for each component.
    flux_array = []
    wfm_x = wt_to_x(wfm_solution.x, s.rmm_i)
    for comp in components:
        flux_array.append(component_membrane_flux(wfm_x[comp.n], s, comp))
    print('wfm = ', wfm_solution.x, '        flux = ', flux_array)

    # Convert membrane surface wt fractions to mole fractions. Calculate J_tot.

    J_tot = boundary_layer_flux(wfm_solution.x[0], s, components[0])
    print(J_tot)
    return flux_array, J_tot


# Solve 2nd order ODE for mass balance in the z-direction; across the width of the membrane, ie perpendicular to the membrane surface.
def component_weight_fraction_derivative(s, c, dwdz, flux):
    # Differential for mass balance, assuming steady-state. Must equal to 0.
    diff = (s.u * dwdz + flux / (s.rho_tot * s.membrane_height)) / c.Dax
    return diff


def mass_differential(s, flux_array, J_tot):
    dmdz = s.membrane_width * J_tot
    return dmdz


def dm_dX(s, flux):
    dmdx = -2 * s.membrane_height * flux
    return dmdx


def odes(z, x, system_data, components):
    # Update system class object.
    system_data.update(new_mass_fractions=[x[i]/sum(x) for i in range(len(x))],
                       new_mass_flow=sum(x),
                       new_T=system_data.temp,
                       new_P=system_data.pressure
                       )

    # Update component class objects.
    for comp in components: comp.update(system_data=system_data)

    # Assign dwdz for x[i] at which i is a weight fraction.
    derivs = np.zeros_like(x)

    # Calculate the component SECOND derivatives.
    flux_array, J_tot = solve_component_flux(s=system_data, components=components)
    for comp in components:
        # Calculate 2nd derivative, assign to dwdz.
        derivs[comp.n] = dm_dX(s=system_data, flux=flux_array[comp.n])

    return derivs


def main():
    # Extract data from templates. Create class objects for system data and component properties.
    system_data, components = parse_templates()

    # Initial weight fractions, and second derivatives are 0.
    x0 = []
    for wt in system_data.mass_fractions:
        x0.append(wt * system_data.total_mass_flow)

    # Solve ODEs.
    sol = solve_ivp(fun=odes, t_span=(0, system_data.length), y0=x0, method='RK45', args=(system_data, components))

    # Extract each component solution.
    df = pd.DataFrame()
    df['z'] = sol.t

    # Plot mass flow of each component, and weight fractions.
    fig = go.Figure()
    fig_1 = go.Figure()
    df_2 = df.drop(columns=['z'])
    df_normalized = df.div(df_2.sum(axis=1), axis=0)

    for comp in components:
        df[comp.label] = sol.y[comp.n]
        # Show as figure.
        fig.add_trace(go.Scatter(x=sol.t, y=sol.y[comp.n], name=comp.label, yaxis='y1', mode='lines'))
        fig_1.add_trace(go.Scatter(x=df['z'], y=df_normalized[comp.label], mode='lines', name=comp.label))

    fig.show()
    fig_1.show()
    print(df)

    return 0


main()