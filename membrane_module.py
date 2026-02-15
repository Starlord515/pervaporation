#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 14 12:54:19 2026

@author: ben

Models the output mass fractions, temperature, and pressure.
"""

import numpy as np


def component_membrane_flux(wfm, T, pp, memflux_args):
    # memflux_args = [q0, ai, B, psat, gamma_L0, T0]
    q = memflux_args[0]
    ai = memflux_args[1]
    B = memflux_args[2]
    psat = memflux_args[3]
    gamma_L0 = memflux_args[4]
    T0 = memflux_args[5]
    
    # Where ww is the weight-fraction of water, ai is the permeance coefficient
    # of i, B is another permeance coefficient of i.
    
    xL0 = wfm * mtot / RMi
    yl = wp * pm / RMi
    
    flux = q * (psat * xL0 * gamma_L0 - pp * yl)
    # Where xL0 is the mole-fraction of the liquid on the membrane surface,
    # and yl is the permeate mole-fraction.
    # psat is the sat vapor pressure of i, gamma_L0 is the liquid activity
    # coefficient of i.
    return flux


def boundary_layer_flux(wfm, wp, wf, boundary_args):
    # boundary_args = [dh, D, rho_tot, k]
    dh = boundary_args[0]
    D = boundary_args[1]
    rho_tot = boundary_args[2]
    k = boundary_args[3]
    
    flux = rho_tot * k * np.log((wfm - wp)/(wf - wp))
    # Where k is the mass-transfer coefficient.
    return flux


# Numerically solve for membrane flux at a point along the length of the module.
# Need to match yl (permeate molar fraction) and weight 
from scipy.optimize import minimize
def solve_component_flux(water_w, T, pp, wp, wf, an, D, B, mem_width, channel_h, rho_tot, Re, Sc, q0, T0, psat, gamma_L0, dh):
    def objective_flux_fun(wfm, T, pp, memflux_args, wp, wf, boundary_args):
        f1 = component_membrane_flux(wfm, T, pp, memflux_args)
        f2 = boundary_layer_flux(wfm, wp, wf, boundary_args)
        return f1 - f2

    # Estimate mass-transfer coefficient
    a1 = an[0]; a2=an[1]; a3=an[2]; a4=an[3]
    Sh = a1*(Re**a2)*(Sc**a3)*((dh/l)**a4)
    k = Sh * D / dh
    # Permenace Q
    q = q0 * water_w**ai * np.exp( - B / R *(1/T0 - 1/T))
        
    # Flux function argument arrays
    boundary_args = [dh, D, rho_tot, k]
    memflux_args = [q, ai, B, psat, gamma_L0, T0]
    
    bounds = tuple((0, 1))
    guess =  [0.9 * wf]
    wfm_solution = minimize(objective_function, guess, args=(water_w, T, pp, memflux_args, wp, wf, boundary_args), method='SLSQP', bounds=bounds)
    component_flux_solution = component_membrane_flux(wfm_solution.x, water_w, T, pp, memflux_args)
    return component_flux_solution


# Solve 2nd order ODE for mass balance, given w at z=0, and dwdz at z=0.
def channel_mass_balance(u, water_w, T, pp, wp, wf, an, D, B, mem_width, channel_h, rho_tot, Re, Sc, q0, T0, psat, gamma_L0, dh, Dax):
    # Calculate flux
    flux = solve_component_flux(water_w, T, pp, wp, wf, an, D, B, mem_width, channel_h, rho_tot, Re, Sc, q0, T0, psat, gamma_L0, dh)
    
    # Differential for mass balance, assuming steady-state. Must equal to 0.
    diff = ( u * dwdz + flux / (rho_tot * channel_h) ) / Dax
    # Where Dax is axial dispersion, u is flow velocity, 
    return diff


def permiate_mass_differential():
    dmdz = dh * flux_array.sum() * m_tot
    return dmdz


def odes(x, z):
    # Assign each weight fraction ODE to an element of x.
    w_ethanol = x[0]
    dwdt_ethanol = x[1]
    w_water = x[2]
    dwdt_water = x[3]
    w_el = x[4]
    dwdt_el = x[5]
    permeate_mass_flow = x[6]
    
    # Constants
    T0 = 300 #K, relative temperature for permeance
    q0 = [] # permeance
    Ai = []
    Bi = []
    psat = [] # Pa, saturation pressure of mixture
    gamma = [] # Activity coefficients
    RMi = [] # kg/kmol, molecular masses
    Di = [] # component diffusivities
    rho = [] # kg/m3, component densities
    a = [] # Sc coefficients used to estimate mass-transfer coefficient for membrane
    total_mass_flow = 1 # kg/s
    
    mem_width = 1 # m, width of the membrane used
    channel_h = 1 # m. height of a channel
    
    pp = 1 #Pa, permeate-side pressure
    T = 330 # K, feed-side temperature
    
    # Calculate needed parameters
    Re = 
    Sc = 
    rho_tot = w_ethanol*rho[0] + w_water*rho[1] + w_el*rho[2]
    u = (total_mass_flow / rho_tot) / (mem_width * channel_h)
    # Axial dispersion coefficient
    Dax = u * dh * (1/(Re*Sc) + Re*Sc/192)
    # Estimate hydraulic diameter
    dh = 2 * (mem_width  * channel_h) / (mem_width + channel_h)
    
    # Calculate the derivatives.
    2nder_ethanol = channel_mass_balance(u, water_w, T, pp, wp, wf, an, D, B, mem_width, channel_h, rho_tot, Re, Sc, q0, T0, psat, gamma_L0, dh, Dax)

    
    
def main(w, Tf, pp, fflow):
    # Determine mass balance across module.
    mass_balance()
    
    return 0