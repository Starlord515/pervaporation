import numpy as np
import pandas as pd
from thermo.nrtl import *

#Todo saturation pressures and viscosities, for Re calculation

class COMPONENT_PROPERTIES:
    def __init__(self, label, q0, T0, psat, A, B, D, rmm, rho, system_data, diffusion_coefficient, viscosity, n):
        # Basics
        self.n = n
        self.label = label

        # Properties
        self.mu = viscosity
        self.psat = psat
        self.rmm = rmm
        self.rho = rho
        self.D = diffusion_coefficient

        # Mass fraction and flow.
        self.fw = system_data.mass_fractions[self.n]
        self.pw = system_data.p_wf[self.n]
        self.initial_mass = self.fw * system_data.total_mass_flow

        # Mole fractions
        self.x = system_data.mole_fractions[self.n]
        self.permeate_y = system_data.p_mole_fractions[self.n]

        # Component stream parameters.
        self.gamma = self.get_gamma(system_data)
        self.Sc = self.calc_Sc()
        self.Sh = self.calc_Sh(system_data)
        self.Dax = self.calc_Dax(system_data)
        self.k = self.calc_k(self.Sh, system_data.dh)

        # Flux model parameters.
        self.q0 = q0
        self.T0 = T0
        self.A = A
        self.B = B
        self.D = D
        self.Q = self.calc_Q(system_data.water_wt, system_data.temp)

    def update(self, system_data):
        # Mass fraction and flow.
        self.fw = system_data.mass_fractions[self.n]
        self.pw = system_data.p_wf[self.n]
        self.initial_mass = self.fw * system_data.total_mass_flow

        # Mole fractions.
        self.x = system_data.mole_fractions[self.n]
        self.permeate_y = system_data.p_mole_fractions[self.n]

        # Recalculate parameters.
        self.gamma = self.get_gamma(system_data)
        self.psat = self.get_psat(system_data.temp)
        self.Sc = self.calc_Sc()
        self.Sh = self.calc_Sh(system_data)
        self.calc_Dax(system_data)
        self.k = self.calc_k(self.Sh, system_data.dh)
        self.Q = self.calc_Q(system_data.water_wt, system_data.temp)
        return 0

    def get_gamma(self, system_data):
        try:
            gamma = system_data.ge.gammas()
            gamma = gamma[self.n]
        except:
            gamma = self.gamma
        return gamma

    def get_psat(self, T):
        # Todo: use model to calculate saturation pressure for pure component.
        return self.psat

    def calc_Q(self, water_wt_frac, T):
        return self.q0 * water_wt_frac ** self.A * np.exp(-self.B/8.31 * (1/self.T0 - 1/T))

    def calc_k(self, Sh, dh):
        return Sh * self.D / dh

    def calc_permeate_mole_fraction(self, system_data):
        y = (self.initial_mass - self.pw * system_data.mass_flow) / (system_data.initial_mol_flow - system_data.mole_flow)

    def calc_Sc(self):
        return self.mu / (self.rho * self.D)

    def calc_Sh(self, system_data):
        return system_data.an[0] * system_data.Re ** system_data.an[1] * self.Sc ** system_data.an[2] * (system_data.dh / system_data.l) ** system_data.an[3]

    def calc_Dax(self, system_data):
        return system_data.u * system_data.dh * (1/(system_data.Re*self.Sc) + system_data.Re*self.Sc/192)


class SYSTEM_PROPERTIES:
    def __init__(self, membrane_width, membrane_height, mass_flow, mass_fractions, rho_i, an, l, permeate_pressure, GE, T, P, rmm_i, length):
        # Membrane measurements.
        self.membrane_width = membrane_width
        self.membrane_height = membrane_height
        self.l = l
        self.dh = self.calc_dh()
        self.length = length

        # Stream temperature and pressure.
        self.temp = T
        self.pressure = P
        self.permeate_pressure = permeate_pressure

        # Component properties.
        self.rho_i = rho_i
        self.rmm_i = rmm_i
        self.N_components = len(mass_fractions)

        # Mass flow, component flow, fractions.
        self.mass_flow = mass_flow
        self.total_mass_flow = mass_flow
        self.initial_component_mass = [mass_fractions[i] * mass_flow for i in range(len(mass_fractions))]
        self.mass_fractions = mass_fractions
        self.water_wt = mass_fractions[0]

        # Mole flows, component flow, fractions.
        self.mole_flow = self.calc_mole_flow(rmm_i, mass_fractions, mass_flow)
        self.initial_mole_flow_i = self.calc_mole_flow_components(rmm_i, mass_fractions, mass_flow)
        self.mole_fractions = wt_to_x(self.mass_fractions, self.rmm_i)

        # Permeate fractions.
        self.p_mass_flow = 0 #self.total_mass_flow - self.mass_flow
        self.p_mass_flow_i = [0, 0] #[self.initial_component_mass[i] - self.mass_flow * self.mass_fractions[i] for i in range(len(mass_fractions))]
        self.p_wf = [0, 0] #[i / self.p_mass_flow for i in self.p_mass_flow_i]
        self.p_mole_flow = [self.p_mass_flow_i[i] / self.rmm_i[i] for i in range(len(mass_fractions))]
        self.p_mole_fractions = wt_to_x(self.p_wf, self.rmm_i)

        # Mixture properties.
        self.rho_tot = self.calc_rho_tot(mass_fractions, rho_i)
        self.u = self.calc_u()
        self.Re = self.calc_Re()
        self.an = an

        # NRTL method for activity coefficients (gamma).
        self.ge = GE

    def update(self, new_mass_fractions, new_mass_flow, new_T, new_P):
        # Update mass flows and fractions.
        #self.mass_flow = new_mass_flow
        self.mass_fractions = new_mass_fractions
        self.water_wt = new_mass_fractions[0]

        # Update mole flows and fractions.
        self.mole_flow = self.calc_mole_flow(self.rmm_i, self.mass_fractions, self.mass_flow)
        self.mole_fractions = wt_to_x(self.mass_fractions, self.rmm_i)

        # Update temperatures and pressures.
        self.temp = new_T
        self.pressure = new_P

        # Recalculate stream properties with new values.
        self.rho_tot = self.calc_rho_tot(new_mass_fractions, self.rho_i)
        self.mass_flow = new_mass_flow
        self.u = self.calc_u()
        self.Re = self.calc_Re()
        self.ge = self.ge.to_T_xs(new_T, self.mole_fractions)

        # Recalculate permeate properties.
        self.p_mass_flow = self.total_mass_flow - self.mass_flow
        if self.p_mass_flow != 0:
            self.p_mass_flow_i = [self.initial_component_mass[i] - self.mass_flow * self.mass_fractions[i] for i in
                                  range(len(self.mass_fractions))]
            self.p_wf = [i / self.p_mass_flow for i in self.p_mass_flow_i]
            self.p_mole_flow = [self.p_mass_flow_i[i] / self.rmm_i[i] for i in range(len(self.mass_fractions))]
            self.p_mole_fractions = wt_to_x(self.p_wf, self.rmm_i)
        else:
            self.p_mass_flow_i = [0 for i in range(self.N_components)] #[self.initial_component_mass[i] - self.mass_flow * self.mass_fractions[i] for i in
                                  #range(len(self.mass_fractions))]
            self.p_wf = self.p_mass_flow_i # [i / self.p_mass_flow for i in self.p_mass_flow_i]
            self.p_mole_flow = self.p_mass_flow_i # [self.p_mass_flow_i[i] / self.rmm_i[i] for i in range(len(self.mass_fractions))]
            self.p_mole_fractions = self.p_mass_flow_i # wt_to_x(self.p_wf, self.rmm_i)


    def calc_mole_fractions(self, rmm_i, mass_fractions, mass_flow):
        mole_flows = [mass_fractions[i] * mass_flow / rmm_i[i] for i in range(len(mass_fractions))]
        sum_up = sum(mole_flows)
        return [mole_flows[i] / sum_up for i in range(len(mole_flows))]

    def calc_mole_flow(self, rmm_i, mass_fractions, mass_flow):
        mole_flows = [mass_fractions[i] * mass_flow / rmm_i[i] for i in range(len(mass_fractions))]
        return sum(mole_flows)

    def calc_mole_flow_components(self, rmm_i, mass_fractions, mass_flow):
        mole_flows = [mass_fractions[i] * mass_flow / rmm_i[i] for i in range(len(mass_fractions))]
        return mole_flows

    def calc_rho_tot(self, mass_fractions, rho_i):
        sum_up = 0
        for i in range(len(self.mass_fractions)):
            sum_up += mass_fractions[i] * rho_i[i]
        return sum_up

    def calc_u(self):
        return self.mass_flow / self.rho_tot / (self.membrane_width * self.membrane_height)

    def calc_dh(self):
        return 2 * self.membrane_width * self.membrane_height / (self.membrane_width + self.membrane_height)

    def calc_Re(self):
        return self.rho_tot * self.u * self.dh / (10**-5) #self.mu Todo: find estimation method for viscosity.



def parse_templates():
    # Parse system data from template.
    data_1 = pd.read_excel('/home/ben/Documents/Y3_dissertation/membrane_model/membrane_module_data_template.xlsx')
    data_1 = data_1.drop(columns=['Property', 'Units', 'Description'])
    sys_prop = data_1['Value'].values

    # Parse component property data.
    data_2 = pd.read_excel('/home/ben/Documents/Y3_dissertation/membrane_model/membrane_simulation_template.xlsx')
    data_2 = data_2.drop(columns=['Property', 'Units', 'Description'])

    # Get gamma operator. Values are for an ethanol-water mixture at standard temperature.
    from scipy.constants import calorie, R
    N = 2
    T = 70.0 + 273.15
    xs = [0.252, 0.748]
    tausA = tausE = tausF = tausG = tausH = alphaD = [[0.0] * N for i in range(N)]
    tausB = [[0, -121.2691 / R * calorie], [1337.8574 / R * calorie, 0]]
    alphaC = [[0, 0.2974], [.2974, 0]]
    ABEFGHCD = (tausA, tausB, tausE, tausF, tausG, tausH, alphaC, alphaD)
    GE = NRTL(T=T, xs=xs, ABEFGHCD=ABEFGHCD)

    # Update operator to initial temperature and composition.
    weights = data_2.iloc[9].to_numpy()
    rmms = data_2.iloc[5].to_numpy()
    xs = wt_to_x(weights, rmms)
    GE = GE.to_T_xs(sys_prop[5], xs)

    # Place system values into class.
    system_data = SYSTEM_PROPERTIES(membrane_width = sys_prop[0],
                                    membrane_height = sys_prop[1],
                                    mass_flow = sys_prop[2],
                                    mass_fractions = data_2.iloc[9].to_numpy(),
                                    rho_i = data_2.iloc[6].to_numpy(),
                                    an = sys_prop[7:11],
                                    l = sys_prop[3],
                                    permeate_pressure = sys_prop[4],
                                    GE = GE,
                                    T = sys_prop[5],
                                    P = sys_prop[6],
                                    rmm_i = rmms,
                                    length = sys_prop[11])

    # Create an array of components, stored as class objects.
    i=0
    components = []
    for species in data_2.columns:
        prop = data_2[species].values
        components.append(COMPONENT_PROPERTIES(label=species,
                                               q0 = prop[0],
                                               T0 = prop[1],
                                               psat=prop[7],
                                               A=prop[2],
                                               B=prop[3],
                                               D=prop[4],
                                               rmm=prop[5],
                                               rho=prop[6],
                                               system_data=system_data,
                                               diffusion_coefficient=prop[4], #Todo: put diffusivity coefficients in template
                                               viscosity=prop[10], #Todo: add viscosities to template
                                               n=i))
        i += 1
    return system_data, components


def wt_to_x(x, rmm):
    if all(i != 0 for i in x):
        x_1 = [x[i] / rmm[i] for i in range(len(x))]
        sum_up = sum(x_1)
        wt = [i / sum_up for i in x_1]
    else:
        wt = [0 for i in range(len(x))]
    return wt