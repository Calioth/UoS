# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate, optimize
from math import sqrt, atan, sin, tan, cos
import csv

'''class IFT:
    def __init__(self, Mach, pressureratio, densityratio, temperatureratio, flowangle, arearatio):
        self.Mach = Mach
        self.pressureratio = pressureratio
        self.densityratio = densityratio
        self.temperatureratio = temperatureratio
        self.flowangle = flowangle
        self.arearatio = arearatio
        
    def __str__(self):
        return f"{self.Mach}, {self.pressureratio}, {self.densityratio}, {self.temperatureratio}, {self.flowangle}, {self.arearatio}"'''
        
M = np.arange(0, 6, 0.02)

def pressure_ratio(M, gamma = 1.4):
    return round((1+((gamma-1)*M**2)/2)**(-gamma/(gamma-1)),4)

def density_ratio(M, gamma = 1.4):
    return round((1+((gamma-1)*M**2)/2)**(-1/(gamma-1)),4)

def temperature_ratio(M, gamma = 1.4):
    return round((1+((gamma-1)*M**2)/2)**-1,4)

def area_ratio(M, gamma =1.4):
    if M == 0:
        return 0
    return round(1/M * ((2/(gamma+1))*(1+((gamma-1)*M**2)/2))**((gamma+1)/(2*(gamma-1))),4)

def flow_angle(M, gamma = 1.4):
    if M < 1:
        return 0
    x = np.linspace(1., M, 10000)
    y = np.sqrt(x**2 - 1)/(x*(1+((gamma-1)*x**2)/2))
    return round(np.degrees(np.trapz(y, x)), 4)

def mach_after_shock(M, gamma = 1.4):
    return round(np.sqrt((2+(gamma-1)*M**2)/(2*gamma*M**2-(gamma-1))), 4)

def p_across_shock(M, gamma = 1.4):
    return round(1 + (2*1.4*(M**2 - 1))/(gamma + 1), 4)

def rho_across_shock(M, gamma = 1.4):
    return round(((gamma + 1)*M**2)/(2 + (gamma - 1)*M**2), 4)

def t_across_shock(M, gamma = 1.4):
    x = p_across_shock(M)
    y = rho_across_shock(M)
    return round(x/y, 4)

def stagp_across_shock(M, gamma = 1.4):
    M2 = mach_after_shock(M)
    x = 1/ pressure_ratio(M2)    
    y = pressure_ratio(M)
    z = p_across_shock(M)
    return round(x*y*z, 4)

def stagp2_p1(M, gamma = 1.4):
    M2 = mach_after_shock(M)
    return round(p_across_shock(M)*(1/pressure_ratio(M2)), 4)

def nst_writer():
    with open("NST.csv", 'w', newline='' ) as file:
        writer = csv.DictWriter(file, fieldnames=['Mn1','Mn2','p2/p1','rho2/rho1', 'T2/T1', 'p02/p01', 'p02/p1'])
        writer.writeheader()
        for i in M:
            if i < 1:
                continue
            else:
                writer.writerow({'Mn1': round(i,4) , 'Mn2': mach_after_shock(i, 1.667), 'p2/p1': p_across_shock(i, 1.667), 'rho2/rho1': rho_across_shock(i, 1.667), 'T2/T1': t_across_shock(i, 1.667), 'p02/p01': stagp_across_shock(i, 1.667), 'p02/p1': stagp2_p1(i, 1.667)})

def main():
    x = input("Which table do you want?").strip()  
    match x:
        case 'IFT':
            ift_writer()
        case 'NST':
            nst_writer()
    
def ift_writer():
    with open("IFT.csv", 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Mach', 'p/p0', 'rho/rho0', 'T/T0', 'v', 'A/A*'])
        writer.writeheader()
        for i in M:
            Mach = round(i,4)
            pressureratio = pressure_ratio(i)
            densityratio = density_ratio(i)
            temperatureratio = temperature_ratio(i)
            flowangle = flow_angle(i)
            arearatio = area_ratio(i)
            writer.writerow({'Mach': Mach, 'p/p0': pressureratio, 'rho/rho0': densityratio, 'T/T0': temperatureratio, 'v': flowangle, 'A/A*': arearatio})

def oblique_shock_relation(M, beta, gamma = 1.4):
# Equation:
    beta = np.radians(beta)
    tan_theta = 2 * (1 / tan(beta)) * (
    (M**2 * sin(beta)**2 - 1) /
    (M**2 * (gamma + cos(2 * beta)) + 2))
# If you want θ itself:
    theta = np.degrees(atan(tan_theta))
    return round(theta,4)



def interpolate(given_val, given_ser, result_ser):
    xp=[]
    fp=[]
    with open('IFT.csv') as file:
        reader = csv.DictReader(file)
        for row in reader:
            xp.append(row[given_ser])
            fp.append(row[result_ser])
        xp = np.asarray(xp, dtype=float)
        fp = np.asarray(fp, dtype=float)
        return np.interp(given_val, xp, fp)
         
            
        
if __name__ == "__main__":
    main()