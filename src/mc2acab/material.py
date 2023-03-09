#! /usr/bin/env python

''' Material object for ESS-Bilbao MCNP suite by
    Miguel Magan and Dr. Octavio Gonzalez Feb 2023'''

# import re
import numpy as np
from . import MCNP_outparser
from .pyhtape3x import atomic_mass

def is_number(string):
    '''Check if a number is a number'''
    try:
        float(string)
        return True
    except ValueError:
        return False

class Mat:
    """
    This is an MCNP material.
    """
    def __init__ (self,number):
        self.number = number
        self.zaid = []  # Isotopes
        self.frac = []  # Abundance

    def n2ro(self, ro):
        """ Rearrange the isotopic composition to match atomic density in at/bn-cm,
        for ro material density"""
        in_ro = sum(self.frac)
        if in_ro > 0:  # isotope composition in fractions, we only need a norm. factor
            if ro < 0:  # Density is in g/cm3, and we need to put it in atomic dens
                atomro = 0
                for i, atom in enumerate(self.zaid):
                    (z,a) = divmod(atom,1E3) #  a is the atomic mass z is the atomic number
                    if a != 0:
                        atomro += a * self.frac[i]
                    else:
                        atomro += atomic_mass(z) * self.frac[i]
                factor = -ro / (atomro/0.6023)
            else:
                factor = ro / in_ro
            self.frac[:]=[i*factor for i in self.frac]
        else: # Isotope composition in weight
            if ro<0: # Density in g/cm3.
                self.frac[:]=[i/in_ro for i in self.frac]  # Now we have weight fraction
                for index,atom in enumerate(self.zaid):
                    (z,a)=divmod(atom,1000)
                    iso_ro=-self.frac[index]*ro # Density of this isotope
                    if a!=0: # Specific isotope
                        self.frac[index]=iso_ro/(a/0.6023)  # This is the atomic density of the isotope
                    else: # Natural abundance
                        self.frac[index]=iso_ro/(atomic_mass(z)/0.6023)  # This is the atomic density of the isotope
            else: # Density is in at/bn-cm
                for index,atom in enumerate(self.zaid):
                    (z,a)=divmod(atom,1000)
                    if a!=0: # Specific  isotope
                        self.frac[index]=self.frac[index]/a
                    else: #Natural abundance
                        self.frac[index]=self.frac[index]/atomic_mass(z)
                # Now we have a proportional atomic composition. Now it is just apply a factor
                in_ro=sum(self.frac) # Recalculate now the input density
                self.frac=[i*ro/in_ro for i in self.frac]

    def normalize(self):
        """Normalize the fractions to 1, either positive or negative"""
        totalfract = sum(self.frac)
        self.frac = [m/totalfract*np.sign(totalfract) for m in self.frac]

    def __eq__(self,other):
        ''' __eq__ overload to compare isotopes and compositions '''
        self.normalize()
        other.normalize()
        return all([self.zaid == other.zaid, self.frac == other.frac])

# ====================================================== #

def oget(infile, number):
    """ Get the material number from MCNP output file infile using material declaration"""
    if number == 0:
        return Mat(0)
    N = []
    M = []
    inputlines = MCNP_outparser.input_finder(infile)
    m_info = []
    for i, line in enumerate(inputlines):
        if line[0] in ['m','M'] and is_number(line.split()[0][1:]):
            if int(line.split()[0][1:]) == number:
                print (f"found material {number}")
                tokens = MCNP_outparser.line_parser(line)
                m_info.extend(tokens[1:])
                init_M = i+1
                break
    else:
        print(f"material {number} not found")
        return None
    line = inputlines[init_M]
    while line[0] in ["c", "C", " "]:
        tokens = MCNP_outparser.line_parser(line)
        m_info.extend(tokens)
        init_M += 1
        line = inputlines[init_M]
    material = Mat(number)
    m_info = [m for m in m_info if m != '']
    for i in range(0, len(m_info), 2):
        N.append(int(m_info[i].split(".")[0]))
        M.append(float(m_info[i+1]))
    material.zaid = N
    material.frac = M
    return material

# TODO: This should be in a separate library

# def mgetall(infile):
#     """ Get all materials info from MCNP output infile using the cells table"""
#     lcell = np.dtype([('cellNum',int),('Material',object),('RoA',float),('RoG',float),('cellID',int),('volume',float)])
#     materials_map = np.zeros((1),dtype=lcell)
#     mat_index = []
#     with open(infile,"r", encoding='utf-8') as outp:
#         lines = outp.readlines()
#         for i, line in enumerate(lines):
#             if re.search('table 60', line) and re.match('1cells', line):
#                 words = lines[i+5].split()
#                 while len(words) != 0:
#                     cID = int(float(words[0])) # Program ID of the cell.
#                     cN = int(float(words[1])) # cell number
#                     vol = float(words[5]) # cell volume
#                     mat_i = int((words[2].replace('s',''))) # material number
#                     if mat_i not in mat_index:
#                         print(f"adding material: {mat_i}\n")
#                         mat_index.append(mat_i)
#                     mater = oget(infile, mat_i)
#                     densidad_at = float(words[3]) # atom density
#                     densidad_gr = float(words[4]) # gram density
#                     if mater is not None:
#                         mater.n2ro(densidad_at)
#                     linea = np.array((cN,mater,densidad_at,densidad_gr,cID,vol),dtype=lcell)
#                     materials_map = np.append(materials_map,linea)
#                     i += 1
#                     words = lines[i+5].split()
#                 break
# #    print(len(mat_index))
#     return materials_map,len(materials_map),len(mat_index)
