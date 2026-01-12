#! /usr/bin/env python

"""
Allows to create from scracht an input for ACAB 2008, 
Requires:
- A txt file with a material definition with format  MCNP  material card
- A txt file with a flux in 2 columns with 211 energy groups corresponding to Vit-J 
- A volume in cm3 for de material cell be required for execution
- A irratiation time in hours that would be required for execution
"""

import os
import numpy as np
import mc2acab.MCNP_ACAB_library as MCNPACAB
import mc2acab.cell as cel
import mc2acab.material as material
import mc2acab.MCNP_outparser as MCNP_outparser
import tally as tal

def is_number(string):
    '''Check if a number is a number'''
    try:
        float(string)
        return True
    except ValueError:
        return False

def get_user_input(flux_file,mat_file):
    """
    Prompts the user to input values for source term and irradiation time,
    read flux file (3 columns, energy_bin/value/error) and returns a tally object
    irradiation time in seconds, and source term value in particles per second.
    """
    # Check if the flux file exists
    if not os.path.exists(flux_file):
        raise FileNotFoundError('flux file not found')
    if not os.path.exists(mat_file):
        raise FileNotFoundError('mat file not found')

    source = MCNPACAB.get_user_source()
    irr_time = MCNPACAB.get_user_time()
    user_input = input('Input cell volume in cm3: ')
    cell_vol = float(user_input)
    user_input = input('Input material number: ')
    mat_ID = int(user_input)
    user_input = input('Input material density in g/cm3: ')
    mat_ro = abs(float(user_input))*-1
    Nuc_lib, id_Egroup = MCNPACAB.get_user_LIB()

    flux_data = np.zeros([212,3], dtype=float)
    with open(flux_file, 'r') as file:
        i = 0
        for linea in file:
            if linea[:2] != "c ":
                for j, val in enumerate(linea.split()):
                    if is_number(val):
                        flux_data[i,j] = val
                i += 1
    
    print('Obtained flux')
    tally_i = tal.Tally(1,1,1,1,211)
    tally_i.comment = 'Fake tally for MCNP_ACAB'
    tally_i.ebins = flux_data[:]
    tally_i.value[0,0,0,0,:] = flux_data[:,1]
    tally_i.error[0,0,0,0,:] = flux_data[:,2]
    tally_i.mass[0,0] = cell_vol
    
    mat_i = get_mat(mat_ID, mat_file)
    
    cell_i = cel.Cell(f'ACAB_m{mat_ID}_flux{source:.1e}part_irrT{irr_time:.1e}s')
    cell_i.density = mat_ro
    cell_i.volume = cell_vol
    print('Obtained cell')
    
    return tally_i, cell_i, mat_i, irr_time, source, Nuc_lib, id_Egroup

def get_mat(mat_ID, mat_file):
    """
    Reads a MCNP format material in a txt file
    """
    N = []
    M = []
    m_info = []
    with open(mat_file, 'r') as file:
        for i, line in enumerate(file):
            if line[0] in ['m','M'] and is_number(line.split()[0][1:]):
                if int(line.split()[0][1:]) == mat_ID:
                    print (f"found material {mat_ID}")
                    tokens = MCNP_outparser.line_parser(line)
                    m_info.extend(tokens[1:])
                    init_M = i+1
                    break
        else:
            print(f"material {mat_ID} not found")
            return None
        file.seek(0)
        inputlines = file.readlines()
    line = inputlines[init_M]
    while line[0] in ["c", "C", " "]:
            tokens = MCNP_outparser.line_parser(line)
            m_info.extend(tokens)
            init_M += 1
            try:
                line = inputlines[init_M]
            except:
                break
    mat_i = material.Mat(mat_ID)
    m_info = [m for m in m_info if m != '']
    for i in range(0, len(m_info), 2):
        N.append(int(m_info[i].split(".")[0]))
        M.append(float(m_info[i+1]))
    mat_i.zaid = N
    mat_i.frac = M
    print('Obtained material')
    return mat_i    


def exec_MCNP_ACAB(flux_file, mat_file, part_i='n', sce_file_i = None, save_i=True):
    """
    Executes MC2ACAB without MCNP files, basically writes down ACAB input and executes all steps to get ACAB results
    """
    tally_i, cell_i, mat_i, irr_time, source_i, nuc_lib_i, id_Egroup_i = get_user_input(flux_file,mat_file)
    outputs = MCNPACAB.MCNP_ACAB_Map(tally0 = tally_i, mater = mat_i, n_id = 0,
                                    irr_cell = cell_i, irr_time = irr_time,
                                    irr_type = part_i, source = source_i,
                                    save = save_i, esc_file = sce_file_i,
                                    passive_sector = None,
                                    id_lib = nuc_lib_i,
                                    id_ILIB = id_Egroup_i,
                                    corte = 0.9)
    # Outputs:
    #     0 Timesets (arrays of times)
    #     1 Decay= Bq as ACAB
    #     2 Gamma= PHOTONS/CCM/SEC (as ACAB)
    #     3 Heat= W/cm3
    #     4 Dose= mSv/h (ACAB is Sv/h)
    #     5 mol = mol
    return outputs
