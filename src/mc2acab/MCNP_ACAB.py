#! /usr/bin/env python

# Version Feb 2023
# LECTURA DE LA SALIDA DE HTAPE3X.
# GENERACION DE INPUT DE ACAB.
# Reescrito en Python para salvaguardar la salud mental de todos los usuarios.
# By Miguel Magan and Octavio Gonzalez

import sys
import os
import MCNP_ACAB_library as MCNPACAB
import material
from multiprocessing import Pool
import cell as cel
import numpy as np
import tally as tal

def MCNP_ACAB_Mapstar(n):
    outputs = MCNPACAB.MCNP_ACAB_Map(tally0 = tally0, mater = mat[n], n_id = n,
                                     irr_cell = irr_cell[n], irr_time = reqs['-irr_time'],
                                     irr_type = reqs['-part'], source = reqs['-st'],
                                     save = options['-save'], esc_file = options['-sce_file'], 
                                     passive_sector = options['-passive_sector'],
                                     id_lib =options['-nuc_lib'], 
                                     id_ILIB = options['-id_Egroup'],
                                     corte = options['-apypa_verge'])
    # Outputs:
    #     0 Timesets (arrays of times)
    #     1 Decay= Bq as ACAB
    #     2 Gamma= PHOTONS/CCM/SEC (as ACAB)
    #     3 Heat= W/cm3
    #     4 Dose= mSv/h (ACAB is Sv/h)
    #     5 mol = mol
    return outputs

print('''
      ***************************************************************************
      *  \033[36m SPALLATION PRODUCT ACTIVATION \033[0m                                        *
      *  This subroutine reads the outp file generated by MCNPX,                *
      *  and the spallation products recorded in the histp file.                *
      *  It generates an input for COLLAPS and then for ACAB, including the     *
      *  spallation products.                                                   *
      *  It executes and post-processes the output of ACAB.                     *
      *  MIGUEL MAGAN ROMERO,Dr. OCTAVIO GONZALEZ                               *
      * \033[36m Targets and Neutronic Applications Group \033[0m                              *
      *                                             \033[31m ESS-BILBAO \033[0m                *
      ***************************************************************************''')

if len(sys.argv) <= 1:
    print ('MCNP2ACAB Arguments [Opcions]')
    print ('\033[31m Mandatory: \033[0m')
    print ('-n    Neutron activation (just n, no histp module required)')
    print ('-p    Calculo de activacion (Modulos histp y flujo protonico)')
    print ('-np   Calculo de activacion (Flujo neutronico y modulo histp)')
    print ('\033[36m Requirements \033[0m')
    print ('-outpfile=filename outp file name defatult outp')
    print ('-tally_num = number for specific tallies')
    print ('-st_units=number Indicates the how indicates source term: 1 particles/s or 2 mA')
    print ('-source_term=number Indicates source term value ')
    print ('-irr_time=number Indicates irradiation time in hours ')
    print ('\033[36m Options \033[0m')
    print('-normal_flux Use normalized flux with FM card')
    print('-sce_file=File Use external irradiation scenario file')
    print('-save=[All,True,False] Modify ACAB_writer output files'
          ' and delete folders after execution')
    print('-threshold=[0 to 1] Indicates cutoff for Apipa')
    print('-decay_times=Set decay times list for ACAB (no spaces)')
    print('-Rotate=n Use cell composition with passive cells'
          ' terminated in n. Does not work for rotary elements')
    print ('')
    sys.exit(1)

# '''Para hacer una ejecucion automatizable'''
reqs ={
    '-part': None,
    '-outpfile': None,
    '-tally_num': None,
    '-st_units': None,
    '-st': None,
    '-irr_time': None,
       }
options = {
    '-normal_flux': False, # useless??
    '-save': True,
    '-sce_file': None,
    '-apypa_verge' : 0.9,
    '-decay_times' : None,
    '-decay_outs' : None,
    '-passive_sector' : None, # so far useless
    '-nuc_lib': 'EAF', 
    '-id_Egroup': 'vitJ+', 
}
def __parse_args(reqs, options, args):
    for arg in args:
        if  arg == '-np' and reqs['-part'] == None:
            reqs['-part'] = 'np'
        elif  arg == '-n' and reqs['-part'] == None:
            reqs['-part'] = 'n'
        elif arg =='-p' and reqs['-part'] == None:
            reqs['-part'] = 'p'
        elif arg.startswith('-outpfile='):
            reqs['-outpfile'] = arg.split('=')[1]
        elif arg.startswith('-tally_num='):
            reqs['-tally_num'] = arg.split('=')[1]
        elif arg.startswith('-st_units='):
            reqs['-st_units'] = arg.split('=')[1]
            while reqs['-st_units'] not in ['1','2']:
                reqs['-st_units'] = input('Warning!!! just 1 or 2:')
        elif arg.startswith('-source_term='):
            reqs['-st'] = float(arg.split('=')[1])
        elif arg.startswith('-irr_time='):
            reqs['-irr_time'] = float(arg.split('=')[1]) * 3600
        elif '-normal_flux' in arg:
            options['-normal_flux'] = True
        elif arg.startswith('-sce_file='):
            options['-sce_file'] = arg.split('=')[1]
            while not os.path.exists(options['-sce_file']):
                options['-sce_file'] = input('Warning!!! No scenenario file!!!'
                                             f" {options['-sce_file']} is not here"
                                             '\nEnter scenario file name: ')
        elif arg.startswith('-threshold='):
            options['-apypa_verge'] = float(arg.split('=')[1])
        elif arg.startswith('-decay_times='):
            dec_times = arg.split('=')[1]
            options['-decay_times'] = [float(time) for time in dec_times.split(',')]
            options['-decay_outs'] = [0]
            for time in options['-decay_times']:
                if float(time) >= 3600:
                    options['-decay_outs'].append(1)
                else:
                    options['-decay_outs'].append(0)
        # elif arg.startswith('-decay_outs='):
        #     dec_outs = arg.split('=')[1]
        #     options['-decay_outs'] =[int(out) for out in dec_outs.split(',')]
        #     while not all options['-decay_outs'] in [0,1]:
        #         dec_outs = input('-decay outs must be a list of 0 (no output) or 1 (output) and must start with a 0')
        #         options['-decay_outs'] =[int(out) for out in dec_outs.split(',')] 
        elif arg.startswith('-rotate='):
            options['-passive_sector='] = int(arg.split('=')[1])
        elif arg.startswith('-nuc_lib='):
            options['-nuc_lib'] = arg.split('=')[1]
        elif arg.startswith('-id_Egroup='):
            options['-id_Egroup'] = arg.split('=')[1]
    # print(reqs)
    # print(options)
    while reqs['-part'] not in ['n','np','p']:
        print('WARNING!!! Unknown activation process')
        reqs['-part'] = input('Indicate activation process type')
    if os.path.exists('outp') and not reqs['-outpfile']:
        reqs['-outpfile'] = 'outp'
    elif not reqs['-outpfile']:
        reqs['-outpfile'] = input('Enter outp file name: ')
    while not os.path.exists(reqs['-outpfile']):
        reqs['-outpfile'] = input(f"Warning!!! No outp file!!! {reqs['-outpfile']}"
                                  ' is not here \nEnter outp file name: ')
    if reqs['-st_units'] == '2':
        reqs['-st'] *= 6.24E15
    if None in reqs.values():
        raise ValueError('Warning!!! input incomplete, further information required:')
    if not options['-sce_file'] and options['-decay_times']:
        print('Generating automatic scenario file')
        MCNPACAB.Escenary_generator(reqs['-irr_time'], options['-decay_times'],
                                    options['-decay_outs'], Sce_name ='Auto_Sce_file')
        options['-sce_file'] = 'Auto_Sce_file'
    return reqs, options

try:
    reqs, options = __parse_args(reqs, options, sys.argv[1:])
    print(f'source_term = {reqs["-st"]:.3e} n/s')
    input_complete = True
    tally0 = tal.oget(reqs['-outpfile'],reqs['-tally_num'])
except Exception as e:
    print(e)
    tally0, reqs['-irr_time'], reqs['-st'], options['-nuc_lib'], options['-id_Egroup'] = MCNPACAB.get_user_input(reqs['-outpfile'])
#TODO
# tally = MCNPACAB.tally_compose(tally0, Passive_sector)
# cmatrix = MCNPACAB.comp_matrix(tally0, Passive_sector)
ncel = [int(cell0) for cell0 in (tally0.cells)]
print('Obtained cell numbers')
vol0 = tally0.mass
irr_cell = [cel.oget(reqs['-outpfile'],ncell_i) for ncell_i in ncel]
print('Obtained cell properties')
# Because there can be quite a lot of cells with the same material, it is interesting to cache them
mat = []
matnumbers = np.zeros(0)
for ncell0 in irr_cell:
    try:
        Mindex = list(matnumbers).index(ncell0.mat)
        mat0 = material.mat(mat[Mindex].number)
        mat0.N = list(mat[Mindex].N)
        mat0.M = list(mat[Mindex].M)
        mat.append(mat0)
    except:
        mat.append(material.oget(reqs['-outpfile'],ncell0.mat))
        matnumbers=np.append(matnumbers,mat[-1].number)
print('Obtained materials')

with Pool() as pool:
    totaldata = pool.map(MCNP_ACAB_Mapstar, range(tally0.ncells))

if not options['-decay_times']:
    t_times = list(totaldata[0][0].index)
else:
    o_times = [1.0]
    for time in options['-decay_times']:
        if time in list(totaldata[0][0].index) and time not in o_times:
            o_times.append(time)
    t_times = o_times
MCNPACAB.summary_table_gen(totaldata,tally0,t_times)
