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

print('''
    ***************************************************************************
    *              \033[36m ACTIVACION CON PRODUCTOS DE ESPALACION \033[0m                   *
    * Esta subrutina lee el archivo outp generado por MCNPX, y los productos  *
    * de espalacion registrados en el archivo histp.                          *
    * Genera un input para COLLAPS y despues para ACAB, incluyendo los prod.  *
    * de espalacion. Ejecuta y postprocesa la salida de ACAB.                 *
    * MIGUEL MAGAN ROMERO     RAUL VIVANCO SANCHEZ    Dr. OCTAVIO GONZALEZ    *
    *    \033[36m   Grupo de Blancos y Aplicaciones Neutronicas\033[0m                       *
    *                                        \033[31m ESS-BILBAO \033[0m                     *
    ***************************************************************************''')


if len(sys.argv) == 1:
    print ('USO: MCNP_ACAB Argumentos [Opciones]')
    print ('ARGUMENTOS')
    print ('-outpfile=filename outp file name defatul outp')
    print ('-st_units=number Indicates the how indicates source term: 1 particles/s or 2 mA')
    print ('-source_term=number Indicates source term value ')
    print ('-time=number Indicates irradiation time in hours ')
    print ('-tally_num = number for specific tallies')
    print ('-np   Calculo de activacion (Flujo neutronico y modulo histp)')
    print ('-n    Calculo de activacion (SOLO n, sin modulo histp)')
    print ('-p    Calculo de activacion (Modulos histp y flujo protonico)')
    print (' OPCIONES')
    print ('-Normal_flux Usar flujo normalizado con card FM')
    print ('-Rotate=n Usar composicion de celdas con celdas pasivas terminadas'
           ' en n. Para elementos rotatorios No funciona')
    print ('-EXT_FILE=Archivo   Usar archivo externo de escenario de irradiacion')
    print ('-save = [All,True,False] Modifica salida de archivos de ACAB_writer'
           ' y elimina las carpetas tras la ejecución')
    print ('-corte Indica el corte para Apipa')
    print ('-decay_times = Set decay times list for ACAB (no spaces)')
    print ('')
    sys.exit(1)

#print(sys.argv)

esc_file0 = None
save0 = True
passive_sector0 = None
t_times = None
corte0 = 0.9
for arg in sys.argv:
    if arg == '-Normal_flux':
        Normal_flux = True
    if '-EXT_FILE' in arg:
        esc_file0 = arg.split('=')[1]
        print("using external scenario file")
    if  arg == '-np':
        particles='np'
    elif  arg == '-n':
        particles = 'n'
    elif arg =='-p':
        particles = 'p'
    if '-save' in arg:
        if len(arg.split('='))>=1:
            save0 = arg.split('=')[1]
        else:
            save0 = True
    if '-Rotate' in arg:
        str_passive_sector = arg.split('=')[1]
        passive_sector0 = int(str_passive_sector)
        print("Using composed results for rotating elements")
    if '-corte' in arg:
        corte0 = float(arg.split('=')[-1])

if particles not in ['n','np','p']:
    print('WARNING!!! type of activation not included')
    sys.exit(1)

# '''Para hacer una ejecucion automatizable'''
input_complete = False
for i, arg_i in enumerate(sys.argv):
    if '-st_units' in arg_i:
        for j, arg_j in enumerate(sys.argv):
            if '-source_term' in arg_j:
                for k, arg_k in enumerate(sys.argv):
                    if '-time' in arg_k or '-EXT_FILE' in arg_k:
                        for l, arg_l in enumerate(sys.argv):
                            if '-outpfile' in arg_l:
                                for m, arg_m in enumerate(sys.argv):
                                    if '-tally_num' in arg_m:
                                        outpfile = str(arg_l.split('=')[-1])
                                        while not os.path.exists(str(outpfile)):
                                            print('Warning!!! No outp file!!! '
                                                  f'{str(outpfile)} is not here')
                                            outpfile = input('Enter outp file name: ')
                                        ntal = int(arg_m.split('=')[-1])
                                        tally0 = tal.oget(outpfile,ntal)
                                        if '-EXT_FILE' not in arg_k:
                                            irr_time0 = float(arg_k.split('=')[-1])*3600
                                        else:
                                            irr_time0 = 1 #its a dummy number...
                                        source0 = float(arg_j.split('=')[-1])
                                        units = arg_i.split('=')[-1]
                                        if units == '2':
                                            source0 *= 6.24E15
                                        print(f'Termino fuente = {source0:.2e} n/s')
                                        nuc_lib = 'EAF' # the only one that works in ACAB
                                        id_Egroup = 'vitJ+' # the only one that works in ACAB
                                        input_complete = True

if input_complete == False:
    print('Warning!!! input required! Complete further info: ')
    if os.path.exists('outp'):
        outpfile = 'outp'
    else:
        outpfile = input('Enter outp file name: ')
    while not os.path.exists(str(outpfile)):
        print(f'Warning!!! No outp file!!! {str(outpfile)} is not here \n')
        outpfile = input('Enter outp file name: ')
    tally0, irr_time0, source0, nuc_lib, id_Egroup = MCNPACAB.get_user_input(outpfile)

# tally = MCNPACAB.tally_compose(tally0, Passive_sector)
# cmatrix = MCNPACAB.comp_matrix(tally0, Passive_sector)

for arg in sys.argv:
    if '-decay_times' in arg and not esc_file0:
        dec_times = arg.split('=')[-1]
        t_times = [float(time) for time in dec_times.split(',')]
        outputs_q = [0]
        for time in t_times:
            if float(time) >= 3600:
                outputs_q.append(1)
            else:
                outputs_q.append(0)
        print('Generating sutomatic scenario file')
        MCNPACAB.Escenary_generator(irr_time0, t_times, outputs_q, Sce_name ='Auto_Sce_file')
        esc_file0 = 'Auto_Sce_file'

ncel = [int(cell0) for cell0 in (tally0.cells)]
print('obtained cell numbers')
irr_cell = [cel.oget(outpfile,ncell_i) for ncell_i in ncel]
print('obtained cell properties')

# Because there can be quite a lot of cells with the same material, it is interesting to cache them
mat=[]
matnumbers=np.zeros(0)
for ncell0 in irr_cell:
    try:
        Mindex = list(matnumbers).index(ncell0.mat)
        mat0 = material.mat(mat[Mindex].number)
        mat0.N = list(mat[Mindex].N)
        mat0.M = list(mat[Mindex].M)
        mat.append(mat0)
    except:
        mat.append(material.oget(outpfile,ncell0.mat))
        matnumbers=np.append(matnumbers,mat[-1].number)
print('obtained materials')
vol = tally0.mass

def MCNP_ACAB_Mapstar(n):
    outputs = MCNPACAB.MCNP_ACAB_Map(tally0=tally0,mater=mat[n],irr_cell=irr_cell[n],
                               irr_time=irr_time0,n_id=n,save=save0,esc_file=esc_file0,
                               passive_sector=passive_sector0,source=source0,
                               id_lib=nuc_lib,id_ILIB=id_Egroup,corte=corte0)
    # Outputs:
    #     0 Timesets (arrays of times)
    #     1 Decay= Bq as ACAB
    #     2 Gamma= PHOTONS/CCM/SEC (as ACAB)
    #     3 Heat= W/cm3
    #     4 Dose= mSv/h (ACAB is Sv/h)
    #     5 mol = mol
    return outputs

with Pool(6) as pool:
    totaldata = pool.map(MCNP_ACAB_Mapstar, range(tally0.ncells))

if not t_times:
    t_times = list(totaldata[0][0].index)
else:
    o_times = [1.0]
    for time in t_times:
        if time in list(totaldata[0][0].index) and time not in o_times:
            o_times.append(time)
    t_times = o_times

MCNPACAB.summary_table_gen(totaldata,tally0,t_times) # generates summary table for shorter time
