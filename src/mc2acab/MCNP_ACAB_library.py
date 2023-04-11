#! /usr/bin/env python

# Version Feb 2023
''' Module to generates ACAB inputs and perform ACAB automatic calculations
    By Miguel Magan Romero and Dr. Octavio Gonzalez-del Moral updating an original
    code from Dr. Fernando Sordo'''

import subprocess
import time
import os
import sys
import shutil
import re
import datetime
import numpy as np
import pandas as pd
import apypa
import tally as tal
from . import pyhtape3x


def __is_number(s):
    '''Check if a number is a number'''
    try:
        float(s)
        return True
    except ValueError:
        return False

intervals = (
    ('years', 31536000), # 60 * 60 * 24 * 365
    ('months', 2592000), # 60 * 60 * 24 * 30
    ('weeks', 604800),  # 60 * 60 * 24 * 7
    ('days', 86400),    # 60 * 60 * 24
    ('hours', 3600),    # 60 * 60
    ('minutes', 60),
    ('seconds', 1),
    )

def __display_time(seconds, granularity=2):
    ''' Chanege any time expresion in seconds to human readable units \n
    granularity: level of detail of the transformation\n'''
    result = []
    for name, count in intervals:
        value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip('s')
            result.append(f"{int(value)}_{name}")
    if not result:
        result.append('Shutdown')
    return '+'.join(result[:granularity])

def __backup_previous(item):
    ''' Check if a file exits and changes its name to save it from being overwrite'''
    if os.path.exists(item):
        date = datetime.datetime.now()
        name_id=str(date.strftime("%Y%m%d%H%M%S"))
        os.replace(item,f'{item}_bk_{name_id}')
        print(f"\nBacking up existing {item} as {item}_bk_{name_id}")

def __get_user_source():
    print('Calculation of source intensity:')
    print('1: Source term units in particles per second (part/s):')
    print('2: Source term units in mA (miliAmperes): ')
    while True:
        units = input('Choose: ')
        if units in ['1','2']:
            break
        print('Please, just 1 or 2')
    user_input = input('Input source term value: ')
    source = float(user_input)
    if units == '2':
        source *= 6.24E15
        print(f'Source term = {source:.2E} n/s')
    return source

def __get_user_time():
    user_input = input('Enter irradiation time in hours: ')
    try:
        irr_time = float(user_input) * 3600  # convert to seconds
        return irr_time
    except ValueError:
       print('Invalid input for irradiation time')
       sys.exit(1)

def __get_user_LIB():
    while True:
        Nuc_lib = input('Enter the Nuclear library (Default = EAF): ')
        if Nuc_lib in ['EAF','']:
            break
        print('Please, just EAF')
    Nuc_lib = 'EAF' if Nuc_lib == '' else 'EAF' # There are one single option....

    while True:
        id_Egroup = input('Enter energy groups ID (vitJ+ or other), default = vitJ+: ')
        if id_Egroup in ['vitJ+', 'VitJ+', 'other', 'Other','']:
            break
        print('Please, just vitJ+ or other')
    id_Egroup = 'vitJ+' if id_Egroup == '' else id_Egroup # default option
    print('Read from user', Nuc_lib, id_Egroup)
    return Nuc_lib, id_Egroup

def get_user_input(outp_name='outp'):
    """
    Prompts the user to input values for source term and irradiation time, reads data from a file type "outp",
    and returns a modified tally object, irradiation time in seconds, and source term value in particles per second.
    """
    # Check if the outp file exists
    if not os.path.exists(outp_name):
        raise FileNotFoundError('outp file not found')

    source = __get_user_source()
    irr_time = __get_user_time()
    Nuc_lib, id_Egroup = __get_user_LIB()

    with open(outp_name,'r', encoding='utf-8') as datafile:
        nps_list = []
        lines = datafile.readlines()
        for line in lines:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0] == 'dump':
                nps_list.append(float(words[8]))
        nps = nps_list[-1]
    print(f'Last nps recorded: {nps}')
    print(tal.olist(outp_name))

    while True:
        user_input = input('Enter the tally number for the cell of interest: ')
        try:
            ntal = int(user_input)
        except ValueError:
            print('Invalid tally number')
            continue
        tally = tal.oget(outp_name, ntal)
        if tally:
            break
        print(f'Tally {ntal:d} is not here...\n Please, check it and re-enter it:')

    return tally, irr_time, source, Nuc_lib, id_Egroup


def collapse(tally, source, **kwargs):
    """
    Escribe un archivo de entrada COLL.inp de ACAB COLLAPS a partir del flujo de tally y lo ejecuta.
    Args:
        tally: objeto de tally
        source: fuente de partículas
        id_lib: identificador de la librería de secciones transversales
        id_ilib: identificador de la librería de grupos de energía
        cell: número de celda
        surface: número de superficie
        cos: número de ángulo coseno
        t_time: tiempo de irradiación
    """

    id_lib = kwargs.get('id_lib','EAF')
    id_ilib = kwargs.get('id_ilib','vitJ+')
    cell = kwargs.get('cell',0)
    surface = kwargs.get('surface',0)
    cos = kwargs.get('cos',0)
    t_time = kwargs.get('t_time',0)

    print("*********** RUNNING ESPECTRO-4-ACAB **********")
    print(f"Using tally {tally.n} with total flux {tally.value[cell, surface, cos, t_time, -1]}")
    print(f"Total N1: {tally.value[cell, surface, cos, t_time, -1]:.3e} parts/cm2 per source particle ")
    print(f"Total N2: {tally.value[cell, surface, cos, t_time, -1] * source:.3e} parts/cm2 s")
    print(f"Volume: {tally.mass[cell][surface]:.3f}")
    print("*********************************************\n")

    num_lines = 16 if id_lib == 'EAF' else 32

    if id_ilib == 'vitJ+':
        ilib = 12
        iesf = 12
        ngroup = 211
    else:
        print('WARNING!!! Unknown energy group structure DO NOT USE!!! ACAB works funny')
        while True:
            user_input = input('Are you sure? You are proposing an unknown energy groups distributuion (Y/N) ')
            if user_input in ['Y','Yes','y','N','No','n']:
                break
            print('Please, Y/N')
        if user_input in ['Y','Yes','y']:
            print(f'Using non-standard energy groups, getting them from tally {tally.n}')
            ilib = 5 # 5 do not exits in ACAB, so it do not work
            iesf = 5
            ngroup = tally.eints - 1
        else:
            return

    with open('COLL.inp', 'w', encoding='utf-8') as outfile:
        outfile.write(f'{ilib} {iesf}\n')  # card 1 ILIB IESF
        outfile.write(f'{num_lines:d}\n')  # card 2  IHEAD
        outfile.write('0 0 0 0\n')  # card 3 ISFIS IGEN ISOCA IBEST if ISFIS = 0 Other values are ignored
        # card 4 not present because ISFIS = 0
        outfile.write(f'-{ngroup} 0\n')  # card 5 NGROUP decreasing energy, FF = 0 units n/cm2-s
        # card 6 not present if IESF != 5
        if iesf == 5:
            ebins_st = ['{:.5e}'.format(float(ebin)) for ebin in tally.ebins[:-1]]
            ebins_st.reverse()
            outfile.write('\n'.join([' '.join(ebins_st[i:i+6]) for i in range(0,tally.eints, 6)]))
            outfile.write('\n')
        # card 7 Flux levels by decreasing energy groups Normalize by n/cm2 s
        str_flux = ['{:.5e}'.format(flux_i*source) for flux_i in tally.value[cell,surface,cos,t_time,:ngroup]]
        str_flux.reverse()
        outfile.write('\n'.join([' '.join(str_flux[i:i+6]) for i in range(0,tally.eints, 6)]))
        outfile.write('\n0\n')  # card 8 IUNC3G
        outfile.write('0\n')  # card 9 ISTOP

    print("*********** RUNNING COLLAPS **********")
    subprocess.run(['collaps_2008'], check=True)


def scenary_generator(irr_time, cooling_times, outputs, **kwargs):
    ''' Generates an automatic ACAB Scenario file '''
    # TODO Think how to improve generator for more complex scenarios... i.e. multiple irradiation cycles
    Sce_name = kwargs.get('Sce_name', None)
    feeds = kwargs.get('feeds', None)
    if not isinstance(cooling_times,list) or not any(__is_number(time) for time in cooling_times):
        raise ValueError('Cooling times must we a list of times (integers in secons)')
    if not isinstance(outputs,list) or any(output not in [0,1] for output in outputs):
        raise ValueError('Outputs times must we a list of 0 (NO output) or 1 (output) '
              'plus an initial 0 correspoding to the irradiation')
    if len(outputs) != len(cooling_times) + 1:
        raise ValueError('Outputs times must we a list of 0 (NO output) or 1 (output) '
              'plus an initial 0 correspoding to the irradiation')
    isfed = 1 if feeds else 0
    inputfile = ''
    inputfile +="<  Bloque de acumulacion y quemado\n"
    inputfile += f"10  10  1  10  1  {isfed:d}  0  0\n"
    irr_times = ['{:.3e}'.format(irr_time * pow(2,-10+1+i)) for i in range(10)]
    inputfile += '\n'.join([', '.join(irr_times)])
    prev_nsteps = 10
    prev_time = 0 # To keep tracks of previous cooling time
    for i, ctime in enumerate(cooling_times):
        inputfile += f"\n< Bloque de decaimiento {i}\n"
        if ctime < 1000:
            isend = 0 if ctime == cooling_times[-1] else 1 
            inputfile += f"0  1  {isend}  {prev_nsteps}  1  0  0  0\n"
            inputfile += f"{ctime:.3E}"
            prev_nsteps = 1
        else:
            isend = 0 if ctime == cooling_times[-1] else 1
            inputfile += f"0  10  {isend}  {prev_nsteps}  1  0  0  0\n"
            dec_times = ['{:.3e}'.format(((ctime - prev_time)*pow(2,i)/pow(2,9))+prev_time) for i in range(10)]
            inputfile += '\n'.join([', '.join(dec_times)])
            prev_nsteps = 10
        prev_time = ctime
    # BLOCK 9
    inputfile += "\n1.0E-25 1.0E+00\n"
    # BLOCK 10
    inputfile += "0 0 0 \n"
    #BLOCK 11:
    inputfile += "<  Block 11 \n"
    # Card 1 IWP IMTX IWDR       IDOSE IPHCUT IDHEAT     IOFFSD ICEDE INEMISS IDAMAGE
    inputfile += "1 0 1  1 1 0  0 0 0 0 \n"
    inputfile += "0 0 1 0 \n" # Card 2
    #  Card 6 //   NOPUL NTSEQ NOTTS NVFL
    j_total = len(cooling_times) + 1
    inputfile += f"0 0 {j_total:d} 1 \n"
    # Card #7
    inputfile += ' '.join(['1'] * j_total)
    # No Card #8 and this is neutron problem, no BLOCK #12
    inputfile += "\n0\n" # Just a 0... 
    # BLOCK 13
    # Card 1 NCYO IFSO
    inputfile += "0 1\n"
    # No Card 2. Notice that ACAB will STILL label next card as Card #2
    str_outputs = [str(output) for output in outputs]
    inputfile += '\n'.join([' '.join(str_outputs[i:i+8]) for i in range(0,j_total, 8)])
    inputfile += '\n'
    if Sce_name != None:
        with open (Sce_name, "w", encoding='utf-8') as writefile:
            writefile.write(inputfile)
    return inputfile

def create_inp(flux,irr_time,mat,vol,**kwargs):
    ''' Create an ACAB imput file inp.5 kwargs can be:
        sce_file: irradiation scenario file. Renders irr_time irrelevant
        feeds: External isotopical feed, typically for proton activation'''
    sce_file = kwargs.get('sce_file', None)
    feeds = kwargs.get('feeds', None)
    if sce_file is not None:   # Scenario file provided
        print("Using Scenario file", repr(sce_file))
        if not os.path.exists(str(sce_file)):
            raise FileNotFoundError(f'Warning!!! No Scenario File!!! {str(sce_file)} is not here')
        else:
            with open(sce_file, 'r', encoding='utf-8') as infile:
                sce_str = infile.read()
    else:
        print("Building up an automatic scenario file")
        # Decaimiento:                 1 h    1day    1week   1month    1year   5years
        cooling_times = [0.1, 1, 10, 100, 3600, 3600*24, 3600*24*7, 3600*24*30, 3600*24*365, 3600*24*365*5]
        outputs = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1] # One more element than Cooling_time (the irradiation!)
        sce_str = scenary_generator(irr_time, cooling_times, outputs, Sce_name=None, feeds=feeds)

    with open ("inp.5", "w", encoding='utf-8') as inputfile:
        libreria = 2232
        # Block 1
        inputfile.write("Entrada generada por MCNP_ACAB\n") # Card 1
        inputfile.write("0\n") # Card 2
        ITMAX = 250000 if flux == 0 else 900000
        isfed = 1 if feeds else 0
        inputfile.write(f'{libreria:d} {ITMAX:d}  0   1  1 0 2    2 '
                        f'{2*isfed}    24  1  0  4  1  1  0  1  0  0  0  0 \n') #  Card 3
        # Block 2
        inputfile.write("{0:E}     1.0000\n".format(vol))
        inputfile.write("1\n")
        Niso=len(mat.zaid)
        # inputfile.write("{0:d}\n".format(Niso)) # Block 2 Card 4: Numero de isotopos
        inputfile.write("{}\n".format(Niso)) # Block 2 Card 4:
        # Card 5 ISOZO
        if feeds is not None:
            inputfile.write("{0:d}\n".format(len(feeds[1])))
        #Card 6
        Elist = ['2.0e+01', '1.4e+01', '1.2e+01', '1.0e+01', '8.0e+00', '6.5e+00',
                 '5.0e+00', '4.0e+00', '3.0e+00', '2.5e+00', '2.0e+00', '1.7e+00',
                 '1.4e+00', '1.2e+00', '1.0e+00', '8.0e-01', '6.0e-01', '4.0e-01',
                 '3.0e-01', '2.0e-01', '1.0e-01', '5.0e-02', '2.0e-02', '1.0e-02', 
                 '0.0e+00']
        inputfile.write('\n'.join([', '.join(Elist[i:i+8]) for i in range(0,len(Elist), 8)]))
        inputfile.write('\n')
        # Card 7
        inputfile.write("0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00   \n")
        # Card 8: Outputs
        inputfile.write("0 0 1   0 0 0   0 0 1   0 0 1   0 0 0   0 0 0   \n")
        #Block 3: Neutron flux
        inputfile.write("{0:E}\n".format(flux))
        # Block 4:
        inputfile.write("0 IREST\n")
        #Block 5:
        inputfile.write("< Isotopia\n")
        str_mats = [str(isotope) for isotope in mat.zaid]
        inputfile.write('\n'.join([' '.join(str_mats[i:i+5]) for i in range(0,len(str_mats), 5)]))
        inputfile.write('\n')
        str_mats = [f'{isotope:.6e}' for isotope in mat.frac]
        inputfile.write('\n'.join([' '.join(str_mats[i:i+5]) for i in range(0,len(str_mats), 5)]))
        inputfile.write('\n')
        #Blocks 6: Feeds.
        if feeds is not None:
            str_mats = [str(isotope) for isotope in feeds[1]]
            inputfile.write('\n'.join([', '.join(str_mats[i:i+5]) for i in range(0,len(str_mats), 5)]))
            inputfile.write('\n')
            str_mats = [str(isotope) for isotope in feeds[2]]
            inputfile.write('\n'.join([', '.join(str_mats[i:i+5]) for i in range(0,len(str_mats), 5)]))
            inputfile.write('\n')
        #Blocks 7-8: Burn out and cooldown scenario
        inputfile.write(sce_str)

def MCNP_ACAB_Map(**kwargs):
    #tally0,mater,irr_cell,irr_time,irr_type,n_id,save,esc_file,passive_sector,source,id_lib,id_ILIB,corte):
    '''Assistant to carry ouy MCNP_ACAB calculations'''
    feeds = None  # Default value
    start_time = time.time()
    tally0 = kwargs.get('tally0')
    mater = kwargs.get('mater')
    irr_cell = kwargs.get('irr_cell')
    irr_time = kwargs.get('irr_time')
    irr_type = kwargs.get('irr_type','n')
    n_id = kwargs.get('n_id')
    source = kwargs.get('source')
    save = kwargs.get('save', False)
    sce_file0 = kwargs.get('esc_file', None)
    id_lib = kwargs.get('id_lib', 'EAF') # the only one that works in ACAB
    id_ILIB = kwargs.get('id_ILIB', 'vitJ+') # the only one that works in ACAB
    corte = kwargs.get('corte', 1E-2) 
    print('particles: ',irr_type)
    Wdir = str(irr_cell.ncell)
    print(f"doing cell {irr_cell.ncell}")
    flux = tally0.value[n_id, 0, 0, 0, -1]*source
    if not mater.zaid:
        print("null material")
        return None
    if flux == 0:
        print("null tally")
        return None
    __backup_previous(Wdir)
    os.mkdir(Wdir)
    os.chdir(Wdir)
    if sce_file0 is not None:
        os.symlink(os.pardir+os.sep+sce_file0, sce_file0)
# Manipulamos el mat para que pueda representar estados excitados
    mater.n2ro(irr_cell.density)
    matfixed_zaid, matfixed_frac = pyhtape3x.unfold_NA(mater.zaid, mater.frac)
    mater.zaid[:] = [10*i for i in matfixed_zaid] # Fix material with nat abundance AND add excited state info
    mater.frac[:] = list(matfixed_frac)
    vol = tally0.mass[n_id, 0]
     # Parte de enlazar *.dat
    DatFiles=["DHEAT.dat","FYBL.dat","af_asscfy.dat","PHOTON.dat","MACOEF.dat","EBEATA.dat","DECAY.dat","WD.dat"]
    Dat_origin_Files=[]
    for datfile in DatFiles:
        dfile = os.environ["ACAB_LB_PATH"]+datfile
        Dat_origin_Files.append(dfile)
    Dat_origin_Files[1] = os.environ["ACAB_LB_PATH"]+"eaf_n_fis_20070"
    Dat_origin_Files[2] = os.environ["ACAB_LB_PATH"]+"eaf_n_asscfy_20070" 
    for index,datfile in enumerate(DatFiles):
        if not os.path.isfile(datfile):    
            os.symlink(Dat_origin_Files[index],datfile)
#    print('\033[31m flux {0}, tally_ncel {1}, n {2}\033[0m'.format(tally.value[n][-1],tally.cells[n],n))
    if 'n' in irr_type:
        xsfile = f"{os.environ['ACAB_LB_PATH']}eaf_n_gxs_211_flt_20070"
    else:
        xsfile = f"{os.environ['ACAB_LB_PATH']}eaf_p_gxs_211_flt_20070"
    if not os.path.isfile('XSBL.dat'):
        os.symlink(xsfile, 'XSBL.dat')
    
    if  re.match(r"[^pn]", irr_type):
        print("particle type not valid")
        os.chdir(os.pardir)
        return None

    if 'p' in irr_type:  # Deal with the isotopical feeds
        __backup_previous("RES_H")
            print("RES_H file not already present")
        pyhtape3x.createRSH(irr_cell.ncell)
        os.symlink("../histp", "./histp") 
        os.system("htape3x int=RSH outt=RES_H")
        feeds=pyhtape3x.get_atom_feed(irr_cell.ncell,"RES_H")
        feeds[2][:]=[source/6.023E23*i for i in feeds[2]]

    collapse(tally0, source, id_lib=id_lib, id_ilib=id_ILIB, cell=n_id)
    create_inp(flux,irr_time,mater,vol,sce_file=sce_file0,feeds=feeds)

    print("*********** RUNNING ACAB 2008 **********")
    subprocess.run('acab_2008',check=True)
    heat = apypa.heat_isotopes_full_pd("fort.6",threshold = corte) 
    gamma = apypa.gammas_full_pd("fort.6")
    dose = apypa.gamma_dose_isotopes_full_pd("fort.6", threshold = corte)
    # timesets = apypa.get_time_sets('fort.6')
    decay = apypa.rad_act_isotopes_full_pd('fort.6', threshold = corte)
    mol = apypa.iso_mol('fort.6', threshold = corte)
    if save in  ['All','all']:
        os.chdir(os.pardir)
    elif save == True:
        print('\nRemoving REACTIONS.dat and XSECTION.dat\n')
        os.remove('REACTIONS.dat')
        os.remove('XSECTION.dat')
        os.chdir(os.pardir)
    else:
        os.chdir(os.pardir)
        print('Removing working directory '+str(Wdir))
        shutil.rmtree(Wdir)
# Calculamos el tiempo de ejecución
    elapsed_time=time.time()-start_time
# Escribimos la linea en el log
    with open('logfile.txt','a', encoding='utf-8') as logfile:
        line = [f'cell/voxel={Wdir}/{irr_cell.ncell}',f'vol={vol:.2e}ccm', f'ro={irr_cell.density*-1:.2f}g/ccm', 
                f'SourceTerm={source/6.24E15:.3e}mA',f'NeutronFlux={flux:.2e}part/s',
                f'IrrTime={__display_time(irr_time,3)}', f'Run=-{str(irr_type)}',
                f'Time={elapsed_time//60:.0f}m {elapsed_time%60:.2f}s']
        [logfile.write(f'{item} ') for item in line]
        logfile.write('\n')
        logfile.close()
    # Outputs:
    #     1 decay= Bq as ACAB
    #     2 gamma= PHOTONS/CCM/SEC (as ACAB)
    #     3 heat= W/cm3
    #     4 dose= mSv/h (ACAB is Sv/h) 
    #     5 mol = mol 
    return decay, gamma, heat ,dose, mol

def summary_table_gen(totaldata_ACAB,tally,**kwargs):
    """ Script de generacion de tablas resumen de MCNP_ACAB """
    t_times = kwargs.get('t_times','All')
    save = kwargs.get('save',True)
    if t_times in ['all','All']:
        for data in totaldata_ACAB:
            if data != None:
                t_times = list(data[0].index)
                break
    panda_item = np.dtype([('cell',int),('vol',float),('decay',object),
                          ('gamma',object),('heat',object),('dose',object),('mol',object)])
    apypas = np.zeros(len(totaldata_ACAB), dtype=panda_item)
    __backup_previous('summary_apypas.npy')
    for i, pd_list in enumerate(totaldata_ACAB):
        if pd_list != None:
            apypas[i] = tally.cells[i], tally.mass[i], pd_list[0], pd_list[1], pd_list[2], pd_list[3], pd_list[4]
    np.save('summary_apypas',apypas)
    for voxel in apypas:
        totals = pd.DataFrame()
        totals.index.name = f'Cell:{voxel[0]} Vol:{float(voxel[1]):.2e}'
        for panda in list(voxel)[2:]:
            if isinstance(panda,pd.DataFrame):
                if 'Total' in panda.columns:
                    totals[f'Total_{panda.columns.name}'] = panda['Total']
                else:
                    totals[f'Total_{panda.columns.name}'] = panda.T['Total']
        totalsT = totals.T
        for time_i in totalsT.columns:
            if time_i not in t_times:
                totalsT = totalsT.drop(columns = time_i)
        totals = totalsT.T
        totals = totals.round(3)
        totals = totals.applymap('{:.3e}'.format)
        __backup_previous('summary_ACAB_{tally.cells[i]}.csv')
        if save == True:
            totals.to_csv(f'summary_ACAB_{tally.cells[i]}.csv',sep='\t',encoding='utf-8')
    # apypas = np.load('summary_apypas.npy', allow_pickle=True)
    return apypas
