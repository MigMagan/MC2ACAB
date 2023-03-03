#! /usr/bin/env python
import re

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def createRSH(cell, RSHfile="RSH"):
    "Creat a RSH file to get the residues for cell"
    with open (RSHfile,"w", encoding="UTF-8") as RSH:
        RSH.write("Entrada de residuos para HTAPE3X\n")
        line = "Calculo en la celda "+str(cell)+"\n"
        RSH.write(line)
#            IOPT, NERG, NTIM, NTYPE, KOPT, NPARM, NFPRM, FNORM, KPLOT,IXOUT, IRS, IMERGE, ITCONV, IRSP, ITMULT/
        RSH.write("8,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0\n")
        line = str(cell)+'\n'
        RSH.write(line)

def create_composedRSH(acell, pcell, RSHfile="RSH"):
    """
    Create a RSH cell for a composed analysis for rotating elements.
    Passivecell is the matching cell of the passive sector
    n1 is the number of elements, and n2 is the volume ratio (passive/active)
    """
    with open (RSHfile,"w", encoding="UTF-8") as RSH:
        RSH.write("Entrada de residuos para HTAPE3X\n")
        line="Calculo en las celdas "+str(acell)+" y "+str(pcell) +"\n"
        RSH.write(line)
#            IOPT, NERG, NTIM, NTYPE, KOPT, NPARM, NFPRM, FNORM, KPLOT,IXOUT, IRS, IMERGE, ITCONV, IRSP, ITMULT/
        RSH.write("8,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0\n")
        line=str(acell)+","+str(pcell)+'\n'
        RSH.write(line)

def get_atom_feed(cell, Resfile):
    "Wrapper of getatomfeeds() for a single cell"
    res = get_atom_feeds([cell], Resfile)
    return res[0]

def get_atom_feeds(cells, Resfile):
    """
    Returns the arrays of residual nuclei of cells in the htape3x output file Resfile.
    cells is either a single value or a list/array
    """
    try:
        cells = list (cells)
    except TypeError:  #  cells is a single value, not iterable
        cells = [cells]
    # Define the overarching structure
    inputstr2 = []
    found = [False]*len(cells)
    searchstr="distribution of residual nuclei in cell"
    sline = [0]*len(cells)
    eline = [0]*len(cells)
    with open(Resfile,"r", encoding="UTF-8") as htape3xfile:
        for i, lines in enumerate(htape3xfile):
            if re.findall(searchstr, lines):
                ncell = int(lines.split()[-1])
                print(f"found cell {ncell}")
                idx = cells.index(ncell) if ncell in cells else None
                if idx is None:
                    continue
                sline[idx] = i
                found[idx] = True
            if re.findall("all z", lines) and found[idx]:
                eline[idx] = i+1
            if re.findall("metastable state", lines):
                print ("metastable isotopes found")
                while lines != "end":
                    lines = next(htape3xfile, "end")
                    inputstr2.append(lines)
        htape3xfile.seek(0, 0)
        inputstrs = [htape3xfile.readlines()[sl:el] for sl, el in zip (sline, eline)]
        residuals = [__readfeed(inputstr) for inputstr in inputstrs]
    residuals = __addmetastable(residuals, inputstr2)
    if residuals==[]:
        print ("Did not find any residual nuclei. Are you sure this is a correct htape3x output?")
        return None
    return residuals

def __readfeed(inputstr):
    "Internal function to read the feeds from a list of lines inputstr"
    Z, A, N, m, l =[[] for i in range(5)]
    for lines in inputstr:
        palabras=lines.split()
        if not palabras: continue # empty line
        if re.findall("residual nuclei in cell", lines):
            ncell = int(palabras[-1])
        if palabras[0]=="z" and palabras[1]=='a' :
            print('RES_H file empty')
            break # RES_H file empty
        if palabras[0]=="z" and palabras[3]=='n' :
            Z.append(int(palabras[2]))
            N.append(int(palabras[5]))
            A.append(Z[-1]+N[-1])
            fraction=float(palabras[6].replace("D","E"))
            m.append(fraction)
            l.append(0)
        if palabras[0]=="n":
            Z.append(Z[-1])
            N.append(int(palabras[2]))
            A.append(Z[-1]+N[-1])
            fraction=float(palabras[3].replace("D","E"))
            m.append(fraction)
            l.append(0)
        if palabras[0]== "all":  # summary line
            if palabras[1]=="z":
                break # we are done
    isotopes, feeds =([], [])
    for index,isotope in enumerate(A):
        value=Z[index]*10000+isotope*10+l[index] # Get the ACAB value of the isotope
        isotopes.append(value)
        feeds.append(m[index])
    residual=[ncell, isotopes, feeds]
    return residual

def __addmetastable(residuals, inputstr):
    """Internal to add metastable isotopes to residuals
    using a list of output strings inputstr"""
    prev_ID=0
    for lines in inputstr:
        palabras=lines.split()
        if not palabras:continue
        if is_number(palabras[0]):
            Z=(int(palabras[0]))
            A=(int(palabras[1]))
            l=(int(palabras[2]))
            if l==0:
                print(f"WARNING: substituting apparently bugged metastable for Z={Z} A={A}")
                l=1
            # N=A-Z
            fraction=float(palabras[5].replace("D","E"))
            ID=Z*10000+A*10+l
            if ID==prev_ID:
                l=l+1
            fundamental_ID=Z*10000+A*10
            prev_ID=ID
            for residual in residuals:
                if fundamental_ID in residual[1]:  # Encontramos el núcleo. SI NO ESTÁ, NO SE TOCA NADA.
                    residual[1].append(ID)
                    i=residual[1].index(fundamental_ID)
                    total=residual[2][i]
                    residual[2].append(total*fraction)  # ADD the excited nuclide.
        if palabras[-1]=="completed":
            break
    # With all the excited states added, correct the inventory of fundamental states.
    for residual in residuals:
        for i,ID in enumerate(residual[1]):
            if ID % 10 !=0:
                fundamental_ID=ID // 10 *10
                j=residual[1].index(fundamental_ID)
                residual[2][j] = residual[2][j]-residual[2][i]
    return residuals




def nat_abun(at_number):
    """ This function returns two lists, one with the mass numbers, and one with the abundances of the isotopes"""

    Z=[3,5,6,12,14,16,17,19,20,22,23,24,26,28,29,30,32,34,37,38,40,42,44,46,50,51,52,56,72,74,76,78,81,82]
    A=[]
    Abun=[]
    A.append([6,7])  # Li
    Abun.append([0.075,0.925])

    A.append([10,11])  #B
    Abun.append([0.199,0.801])

    A.append([12,13])  #C
    Abun.append([0.989,0.011])

    A.append([24,25,26])  #Mg
    Abun.append([0.7899,0.1,0.1101])

    A.append([28,29,30])  #Si
    Abun.append([0.9223,0.0467,0.0310])

    A.append([32,33,34,36])  #S
    Abun.append([0.9502,0.0075,0.0421,2E-4])

    A.append([35,37])  #Cl
    Abun.append([0.7577,0.2423])

    A.append([39,40,41])  #K
    Abun.append([0.9326,1.2E-4,0.0673])

    A.append([40,42,43,44,46,48])  #Ca
    Abun.append([0.96941,0.00647,0.00135,0.02086,4E-5,0.00187])

    A.append([46,47,48,49,50])  #Ti
    Abun.append([0.0825,0.0744,0.7372,0.0541,0.0518])

    A.append([50, 51])  #V
    Abun.append([0.0025,0.9975])

    A.append([50,52,53,54])  #Cr
    Abun.append([0.04345,0.83789,0.09501,0.02365])

    A.append([54,56,57,58])  #Fe
    Abun.append([0.05845,0.9172,0.02119,0.00282])

    A.append([58,60,61,62,64])  #Ni
    Abun.append([0.68077,0.26223,0.0114,0.03634,0.00926])

    A.append([63,65])  #Cu
    Abun.append([0.6917,0.3083])

    A.append([64,66,67,68,70])  #Zn
    Abun.append([0.4863,0.279,0.041,0.1875,0.0062])

    A.append([70,72,73,74,76])  #Ge
    Abun.append([0.2123,0.2766,0.0773,0.3594,0.0744])

    A.append([74,76,77,78,80,82])  #Se
    Abun.append([0.0087,0.0936,0.0763,0.2378,0.4961,0.0873])

    A.append([85,87])  #Rubidio
    Abun.append([0.72168,0.27835])

    A.append([84,86,87,88])  #Sr
    Abun.append([0.0056,0.0986,0.07,0.8258])

    A.append([90,91,92,94,96])  #Zr
    Abun.append([0.5145,0.1122,0.1715,0.1738,0.028])

    A.append([92,94,95,96,97,98,100])  #Mo
    Abun.append([0.1484,0.0925,0.1592,0.1668,0.0955,0.2413,0.0963])

    A.append([96,98,99,100,101,102,104])  #Rutenio
    Abun.append([0.0552,0.0188,0.127,0.126,0.170,0.316,0.187])

    A.append([102,104,105,106,108,110])  #Pa
    Abun.append([0.0102,0.1114,0.2233,0.2733,0.2646,0.1172])

    A.append([112,114,115,116,117,118,119,120,122,124])  #Sn
    Abun.append([0.0097,0.0066,0.0034,0.1454,0.0768,0.2422,0.0859,0.3258,0.0463,0.0579])

    A.append([121,123])  #Sb
    Abun.append([0.5746,0.4264])

    A.append([120,122,123,124,125,126,128,130])  #Teluro
    Abun.append([9E-4,0.0255,0.0089,0.0474,0.0705,0.1884,0.3174,0.3408])

    A.append([130,132,134,135,136,137,138])  #Ba
    Abun.append([0.00106,0.00101,0.02417,0.06592,0.07854,0.11232,0.71698])

    A.append([174,176,177,178,179,180])  #Hf
    Abun.append([0.0016,0.0526,0.186,0.2728,0.1362,0.3508])

    A.append([180,182,183,184,186])  #W
    Abun.append([0.0012,0.265,0.1431,0.3064,0.2846])

    A.append([184,186,187,188,189,190,192])  #Os
    Abun.append([2E-4,0.0159,0.0196,0.1324,0.1615,0.2626,0.4078])

    A.append([191,193])  #Ir
    Abun.append([0.373,0.627])

    A.append([190,192,194,195,196,198])  #Pt
    Abun.append([0.00014,0.00782,0.32967,0.33832,0.25242,0.07163])

    A.append([203,205])  #Ta
    Abun.append([0.29524,0.70476])

    A.append([204,206,207,208])  #Os
    Abun.append([0.014,0.241,0.221,0.524])

    if at_number in Z:
        i=Z.index(at_number)
        isotopes=A[i]
        Abundance=Abun[i]
    else:
        print("Isotope Z={0} not found".format(at_number))
        return -1
    return isotopes,Abundance
# =======================================================
def atomic_mass(Z):
    """This function returns the average atomic mass of isotope Z using its natural composition"""
    mass=0
    A=nat_abun(Z)
    if A==-1:
        return -1 # Z not found

    for index,isotope in enumerate(A[0]):
        mass=mass+isotope*A[1][index]
    return mass

def unfold_NA(isotopes,feeds):

    """ This function takes an array of isotopes with its concentrations (or feeds) and gives back an array with the natural abundances 'unfolded' """
    import operator

    unfolded_isotopes=[isotope for isotope in isotopes if operator.mod(isotope,1000)!=0]
    unfolded_feeds=[feed for index,feed in enumerate(feeds) if operator.mod(isotopes[index],1000)!=0 ]
    for index,isotope in enumerate(isotopes):
        Z1,resto=divmod(isotope,1000)
        if resto==0:
            A=nat_abun(Z1)
            if A==-1:
                print("Could not find natural abundance for  z={0}".format(Z1))
                unfolded_isotopes.append(isotope)
                unfolded_feeds.append(feeds[index])
            else:
                for i,j in enumerate(A[0]):
                    n=Z1*1000+j
                    unfolded_isotopes.append(n)
                    n=A[1][i]*feeds[index]
                    unfolded_feeds.append(n)

    return unfolded_isotopes,unfolded_feeds
