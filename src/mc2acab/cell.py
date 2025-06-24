#! /usr/bin/env python
# coding: utf-8
import re
from mc2acab import MCNP_outparser

def __is_number(s):
    "Internal to check if a token is a number"
    try:
        float(s)
        return True
    except ValueError:
        return False

class Cell:
    """
    This is an MCNP cell object.
    """

    def __init__ (self, ncell):
        self.ncell=ncell
        self.mat=0
        self.density=0
        self.volume=0
        self.NIMP=0
        self.PIMP=0
        self.EIMP=0
        self.HIMP=0
# ====================================================== #

def __parse_cell(inp):
    "Internal to get a cell from the MCNP input lines inp that define it"
    impn, impp, imph, impe=(1, 1, 1, 1)
    celldef = ''.join(list(inp))
    tokens = MCNP_outparser.line_parser(celldef)
    n = int(tokens[0])
    mat=tokens[1]
    if int(mat)!=0:
        ro = tokens[2]
    else:
        ro = 0
    for i,t in enumerate(tokens):
        if re.match(r'imp:[p,h,d,n,t,e,\/,\|]*n[p,h,d,n,t,e,\/,\|]*', t,
            flags=re.IGNORECASE):
            impn=tokens[i+1]
        if re.match(r'imp:[p,h,d,n,t,e,\/,\|]*h[p,h,d,n,t,e,\/,\|]*', t,
            flags=re.IGNORECASE):
            imph=tokens[i+1]
        if re.match(r'imp:[p,h,d,n,t,e,\/,\|]*e[p,h,d,n,t,e,\/,\|]*', t,
            flags=re.IGNORECASE):
            impe=tokens[i+1]
        if re.match(r'imp:[p,h,d,n,t,e,\/,\|]*p[p,h,d,n,t,e,\/,\|]*', t,
            flags=re.IGNORECASE):
            impp=tokens[i+1]
    Cel = Cell(n)
    Cel.mat = int(mat)
    Cel.density = float(ro)
    Cel.NIMP = int(float(impn))
    Cel.EIMP = int(float(impe))
    Cel.HIMP = int(float(imph))
    Cel.PIMP = int(float(impp))
    # print (mat,ro,impn,impe,imph)
    return Cel

def oget(infile,n):
    """ Get the cell n from MCNP output infile"""
    inp = MCNP_outparser.input_finder(infile)
    nstr=str(n)
    found_sline = False
    for i, lines in enumerate(inp):
        tokens = MCNP_outparser.line_parser(lines)
        if not tokens:  # Empty line
            continue
        if lines[0]==' ':  #Not a cell line
            continue
        ncell = tokens[0]
        if ncell == nstr:
            sline = i  # starting line
            found_sline = True
        elif found_sline:
            eline = i  # finish line
            break
    else:
        print("Cell not found.")
        return None
    return __parse_cell(inp[sline:eline])

def _table60(infile):
    '''confirm if infile has table 60'''
    with open(infile, 'r', encoding='utf-8') as inp_file:
        if 'print table 60' in inp_file.read():
            return True
        else:
            return False

def _find_table60(infile):
    table60 = []
    dentrotabla = False
    with open(infile, 'r', encoding='utf-8') as inp_file:
        for linea in inp_file:
            if 'print table 60' in linea:
                dentrotabla = True
                continue
            if 'total' in linea and dentrotabla:
                break
            if dentrotabla:
                table60.append(linea)
    return table60

def _read_table60(cel, table):
    ''' Read table 60 for specific cell to get all useful info'''
    for line in table:
        if len(line.split()) > 0 and str(cel.ncell) == line.split()[1]:
            cel.mat = float(line.split()[2])
            if cel.density < 0:
                cel.density = float(line.split()[4])*-1  # keep negative as readed
            cel.volume = float(line.split()[5])
            break
    return cel

def ogetall(infile):
    """ Get an array of all cells of MCNP output infile """
    inp = MCNP_outparser.input_finder(infile)
    cellist = []
    nlines = []  # List of lines where a line is defined
    for i, lines in enumerate(inp[1:]):
        if re.match("read file", lines, flags=re.IGNORECASE) is not None:
            inp.pop(i+1)  # remove read file lines
    for i, lines in enumerate(inp[1:]):
        tokens = MCNP_outparser.line_parser(lines)
        if not tokens:
            continue
        if lines[0]==' ':
            continue
        if all([tok == '' for tok in tokens]):
            nlines.append(i+1)
            break
        nlines.append(i+1)  # since we skipped the 1st title line
    for i, nline in enumerate(nlines[:-1]):  # To exclude the line that marks end of cell definition
        cel = __parse_cell(inp[nline:nlines[i+1]])
        cellist.append(cel)
    if _table60(infile):
        table60 = _find_table60(infile)
        print('Table 60 readed')
        for cel in cellist:
            cel = _read_table60(cel, table60)
    return cellist
