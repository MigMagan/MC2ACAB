#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 16:54:35 2023

@author: mmagan
"""
import re
from math import radians, cos
from numpy import linalg, cross, array, zeros, transpose, matmul

def input_finder(infile):
    "Get the input from MCNP output infile"
    inputlines = []
    with open(infile,"r", encoding='utf-8') as outp:
        nlineoffset=0 # This is the # of chars the code reserves to print
        #the line number in the output (plus the '-'). 6 for MCNPX, 11 for MCNP6
        initline=outp.readline()
    #        print(initline)
        if initline.find('mcnpx')>0:
            nlineoffset=6
        if initline.find('MCNP6')>0:
            nlineoffset=11
        if initline.find('MCNP_6.20')>0:
            nlineoffset=11
        if nlineoffset==0:
            print("Can't find MCNP version in header. If it is indeed an MCNP output,\
                  ensure mct in the PRDMP card is not set to -1 ")
        for line in outp:
            if re.findall("[0-9]+- ",line[:nlineoffset+7]):
                inputlines.append(line[nlineoffset+7:])
    return inputlines

def line_parser(line):
    "Take a MCNP input line, and return its tokens, removing comments"
    tokens = []
    if line[0] not in ["c","C",]:
        for token in re.split(' +|=|\n', line):
            if token in ["&", "$"]:
                break
            tokens.append(token) # comment
    return tokens

def __interval_unfold(tokenlist, token_type=float):
    """"Unfold a list of intervals in the shape
    X ni Y. Does this for all intervals found in tokenlist 
    token_type indicate the type of tokens used. It is float by default, use
    ints to get the interval unfolded as integers. BEWARE, it does not check
    validity"""
    uf_tokenlist=[]
    for i, token in enumerate(tokenlist):
        if re.match("(\d+)i", token, re.IGNORECASE):
            n = int(token[:-1])+1  # Not very happy with this, would rather use the Match object
            starti = token_type(tokenlist[i-1])
            stopi = token_type(tokenlist[i+1])
            ilist =[str(starti + (stopi-starti)//n*(j+1)) for j in range(n-1)]
            uf_tokenlist.extend(ilist)
        else:
            uf_tokenlist.append(token)
    return uf_tokenlist


def __build_TR(tokens):
    """Internal to build a TR Matrix using the tokens from the TR Card, taking into
    account the multiple cases possible"""
    while '' in tokens:
        tokens.remove('')
    if '*' in tokens[0]:
        for i, token in enumerate(tokens[4:13]):
            tokens[i+4] = f'{cos(radians(float(token))):.4f}' # to convert from angles to cosene of the angles
    match len(tokens):
        case 4:  # Pure traslation
            TR = array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [float(tk) for tk in tokens[0:3]]])
        case 7:  # One vector
            print("WARNING: using a single-vector defined rotation at MCNP results in arbitrary\
                   vectors and is a borderline bad practice. Returning because we don't know how to\
                   reproduce MCNPs fickle way to define the TR card")
            TR = array[[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]]
        case 9:  # Reconstitute with Euler vector.
            print("WARNING: Eulerian angles does not apparently work well in MCNP and is not\
                  implemented here. Please make a full or 2-vector definition of the rotation matrix")
        case 10:  # Two vectors, cross product.
            TRx = array([float(token) for token in tokens [4:7]])
            TRy = array([float(token) for token in tokens [7:10]])
            TRx/=linalg.norm(TRx)
            TRy/=linalg.norm(TRy)
            TRz = cross(TRx,TRy)
            TRxv = array([TRx[0], TRy[0], TRz[0]])
            xx = linalg.norm(TRxv)**2
            xy = TRx[0]*TRx[1]+TRy[0]*TRy[1]+TRz[0]*TRz[1]
            TRxv/=linalg.norm(TRxv)
            TRyv = array([TRx[1]*xx-TRx[0]*xy, TRy[1]*xx-TRy[0]*xy, TRz[1]*xx-TRz[0]*xy])
            TRyv/=linalg.norm(TRyv)
            TRzv = cross(TRxv, TRyv)
            TR = zeros((4,3))
            TR[0] = [float(t) for t in tokens[1:4]]
            TR[1:4] = transpose([TRxv, TRyv, TRzv])
        case 13|14:  # Fully described matrix
            TRx = array([float(token) for token in tokens [4:7]])
            TRy = array([float(token) for token in tokens [7:10]])
            TRz = array([float(token) for token in tokens [10:13]])
            for T in [TRx, TRy, TRz]:
                T/=linalg.norm(T)
            TRxv = array([TRx[0], TRy[0], TRz[0]])
            xx = linalg.norm(TRxv)**2
            xy = TRx[0]*TRx[1]+TRy[0]*TRy[1]+TRz[0]*TRz[1]
            TRxv/=linalg.norm(TRxv)
            TRyv = array([TRx[1]*xx-TRx[0]*xy, TRy[1]*xx-TRy[0]*xy, TRz[1]*xx-TRz[0]*xy])
            TRyv/=linalg.norm(TRyv)
            TRzv = cross(TRxv, TRyv)
            TR = zeros((4,3))
            TR[0] = array([float(t) for t in tokens[1:4]])
            TR[1:4] = transpose([TRxv, TRyv, TRzv])
    if len(tokens) == 14:
        if tokens[13] == '-1':
            TR[0] = matmul(linalg.inv(TR[1:4]),TR[0])
    return TR

def get_TR(trID, outpinfile):
    "get the transformation matrix from a TR Card from an MCNP output outpinfile"
    in_outp = input_finder(outpinfile)
    tokens = []
    for i, line in enumerate(in_outp):
        if re.match(rf'\*?TR{trID}', line, re.IGNORECASE):
        # if line.startswith('tr'+str(trID)) or line.startswith('*tr'+str(trID)):
            l_fmesh = i
            j = 0
            tokens.extend(line_parser(in_outp[l_fmesh]))
            while in_outp[l_fmesh+1+j].startswith(' '*5):
                tokens.extend(line_parser(in_outp[l_fmesh+1+j]))
                j += 1
            return __build_TR(tokens)
    else:
        print(f"TR card {trID} not found")
        return array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

def get_histpcells(outpinfile):
    '''Look into an outp file for histp entry and return the cell array'''
    in_outp = input_finder(outpinfile)
    tokens = []
    for i, line in enumerate(in_outp):
        if re.match('histp', line, re.IGNORECASE):
            l_histp = i
            j = 0
            tokens.extend(line_parser(in_outp[l_histp])[1:])
            while in_outp[l_histp+1+j].startswith((' '*5, 'c', 'C')):
                tokens.extend(line_parser(in_outp[l_histp+1+j]))
                j += 1
    while '' in tokens:
        tokens.remove('')
    if not tokens:
        print("Card HISTP not found")
        return []
    if int(tokens[0]) < 0:
        tokens.pop(0)
    uf_tokens = __interval_unfold(tokens, token_type=int)
    histp_cells = array(uf_tokens,dtype=int)
    return histp_cells