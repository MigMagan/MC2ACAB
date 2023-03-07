#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 16:54:35 2023

@author: mmagan
"""
import re

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
