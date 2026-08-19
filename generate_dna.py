# -*- coding: utf-8 -*-
from math import sin, cos, pi
from geometry import A,G,T,C

def rot_trans(L, dphi, dd):
    # rotate L by an angel phi and shift along z by a distance d
    L1 = []
    for i in L:
        i1 = []
        r = i[2]
        phi = i[3] + dphi
        phi = round(phi, 3)
        z = i[4] + dd
        z = round(z, 3)
        i1.append(i[0])
        i1.append(i[1])
        i1.append(r)
        i1.append(phi)
        i1.append(z)
        L1.append(i1)
    return L1
 

def cc2car(L):
    # convert cylindrical coordinates to Cartesian
    L1 = []
    for i in L:
        i1 = []
        r = i[2]
        phi = i[3]
        phi = phi / 180.0 * pi
        x1 = r * cos(phi)
        x1 = round(x1, 3)
        y1 = r * sin(phi)
        y1 = round(y1, 3)
        z1 = i[4]
        i1.append(i[0])
        i1.append(i[1])
        i1.append(x1)
        i1.append(y1)
        i1.append(z1)
        L1.append(i1)
    return L1


def add2file(file, L, atom_seq, res_seq):
    # file: the open file to be written to
    L1 = cc2car(L)
    res_seq += 1
    for i in L1:
        atom_seq += 1
        atom_name = i[0]
        if len(atom_name) <= 3: #这一点很重要，否则占位不对
            atom_name = ' ' + atom_name
        element = atom_name.strip()[0]
        res_name = 'D' + i[1]
        x = i[2]
        y = i[3]
        z = i[4]
        line = "%4s  %5i %-4s %3s  %4i    %8.3f%8.3f%8.3f                      %2s  " \
               % ("ATOM", atom_seq, atom_name, res_name, res_seq, x, y, z, element)
        file.write(line + '\n')
    file.write('TER\n')
    return atom_seq, res_seq


def ssDNA(file, seq, phi_inc, d_inc, atom_seq, res_seq):
    # for B-form DNA, phi_inc = 36.0, d_inc = 3.38
    nres = len(seq)
    for i in range(nres):
        L = seq[i]
        if L == 'A':
            L = A
        if L == 'G':
            L = G
        if L == 'C':
            L = C
        if L == 'T':
            L = T
        phi = phi_inc * (res_seq)
        d = d_inc * (res_seq)
        L = rot_trans(L, phi, d)
        atom_seq, res_seq = add2file(file, L, atom_seq, res_seq)
    return atom_seq, res_seq

'''
def ssDNA_complement(file, seq1, phi_inc, d_inc, atom_seq, res_seq):
    # for B-form DNA, phi_inc = 36.0, d_inc = 3.38
    seq2 = ''
    for L in seq1:
        if L == 'A':
            seq2 += 'T'
        if L == 'G':
            seq2 += 'C'
        if L == 'C':
            seq2 += 'G'
        if L == 'T':
            seq2 += 'A'
    nres = len(seq2)
    for i in range(nres):
        L = seq2[i]
        if L == 'A':
            L = A
        if L == 'G':
            L = G
        if L == 'C':
            L = C
        if L == 'T':
            L = T
        phi = phi_inc * (res_seq - nres) + 180.0
        d = d_inc * (res_seq - nres)
        L = rot_trans(L, phi, d)
        atom_seq, res_seq = add2file(file, L, atom_seq, res_seq)
    return atom_seq, res_seq
'''

def gen_DNA(filename, seq, phi_inc, d_inc, single):
    file = open(filename, 'w')
    atom_seq, res_seq = ssDNA(file, seq, phi_inc, d_inc, 0, 0)
    #if not single:
        #atom_seq, res_seq = ssDNA_complement(file, seq, phi_inc, d_inc, atom_seq, res_seq)
    file.close()
    return atom_seq, res_seq

if __name__ == '__main__':
    seq = 'C'
    natom, nres = gen_DNA('single_DNA.pdb', seq, 36.0,3.38, single = True)