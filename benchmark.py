import os
import sys
import shlex
import subprocess
import re
import argparse
import time
from math import factorial, log2, ceil
import mpmath as mp
import sum_methods as sm
import DLG_alg as dlg

#SETTINGS
COEFF_PREC = 200
TEST_PREC = 200

#Polynomial class to store function details
class Polynomial:
    def __init__(self, coeff_data):
        self.c = 0
        if type(coeff_data) == list:
            self.fromList(coeff_data)
        else:
            self.fromFile(coeff_data)
    
    def fromList(self, p_coeffs):
    # coeffs are in dense format in the order of [lc -> tc]
    # this is in the reverse order of the mpsolve input format
    # but in line with the numpy poly format
        deg = len(p_coeffs)-1
        p = lambda x: mp.polyval(p_coeffs, mp.fsub(x,self.c))

        self.deg = deg
        self.p = p
        self.dp = lambda x: mp.diff(self.p, x)
        self.p_rev = lambda x: mp.fmul(self.p(mp.power(x,-1)), mp.power(x, self.deg))
        self.dp_rev = lambda x: mp.diff(self.p_rev, x)
        self.coeffs = {deg-i: p_coeffs[i] for i in range(deg+1)}
    
    def fromFile(self, pol_file):
    # Takes in MPSolve input files in the following formats:
    # dri drq drf dci dcq sci sri srf
        with open(pol_file) as f:
            L = f.read().split('\n')
    
        L = [l for l in L if l.strip() != '']
        while L[0][0] == '!': L.pop(0)
        pol_type = L.pop(0)
        L.pop(0) # 0 line. 
        deg = int(L.pop(0))
    
        if pol_type[2] == 'i': 
            num_lines = 1
            num_fcn = lambda x: mp.mpf(x)
        elif pol_type[2] == 'f':
            num_lines = 1
            num_fcn = lambda x: mp.mpf(x)
        elif pol_type[2] == 'q':
            num_lines = 2
            num_fcn = lambda x: mp.mpf(x)
    
        if pol_type[1] == 'r':
            num_comps = 1
        elif pol_type[1] == 'c':
            num_comps = 2
    
        if pol_type[0] == 's':
            num_terms = int(L.pop(0))
            p_degs = []
            for i in range(num_terms):
                p_degs.append(int(L.pop(i*(num_lines*num_comps+1)-i)))
    
        #L left over is p_coeffs
        p_coeffs = list(map(num_fcn, L))
        if pol_type[2] == 'q':
        #coeffs are rationals
            p_coeffs = [ mp.fdiv(x,y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
        
        if pol_type[1] == 'c':
        #coeffs are complex numbers
            p_coeffs = [ mp.mpc(real=x, imag=y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
        
        if pol_type[0] == 'd':
        #dense format
            p_coeffs.reverse() #the input files of mpsolve are in deg 0 -> deg d order
            p_degs = [deg-i for i in range(deg+1)]
            p = lambda x: mp.polyval(p_coeffs, mp.fsub(x, self.c))
        
        elif pol_type[0] == 's':
        #sparse format giving coeff/degree of the monomials
            p = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(mp.fsub(x, self.c),z)), p_coeffs, p_degs)))
    
        self.deg = deg
        self.p = p
        self.dp = lambda x: mp.diff(self.p, x)
        self.p_rev = lambda x: mp.fmul(self.p(mp.power(x,-1)), mp.power(x, self.deg))
        self.dp_rev = lambda x: mp.diff(self.p_rev, x)
        self.coeffs = {p_degs[i]: p_coeffs[i] for i in range(len(p_degs))}
    
    def get_trailing_coeffs(self, N, rev=False):
    # returns the coefficients of N+1 trailing coefficients p_0, ..,  p_N
        if rev == False:
            tcs =  [(self.coeffs[i] if i in self.coeffs.keys() else 0) for i in range(min(N+1, self.deg+1))] + [0 for i in range(self.deg+1, N+1)]
            return tcs
        #reverse poly
        tcs =  [(self.coeffs[self.deg-i] if (self.deg-i) in self.coeffs.keys() else 0) for i in range(min(N+1, self.deg+1))] + [0 for i in range(self.deg+1, N+1)]
        return tcs

    def set_shift(self, c):
        self.c = c


# runs mpsolve on the given input file and finds roots and max/min root radii
# returns the roots and the extremal root radii
def get_root_radii(pol_file):
    cmd_str = 'mpsolve '+ pol_file
    Cmd = shlex.split(cmd_str)
    start_time = time.time()
    output = subprocess.check_output(Cmd).strip().split(b'\n')
    mpsolve_time = time.time() - start_time
    num_pat = b'\((.+), (.+)\)'
    roots = []
    radii = []
    orig_radii = []

    for i in range(len(output)):
        r = re.match(num_pat, output[i]).groups()
        rc = mp.mpc(real=float(r[0]), imag=float(r[1]))
        radii.append(mp.fabs(rc))

    r_min = min(radii) #min radii
    r_max = max(radii) #max radii

    return roots, r_min, mpsolve_time

def get_err(x, y):
    if y == 0: return "inf"
    return mp.fabs(mp.fdiv(mp.fsub(x,y), y))

def get_numbers(poly, rootsum, r_min, d, k):
    s = mp.fabs(rootsum)
    if s == 0:
        bound_rd = "inf"
        bound_rd_err = "-"  #1.0
    else:
        bound_rd = mp.root(mp.fabs(mp.fdiv(d, s)), k)
        bound_rd = mp.fabs(mp.fadd(bound_rd, mp.fabs(poly.c)))
        bound_rd_err = get_err(bound_rd, r_min)

    return bound_rd, bound_rd_err



#accuracy test against MPSolve
def get_bound(poly, l, roots, r_min, mpsolve_time, pol_family=None):
    start_time = time.time()
  
    #default DLG eval precision
    e = int(mp.mp.prec*0.95)
    angle = mp.rand()
    x = mp.expjpi(angle*2)
    x = mp.fsub(0,mp.fmul(x,mp.power(2,-e)))

    dlg_out = dlg.DLG(poly.p, poly.dp, x, l, e) #sum 1/x_i, use for x_d. 
    dlg_out = mp.fabs(dlg_out)
    rd, rd_err = get_numbers(poly, dlg_out, r_min, poly.deg, 2**l)
    dlg_time = time.time() - start_time
    
    ##simple refinement, assuming there is a convenient way to obtain the LC.
    #lc = poly.coeffs[poly.deg]
    #new_rd = mp.fabs(mp.root(mp.fdiv(poly.p(0), lc), poly.deg))
    #new_rd_err = get_err(new_rd, r_min)

    return rd, rd_err

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="polynomial file in mpsolve input format")
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    pol_file = args.infile
    mp.mp.dps = COEFF_PREC
    poly = Polynomial(pol_file)
    roots, r_min, mpsolve_time = get_root_radii(pol_file)

    ##for precision, adjust this following line
    mp.mp.dps = TEST_PREC
    print("decimal precision used: %s" % mp.mp.dps)
    print("bit precision used: %s" % mp.mp.prec)

    l = ceil(log2(log2(poly.deg)))
  
    rd_bd, rd_err = get_bound(poly, l, roots, r_min)

if __name__ == '__main__':
    main()
