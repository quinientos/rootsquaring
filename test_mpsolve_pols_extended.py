import os
import sys
import shlex
import subprocess
import re
import argparse
from DLG_alg_mpmath import DLG_rational_form, DLG, add, sub, mul, div, precision, mpc, mpf, pod
import newton_ids_alg as ni
import mpmath as mp
import math
import time

#mp.precision = 30

class Polynomial:
    def __init__(self, coeff_data):
        if type(coeff_data) == list:
            self.fromList(coeff_data)
        else:
            self.fromFile(coeff_data)
    
    def fromList(self, p_coeffs):
    # coeffs are in dense format in the order of [lc -> tc]
    # this is in the reverse order of the mpsolve input format
    # but in line with the numpy poly format
        deg = len(p_coeffs)-1
        p = lambda x: mp.polyval(p_coeffs, x)
        p_degs = [deg-i for i in range(deg+1)]
        p = lambda x: mp.polyval(p_coeffs, x)
        dp_coeffs = [p_coeffs[i]*(deg-i) for i in range(len(p_coeffs)-1)]
        dp = lambda x: mp.polyval(dp_coeffs, x)

        #not strictly needed but for posterity
        dp_degs = [d-1 for d in p_degs] 
        if -1 in dp_degs: 
            loc = dp_degs.index(-1) 
            dp_degs.pop(loc)
    
        #get rev polys
        p_rev_degs = [deg - j for j in p_degs]
        p_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_rev_degs)))
        dp_rev_degs = [d-1 for d in p_rev_degs]
        dp_rev_coeffs = [p_coeffs[i]*p_rev_degs[i] for i in range(len(p_rev_degs))]

        if -1 in dp_rev_degs:
            loc = dp_rev_degs.index(-1)
            dp_rev_degs.pop(loc)
            dp_rev_coeffs.pop(loc)
        dp_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_rev_coeffs, dp_rev_degs)))

        self.deg = deg
        self.p = p
        self.dp = dp
        self.p_rev = p_rev
        self.dp_rev = dp_rev
        #print("degs: %s" % p_degs)
        #print("coeffs: %s" % p_coeffs)
        self.coeffs = {p_degs[i]: p_coeffs[i] for i in range(len(p_degs))}

    def fromFile(self, pol_file):
        # Should handle the following formats:
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
            p_coeffs = [ div(x,y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
        
        if pol_type[1] == 'c':
            p_coeffs = [ mp.mpc(real=x, imag=y) for (x,y) in zip(p_coeffs[0::2], p_coeffs[1::2])]
        
        if pol_type[0] == 'd':
            p_coeffs.reverse() #the input files of mpsolve are in deg 0 -> deg d order
            p_degs = [deg-i for i in range(deg+1)]
            p = lambda x: mp.polyval(p_coeffs, x)
            dp_coeffs = [p_coeffs[i]*(deg-i) for i in range(len(p_coeffs)-1)]
            dp = lambda x: mp.polyval(dp_coeffs, x)
    
            #not strictly needed but for posterity
            dp_degs = [d-1 for d in p_degs] 
            if -1 in dp_degs: 
                loc = dp_degs.index(-1) 
                dp_degs.pop(loc)
        
        elif pol_type[0] == 's':
            p = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_degs)))
            dp_degs = [d-1 for d in p_degs] 
            dp_coeffs = [p_coeffs[i]*p_degs[i] for i in range(num_terms)]
    
            if -1 in dp_degs:
                loc = dp_degs.index(-1)
                dp_degs.pop(loc)
                dp_coeffs.pop(loc)
            dp = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_coeffs, dp_degs)))
    
        #get rev polys
        p_rev_degs = [deg - j for j in p_degs]
        p_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), p_coeffs, p_rev_degs)))
        dp_rev_degs = [d-1 for d in p_rev_degs]
        dp_rev_coeffs = [p_coeffs[i]*p_rev_degs[i] for i in range(len(p_rev_degs))]
    
        if -1 in dp_rev_degs:
            loc = dp_rev_degs.index(-1)
            dp_rev_degs.pop(loc)
            dp_rev_coeffs.pop(loc)
        dp_rev = lambda x: mp.fsum(list(map( lambda y,z: mp.fmul(y, mp.power(x,z)), dp_rev_coeffs, dp_rev_degs)))
    
        self.deg = deg
        self.p = p
        self.dp = dp
        self.p_rev = p_rev
        self.dp_rev = dp_rev
        #print("degs: %s" % p_degs)
        #print("coeffs: %s" % p_coeffs)
        self.coeffs = {p_degs[i]: p_coeffs[i] for i in range(len(p_degs))}
        #self.coeffs = p_coeffs
        #self.degs = p_degs
    
    def get_trailing_coeffs(self, N, rev=False):
    # returns the coefficients of N+1 trailing coefficients p_0, ..,  p_N
        if rev == False:
            tcs =  [(self.coeffs[i] if i in self.coeffs.keys() else 0) for i in range(min(N+1, self.deg+1))] + [0 for i in range(self.deg+1, N+1)]
            return tcs
        #reverse poly
        #print("deg = %d, N = %d" % (self.deg, N))
        tcs =  [(self.coeffs[self.deg-i] if (self.deg-i) in self.coeffs.keys() else 0) for i in range(min(N+1, self.deg+1))] + [0 for i in range(self.deg+1, N+1)]
        return tcs
            

# runs mpsolve on the given input file and finds roots and max/min root radii
# returns the roots and the extremal root radii
def get_root_radii(pol_file):
    cmd_str = 'mpsolve '+ pol_file
    Cmd = shlex.split(cmd_str)
    output = subprocess.check_output(Cmd).strip().split(b'\n')
    num_pat = b'\((.+), (.+)\)'
    roots = []
    radii = []

    for i in range(len(output)):
        r = re.match(num_pat, output[i]).groups()
        rc = mp.mpc(real=float(r[0]), imag=float(r[1]))
        roots.append(rc)
        rad = mp.fabs(rc)
        radii.append(rad)

    r_min = min(radii) #min radii
    r_max = max(radii) #max radii
    print("min rt radius = %s" % mp.nstr(r_min, 3))
    print("max rt radius = %s" % mp.nstr(r_max, 3))
    return roots, r_min, r_max


# runs tests using the DLG algorithm via recurrence and tcs algorithms
#def run_tests(deg, p, dp, p_rev, dp_rev, roots, r_min, r_max):
def run_tests(poly, roots, r_min, r_max):
    deg = poly.deg
    p = poly.p
    dp = poly.dp
    p_rev = poly.p_rev
    dp_rev = poly.dp_rev

    l = int(math.log2(deg))+ 4
    #l = int(math.log2(math.log2(deg))) +1
    #x is the point which defines a line to 0 on which we are taking a limit
    angle = mp.rand()
    x = mp.expjpi(angle*2)

    #star_min = time.time()

    print("l=%s" % l)
    #approx = div(deg,DLG(p,dp,sub(0,mul(x,mp.power(2,-e))),l,e))
    rat = ni.get_root_radii_via_Ris(poly, l)
    print("rat = ", rat)
    if rat == 0:
        #l = l+1
        #print("trying again with l=%d due to 0" % l)
        #rd = ni.get_root_radii_via_Ris(poly, l)
        #print("rd = ", rd)
        rd = "inf"
        #print("(p'/p)^(\ell) = 0. Skipping analysis")
        min_rel_err = "inf"
    else:
        approx = div(deg, rat)
        rd = mp.root(mp.fabs(approx), mp.power(2,l))
        #approx = div(deg, ni.get_root_radii_via_Ris(poly, l))
        print("approx=",mp.fabs(approx))
        real = mp.power(r_min,mp.power(2,l))
        #print("radius=", mp.fabs(real))
        #print("error=", mp.fabs(div(sub(real,approx),(mp.fabs(real)))))
        #print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
        #print("radius_root=", mp.fabs(r_min))
        #print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_min))))
        #print("rel_error_root=", int(mul(100,div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_min))),mp.fabs(r_min)))),"%")
        min_rel_err = mp.nstr(div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_min))),mp.fabs(r_min)),3)

    #l_max = int(math.log2(math.log2(deg)))
    #approx = div(DLG(p_rev,dp_rev,sub(0,mul(x,mp.power(2,-e))),l,e),deg)
    rat = ni.get_root_radii_via_Ris(poly, l, rev=True)
    print("rat_rev = ", rat)
    if rat == 0:
        print("(p_rev'/p_rev)^(\ell) = 0.")
        #max_rel_err = "-"
        #l_max += 1
        #print("trying again with l=%d due to 0" % l_max)
        #rd = ni.get_root_radii_via_Ris(poly, l_max, rev=True)
        #print("rd = ", rd)
    #else:
    approx = mp.fabs(div(rat, deg))
    r1 = mp.root(mp.fabs(approx), mp.power(2,l))
    #approx = div(ni.get_root_radii_via_Ris(poly, l, rev=True), deg)
    #print("approx=",mp.fabs(approx))
    real = mp.power(r_max,mp.power(2,l))
    #print("radius=", mp.fabs(real))
    #print("error=", mp.fabs((mp.fabs(real))-(mp.fabs(approx)))/(mp.fabs(real)))
    #print("approx_root=",mp.root(mp.fabs(approx), mp.power(2,l)))
    #print("radius_root=", mp.fabs(r_max))
    #print("error_root=", mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_max))))
    #print("rel_error_root=", int(mul(100,div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_max))),mp.fabs(r_max)))),"%")
    max_rel_err = mp.nstr(div(mp.fabs(sub(mp.fabs(mp.root(mp.fabs(approx), mp.power(2,l))), mp.fabs(r_max))),mp.fabs(r_max)),3)
    
    #duration = time.time()-star_min
    #print("time=",time.time()-star_min)
    
    #print("%s & %s & %s & %s & %s & %s & %s & %s & $[%s, %s]$\\\\" % (deg, l, mp.mp.dps, mp.nstr(rd,3), min_rel_err, mp.nstr(r1,3), max_rel_err, round(duration,2), mp.nstr(r_min,3), mp.nstr(r_max,3) ))
    #print("%s & %s & %s & %s & %s & %s & $[%s, %s]$\\\\" % (deg, l, mp.nstr(rd,3), min_rel_err, mp.nstr(r1,3), max_rel_err, mp.nstr(r_min,3), mp.nstr(r_max,3) ))
    print("%s & %s & %s & %s & $[%s, %s]$\\\\" % (deg, l, min_rel_err, max_rel_err, mp.nstr(r_min,3), mp.nstr(r_max,3) ))

# reads in given arguments
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=str, help="polynomial file in mpsolve input format")
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    infile = args.infile
    #deg, p, dp, p_rev, dp_rev = get_pols(infile)
    poly = Polynomial(infile)

    #l = int(math.log2(poly.deg))+3
    l = int(math.log2(math.log2(poly.deg)))
    extra_precision = int(l/3)
    mp.mp.dps = precision + 2**extra_precision + 100

    roots, r_min, r_max = get_root_radii(infile)
    #run_tests(deg, p, dp, p_rev, dp_rev, roots, r_min, r_max)
    run_tests(poly, roots, r_min, r_max)


if __name__ == '__main__':
    main()
