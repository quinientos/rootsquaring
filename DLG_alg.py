import mpmath as mp
from math import log2

mp.mp.dps = 30
precision = 30

add = lambda x,y : mp.fadd(x, y, exact = True)
sub = lambda x,y : mp.fsub(x, y, exact = True)
mul = lambda x,y : mp.fmul(x, y, exact = True)
div = lambda x,y : mp.fdiv(x, y) #dps = precision)
pod = lambda x,y : mp.power(x, y)
mpf = lambda x	 : mp.mpf(x)		#dps = precision)
mpc = lambda x,y : mp.mpc( mpf(x), mpf(y))

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} so that
#circ_sq_root(p,q) = r,s and y = e^{2*pi*i*(r/s)}
#implies y^2 = x
def circ_sq_root(p,q):
	if p%q == 0:
		return 1,1
	else:
		return p, 2*q

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} so that
#circ_neg(p,q) = t,u and y = e^{2*pi*i*(t/u)}
#implies -y = x
def circ_neg(p,q):
	if p%q == 0:
		return 1, 2
	else:
		return 2*p+q, 2*q

#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} and an integer l
#so that circ_roots_rational_form(p,q,l) = r_1,s_1;...;r_{2^l},s_{2^l}
#and y_i = e^{2*pi*i*(r_i/s_i)} implies y_i^{2^l} = x
#with the r_i,s_i ordered according to the DLG recursion
def circ_roots_rational_form(p,q,l):
	r, s	= circ_sq_root(p,q)
	t, u	= circ_neg(r,s)
	if l == 1:
		return [(r,s),(t,u)]
	elif l != 0:
		left	= circ_roots_rational_form(r,s,l-1)
		right = circ_roots_rational_form(t,u,l-1)
		left.extend(right)
		return left
	else:
		return [(p,q)]

circ_roots_rational_form_memo_tree = {(1,1,1):((1,1))}

#A memoization of circ_roots_rational_form
def circ_roots_rational_form_memo(p,q,l):
	if (p,q,l) in circ_roots_rational_form_memo_tree:
		return list(circ_roots_rational_form_memo_tree[(p,q,l)])
	else:
		new_roots = circ_roots_rational_form(p,q,l)
		circ_roots_rational_form_memo_tree[(p,q,l)] = tuple(new_roots)
		return new_roots


#the input to this function is a representation
#of a complex number x = e^{2*pi*i*(p/q)} and an integer l
#so that circ_roots(p,q,l) = y_1,...,y_{2^l} where y_i^{2^l} = x
#with the y_i ordered according to the DLG recursion
def circ_roots(p,q,l):
	#Uncomment here to turn off memoization
	#roots = circ_roots_rational_form(p,q,l)
	roots = circ_roots_rational_form_memo(p,q,l)
	return [mp.expjpi(mul(2,div(r,s))) for r,s in roots]

#the input to this function is a representation of a complex number
#x = r*e^{2*pi*i*(t/u)}, and an integer l
#so that roots(r,t,u,l) = y_1,...,y_{2^l} where y_i^{2^l} = x
#with the y_i ordered according to the DLG recursion
def roots(r,t,u,l):
	circ_root				= circ_roots(t,u,l)
	r_root					= mp.root(r,2**l)
	return [mul(r_root, root) for root in circ_root]

#the input to this function is a black box polynomial p, its derivative p',
#a representation of a complex number x = r*e^{2*pi*i*(t/u)}, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG_rational_form(p,dp,r,t,u,l):
	root			 = roots(r,t,u,l)

	#alternate: divided diff approach
	pvals = [p(root[i]) for i in range(len(root))]
	delta = 2**-(mp.mp.prec//(2*log2(len(root)-1)))
	base_step = [div(sub(pvals[i], p(sub(root[i],delta))), mul(pvals[i], delta)) for i in range(len(root))]
	derivs		 = [base_step]
	for i in range(l):
		derivs.append([])
		for j in range(2**(l-i-1)):
			derivs[i+1].append(div(sub(derivs[i][2*j], derivs[i][2*j+1]),mul(root[2*j],2)))
		root = roots(r,t,u,l-1-i)
	return derivs[l][0]

#the input to this function is a black box polynomial p, its derivative p',
#a complex number x, and an integer l
#so that circ_DLG(p,dp,r,t,u,l) = (1/2z) p'_i(z)/p_i(z) - p'_i(z)/p_i(z) were
#z = x^(-1/2); i.e., circ_DLG(p,dp,r,t,u,l) = "the l^th difference of p'/p at x"
def DLG(p,dp,x,l,epsilon):
	angle = mp.arg(x)
	u			= pod(2,epsilon)
	t			= mp.fmod(angle,pod(2,-epsilon))
	r			= mpf(mp.fabs(x))
	return DLG_rational_form(p,dp,r,t,u,l)
