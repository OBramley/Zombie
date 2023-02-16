# from lib2to3.pygram import Symbols
# from sympy import symbols, summation
from sympy import *
from sympy.parsing.latex import parse_latex
# from sympy.parsing.latex import parse_latex

init_printing()
# symbols('zs:ndet')
# # i,j,k,l,m,n,ndet,d,w1,w2,zetaa,zetab = symbols('i,j,k,l,m,n,ndet,d,w1,w2,zetaa,zetab')
# # var('zetaa,zetab')
# a=
# x=Symbol('x')
# expr=
# expr = integrate(x**x,x)
# expr

latex_string = r"\Delta(k) = \frac{\rho_1-\rho_2}{\rho_1 + \rho_2} gk + \frac{\gamma k^3}{\rho_1 + \rho_2} - \frac{\rho_1 \rho_2}{(\rho_1 + \rho_2)^2} U^2 k^2"
equation = parse_latex(latex_string)
pprint(equation)


# # equation = parse_latex(latex_string)
# pprint(expr)

# i, n = symbols('i n')
# s, x = symbols('s x', cls=Function)
# s = summation(x(i), (i, 1, n))
# frac = x(i)/s
# diff(frac,x(i))