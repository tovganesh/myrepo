# mobihf v0.1, Quantum Chemistry on mobile!, (c) V.Ganesh, GPL
# Based on quantumj (http://quantumj.dev.java.net)
#      and PyQuante (http://pyquante.sf.net)
# Python 3 syntax adaptation 
from math import sin,cos,acos,pi,sqrt,floor,pow,exp,log,erf
import math
import string
from string import *

radii = {"H":0.23, "C":0.77, "N":0.75, "O":0.73}
atmnum = {"H":1, "C":6, "N":7, "O":8}
    
basis_sto = {
    'H': [('S',[
          (3.425251, 0.154329),
          (0.623914, 0.535328),
          (0.168855, 0.444635)])],  
    'C': [('S',[
          (71.616837, 0.154329),
          (13.045096, 0.535328),
          (3.530512, 0.444635)]),
        ('S',[
          (2.941249, -0.099967),
          (0.683483, 0.399513),
          (0.222290, 0.700115)]),
        ('P',[
          (2.941249, 0.155916),
          (0.683483, 0.607684),
          (0.222290, 0.391957)])],
    'N': [('S',[
          (99.106169, 0.154329),
          (18.052312, 0.535328),
          (4.885660, 0.444635)]),
        ('S',[
          (3.780456, -0.099967),
          (0.878497, 0.399513),
          (0.285714, 0.700115)]),
        ('P',[
          (3.780456, 0.155916),
          (0.878497, 0.607684),
          (0.285714, 0.391957)])],
    'O': [('S',[
          (130.709321, 0.154329),
          (23.808866, 0.535328),
          (6.443608, 0.444635)]),
        ('S',[
          (5.033151, -0.099967),
          (1.169596, 0.399513),
          (0.380389, 0.700115)]),
        ('P',[
          (5.033151, 0.155916),
          (1.169596, 0.607684),
          (0.380389, 0.391957)])]
    }

curbas = basis_sto

powers = {
    'S': [(0, 0, 0)],
    'P': [(1, 0, 0), (0, 1, 0), (0, 0, 1)]    
    }

class molecule:
    def __init__(self):
        self.noOfAtoms = 0
        self.title = ""
        self.atoms = []


def boys_function(m, T):
    F = [0.0] * (m + 1)
    if T < 1e-9:
        for i in range(m + 1):
            F[i] = 1.0 / (2 * i + 1) - T / (2 * i + 3)
        return F
        
    sqrt_T = sqrt(T)
    exp_T = exp(-T)
    
    F[0] = sqrt(pi) * (0.0 if sqrt_T == 0 else (1.0 if sqrt_T == 0 else erf(sqrt_T))) / (2 * sqrt_T)
    
    for i in range(m):
        F[i+1] = ((2 * i + 1) * F[i] - exp_T) / (2 * T)
        
    return F

class OSEvaluator:
    def __init__(self):
        pass
        
    def run(self, A, normA, pwrA, alphaA, B, normB, pwrB, alphaB, 
                 C, normC, pwrC, alphaC, D, normD, pwrD, alphaD):
        zeta = alphaA + alphaB
        eta = alphaC + alphaD
        rho = (zeta * eta) / (zeta + eta)
        
        P = [(alphaA * A[i] + alphaB * B[i]) / zeta for i in range(3)]
        Q = [(alphaC * C[i] + alphaD * D[i]) / eta for i in range(3)]
        
        RPQ2 = sum((P[i] - Q[i])**2 for i in range(3))
        T = rho * RPQ2
        
        L_total = sum(pwrA) + sum(pwrB) + sum(pwrC) + sum(pwrD)
        F = boys_function(L_total, T)
        
        RAB2 = sum((A[i] - B[i])**2 for i in range(3))
        RCD2 = sum((C[i] - D[i])**2 for i in range(3))
        
        K_AB = exp(-alphaA * alphaB * RAB2 / zeta)
        K_CD = exp(-alphaC * alphaD * RCD2 / eta)
        
        prefactor = (2 * pow(pi, 2.5)) / (zeta * eta * sqrt(zeta + eta)) * K_AB * K_CD
        
        self.memo = {}
        self.P = P
        self.Q = Q
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.zeta = zeta
        self.eta = eta
        self.rho = rho
        self.F = F
        self.prefactor = prefactor
        
        return self.vrr(pwrA, pwrB, pwrC, pwrD, 0) * normA * normB * normC * normD

    def vrr(self, a, b, c, d, m):
        key = (a, b, c, d, m)
        if key in self.memo:
            return self.memo[key]
            
        if sum(a) == 0 and sum(b) == 0 and sum(c) == 0 and sum(d) == 0:
            return self.prefactor * self.F[m]
            
        if sum(a) > 0:
            i = 0
            if a[1] > 0: i = 1
            if a[2] > 0: i = 2
            
            a_minus = list(a)
            a_minus[i] -= 1
            a_minus = tuple(a_minus)
            
            a_minus_2 = None
            if a[i] >= 2:
                tmp = list(a)
                tmp[i] -= 2
                a_minus_2 = tuple(tmp)
                
            PA = self.P[i] - self.A[i]
            val = PA * self.vrr(a_minus, b, c, d, m)
            
            coeff2 = (self.eta / (self.zeta + self.eta)) * (self.Q[i] - self.P[i])
            val += coeff2 * self.vrr(a_minus, b, c, d, m + 1)
            
            if a_minus[i] > 0: # a[i] >= 2
                c1 = (a[i] - 1) / (2 * self.zeta)
                val += c1 * (self.vrr(a_minus_2, b, c, d, m) - (self.rho / self.zeta) * self.vrr(a_minus_2, b, c, d, m + 1))
                
            if b[i] > 0:
                b_minus = list(b)
                b_minus[i] -= 1
                b_minus = tuple(b_minus)
                c2 = b[i] / (2 * self.zeta)
                val += c2 * (self.vrr(a_minus, b_minus, c, d, m) - (self.rho / self.zeta) * self.vrr(a_minus, b_minus, c, d, m + 1))
                
            if c[i] > 0:
                c_minus = list(c)
                c_minus[i] -= 1
                c_minus = tuple(c_minus)
                c3 = c[i] / (2 * (self.zeta + self.eta))
                val += c3 * self.vrr(a_minus, b, c_minus, d, m + 1)
                
            if d[i] > 0:
                d_minus = list(d)
                d_minus[i] -= 1
                d_minus = tuple(d_minus)
                c4 = d[i] / (2 * (self.zeta + self.eta))
                val += c4 * self.vrr(a_minus, b, c, d_minus, m + 1)
                
            self.memo[key] = val
            return val
            
        elif sum(b) > 0:
            i = 0
            if b[1] > 0: i = 1
            if b[2] > 0: i = 2
            
            b_minus = list(b)
            b_minus[i] -= 1
            b_minus = tuple(b_minus)
            
            a_plus = list(a)
            a_plus[i] += 1
            a_plus = tuple(a_plus)
            
            AB = self.A[i] - self.B[i]
            val = self.vrr(a_plus, b_minus, c, d, m) + AB * self.vrr(a, b_minus, c, d, m)
            self.memo[key] = val
            return val
            
        elif sum(c) > 0:
            i = 0
            if c[1] > 0: i = 1
            if c[2] > 0: i = 2
            
            c_minus = list(c)
            c_minus[i] -= 1
            c_minus = tuple(c_minus)
            
            c_minus_2 = None
            if c[i] >= 2:
                tmp = list(c)
                tmp[i] -= 2
                c_minus_2 = tuple(tmp)
                
            QC = self.Q[i] - self.C[i]
            WQ = (self.zeta / (self.zeta + self.eta)) * (self.P[i] - self.Q[i])
            
            val = QC * self.vrr(a, b, c_minus, d, m)
            val += WQ * self.vrr(a, b, c_minus, d, m + 1)
            
            if c_minus_2 is not None:
                c1 = (c[i] - 1) / (2 * self.eta)
                val += c1 * (self.vrr(a, b, c_minus_2, d, m) - (self.rho / self.eta) * self.vrr(a, b, c_minus_2, d, m + 1))
                
            if d[i] > 0:
                d_minus = list(d)
                d_minus[i] -= 1
                d_minus = tuple(d_minus)
                c2 = d[i] / (2 * self.eta)
                val += c2 * (self.vrr(a, b, c_minus, d_minus, m) - (self.rho / self.eta) * self.vrr(a, b, c_minus, d_minus, m + 1))
                
            self.memo[key] = val
            return val
            
        elif sum(d) > 0:
            i = 0
            if d[1] > 0: i = 1
            if d[2] > 0: i = 2
            
            d_minus = list(d)
            d_minus[i] -= 1
            d_minus = tuple(d_minus)
            
            c_plus = list(c)
            c_plus[i] += 1
            c_plus = tuple(c_plus)
            
            CD = self.C[i] - self.D[i]
            val = self.vrr(a, b, c_plus, d_minus, m) + CD * self.vrr(a, b, c, d_minus, m)
            self.memo[key] = val
            return val
            
        return 0.0


class OneElectronOS:
    def __init__(self):
        pass
        
    def overlap(self, A, normA, pwrA, alphaA, B, normB, pwrB, alphaB):
        zeta = alphaA + alphaB
        P = [(alphaA * A[i] + alphaB * B[i]) / zeta for i in range(3)]
        RAB2 = sum((A[i] - B[i])**2 for i in range(3))
        
        prefactor = (pi / zeta)**1.5 * exp(-alphaA * alphaB * RAB2 / zeta)
        
        self.memo_s = {}
        self.P = P
        self.A = A
        self.B = B
        self.zeta = zeta
        self.prefactor_s = prefactor
        
        return self.vrr_s(pwrA, pwrB) * normA * normB

    def vrr_s(self, a, b):
        key = (a, b)
        if key in self.memo_s:
            return self.memo_s[key]
            
        if sum(a) == 0 and sum(b) == 0:
            return 1.0 * self.prefactor_s
            
        if sum(a) > 0:
            i = 0
            if a[1] > 0: i = 1
            if a[2] > 0: i = 2
            
            a_minus = list(a)
            a_minus[i] -= 1
            a_minus = tuple(a_minus)
            
            PA = self.P[i] - self.A[i]
            val = PA * self.vrr_s(a_minus, b)
            
            if a_minus[i] > 0:
                a_minus_2 = list(a)
                a_minus_2[i] -= 2
                a_minus_2 = tuple(a_minus_2)
                val += (a[i] - 1) / (2 * self.zeta) * self.vrr_s(a_minus_2, b)
                
            if b[i] > 0:
                b_minus = list(b)
                b_minus[i] -= 1
                b_minus = tuple(b_minus)
                val += b[i] / (2 * self.zeta) * self.vrr_s(a_minus, b_minus)
                
            self.memo_s[key] = val
            return val
            
        elif sum(b) > 0:
            i = 0
            if b[1] > 0: i = 1
            if b[2] > 0: i = 2
            
            b_minus = list(b)
            b_minus[i] -= 1
            b_minus = tuple(b_minus)
            
            PB = self.P[i] - self.B[i]
            val = PB * self.vrr_s(a, b_minus)
            
            if b_minus[i] > 0:
                b_minus_2 = list(b)
                b_minus_2[i] -= 2
                b_minus_2 = tuple(b_minus_2)
                val += (b[i] - 1) / (2 * self.zeta) * self.vrr_s(a, b_minus_2)
                
            self.memo_s[key] = val
            return val
            
        return 0.0

    def kinetic(self, A, normA, pwrA, alphaA, B, normB, pwrB, alphaB):
        zeta = alphaA + alphaB
        P = [(alphaA * A[i] + alphaB * B[i]) / zeta for i in range(3)]
        RAB2 = sum((A[i] - B[i])**2 for i in range(3))
        self.prefactor_s = (pi / zeta)**1.5 * exp(-alphaA * alphaB * RAB2 / zeta)
        self.memo_s = {}
        self.P = P
        self.A = A
        self.B = B
        self.zeta = zeta
        
        s_ab = self.vrr_s(pwrA, pwrB)
        
        term_plus = 0.0
        for i in range(3):
            b_plus_2 = list(pwrB)
            b_plus_2[i] += 2
            term_plus += self.vrr_s(pwrA, tuple(b_plus_2))
            
        term_minus = 0.0
        for i in range(3):
            if pwrB[i] >= 2:
                b_minus_2 = list(pwrB)
                b_minus_2[i] -= 2
                term_minus += pwrB[i] * (pwrB[i] - 1) * self.vrr_s(pwrA, tuple(b_minus_2))
                
        lb = sum(pwrB)
        T = -0.5 * (4 * alphaB**2 * term_plus - 2 * alphaB * (2 * lb + 3) * s_ab + term_minus)
        
        return T * normA * normB

    def nuclear_attraction(self, A, normA, pwrA, alphaA, B, normB, pwrB, alphaB, C):
        zeta = alphaA + alphaB
        P = [(alphaA * A[i] + alphaB * B[i]) / zeta for i in range(3)]
        RAB2 = sum((A[i] - B[i])**2 for i in range(3))
        RPC2 = sum((P[i] - C[i])**2 for i in range(3))
        
        prefactor = (2 * pi / zeta) * exp(-alphaA * alphaB * RAB2 / zeta)
        
        T = zeta * RPC2
        L_total = sum(pwrA) + sum(pwrB)
        F = boys_function(L_total, T)
        
        self.memo_v = {}
        self.P = P
        self.A = A
        self.B = B
        self.C = C
        self.zeta = zeta
        self.F = F
        self.prefactor_v = prefactor
        
        return self.vrr_v(pwrA, pwrB, 0) * normA * normB

    def vrr_v(self, a, b, m):
        key = (a, b, m)
        if key in self.memo_v:
            return self.memo_v[key]
            
        if sum(a) == 0 and sum(b) == 0:
            return self.prefactor_v * self.F[m]
            
        if sum(a) > 0:
            i = 0
            if a[1] > 0: i = 1
            if a[2] > 0: i = 2
            
            a_minus = list(a)
            a_minus[i] -= 1
            a_minus = tuple(a_minus)
            
            PA = self.P[i] - self.A[i]
            PC = self.P[i] - self.C[i]
            
            val = PA * self.vrr_v(a_minus, b, m) - PC * self.vrr_v(a_minus, b, m + 1)
            
            if a_minus[i] > 0:
                a_minus_2 = list(a)
                a_minus_2[i] -= 2
                a_minus_2 = tuple(a_minus_2)
                val += (a[i] - 1) / (2 * self.zeta) * (self.vrr_v(a_minus_2, b, m) - self.vrr_v(a_minus_2, b, m + 1))
                
            if b[i] > 0:
                b_minus = list(b)
                b_minus[i] -= 1
                b_minus = tuple(b_minus)
                val += b[i] / (2 * self.zeta) * (self.vrr_v(a_minus, b_minus, m) - self.vrr_v(a_minus, b_minus, m + 1))
                
            self.memo_v[key] = val
            return val
            
        elif sum(b) > 0:
            i = 0
            if b[1] > 0: i = 1
            if b[2] > 0: i = 2
            
            b_minus = list(b)
            b_minus[i] -= 1
            b_minus = tuple(b_minus)
            
            PB = self.P[i] - self.B[i]
            PC = self.P[i] - self.C[i]
            
            val = PB * self.vrr_v(a, b_minus, m) - PC * self.vrr_v(a, b_minus, m + 1)
            
            if b_minus[i] > 0:
                b_minus_2 = list(b)
                b_minus_2[i] -= 2
                b_minus_2 = tuple(b_minus_2)
                val += (b[i] - 1) / (2 * self.zeta) * (self.vrr_v(a, b_minus_2, m) - self.vrr_v(a, b_minus_2, m + 1))
                
            self.memo_v[key] = val
            return val
            
        return 0.0

class Integral:
    def __init__(self):
        self.MAX_ITERATION = 100
        self.EPS = 3.0e-7
        self.SMALL = 0.00000001
        self.FPMIN = 1.0e-30
        self.cof = [76.18009172947146, -86.50532032941677, \
        24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
        self.os_evaluator = OSEvaluator()
        self.one_electron_os = OneElectronOS()
    
    def factorial(self, n):
        value = 1        
        while(n > 1): value = value * n; n-=1
        return value
    
    def factorial2(self, n):
        value = 1        
        while(n > 0): value = value * n; n-=2
        return value
    
    def factorialRatioSquared(self, a, b):
        return (self.factorial(a) / self.factorial(b) / self.factorial(a-2*b))
    
    def binomial(self, i, j):
        return (self.factorial(i) / self.factorial(j) / self.factorial(i-j))

    def mpow(self, x, y):
        if (round(y,10)==0.0): return 1.0
        if (round(x,10)==0.0): return 0.0
        return pow(x, y)
        
    def binomialPrefactor(self, s, ia, ib, xpa, xpb):
        sum = 0.0        
        for t in range(0, s+1):
            if (((s-ia) <= t) and (t <= ib)):
                sum += self.binomial(ia, s-t) * self.binomial(ib, t) \
                        * self.mpow(xpa, ia-s+t) * self.mpow(xpb, ib-t)
        return sum
    
    def dist2(self, a, b):
        x=a[0]-b[0]; y=a[1]-b[1]; z=a[2]-b[2]
        return (x*x+y*y+z*z)

    def dist(self, a, b):
        return (sqrt(self.dist2(a, b)))
    
    def ovrlp(self, alp1, pwr, a, alp2, pwr2, b):
        # Note: original ovrlp didn't take norms, but pg.ovrlp multiplies by norms.
        # Here we need to be careful.
        # The original Integral.ovrlp returned the integral of primitives WITHOUT normalization constants?
        # Let's check pg.ovrlp:
        # return (self.norfac * og.norfac * integral.ovrlp(self.exp, self.pwr, self.org, og.exp, og.pwr, og.org))
        # So Integral.ovrlp computed unnormalized integral.
        # My OneElectronOS.overlap takes norms and returns normalized integral.
        # So I should pass 1.0 as norms here?
        # Yes.
        return self.one_electron_os.overlap(a, 1.0, pwr, alp1, b, 1.0, pwr2, alp2)

    def nuc(self, a, nrm1, pwr1, alp1, b, nrm2, pwr2, alp2, c):
        # Original returned negative value.
        # My OneElectronOS returns positive integral.
        # So I return negative.
        # Also original took norms.
        return -self.one_electron_os.nuclear_attraction(a, nrm1, pwr1, alp1, b, nrm2, pwr2, alp2, c)

    def kinetic(self, alp1, pwr1, a, alp2, pwr2, b):
        # Original didn't take norms.
        return self.one_electron_os.kinetic(a, 1.0, pwr1, alp1, b, 1.0, pwr2, alp2)


        
    def coulombRepulsion(self, a, aNorm, aPower, aAlpha, \
          b, bNorm, bPower, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dPower, dAlpha):
        return self.os_evaluator.run(a, aNorm, aPower, aAlpha, \
                                     b, bNorm, bPower, bAlpha, \
                                     c, cNorm, cPower, cAlpha, \
                                     d, dNorm, dPower, dAlpha)
        
    def coulomb(self, a, b, c, d):
        jij = 0.0        
        aExps = a.exps; aCoefs = a.cofs; aNorms = a.primnrm
        aOrigin = a.org; aPower = a.pwr

        bExps = b.exps; bCoefs = b.cofs; bNorms = b.primnrm
        bOrigin = b.org; bPower = b.pwr

        cExps = c.exps; cCoefs = c.cofs; cNorms = c.primnrm
        cOrigin = c.org; cPower = c.pwr        

        dExps = d.exps; dCoefs = d.cofs; dNorms = d.primnrm
        dOrigin = d.org; dPower = d.pwr
                  
        for i in range(0, len(aExps)):
            iaCoef = aCoefs[i]; iaExp = aExps[i]; iaNorm = aNorms[i]
            for j in range(0, len(bExps)):
                jbCoef = bCoefs[j]; jbExp = bExps[j]; jbNorm = bNorms[j]
                for k in range(0, len(cExps)):
                    kcCoef = cCoefs[k]; kcExp = cExps[k]; kcNorm = cNorms[k]            
                    for l in range(0, len(dExps)):
                        repulsionTerm = self.coulombRepulsion( \
                                        aOrigin, iaNorm, aPower, iaExp, \
                                        bOrigin, jbNorm, bPower, jbExp, \
                                        cOrigin, kcNorm, cPower, kcExp, \
                                        dOrigin, dNorms[l], dPower, dExps[l])                        
                        jij += iaCoef*jbCoef*kcCoef*dCoefs[l]*repulsionTerm
        
        return (a.norfac*b.norfac*c.norfac*d.norfac*jij)
        
    def ijkl2intindex(self, i, j, k, l):
        temp = 0
        
        if (i<j):
            temp=i; i=j; j=temp
        if (k<l):
            temp=k; k=l; l=temp
        
        ij = i*(i+1)/2+j
        kl = k*(k+1)/2+l
        
        if (ij < kl):
            temp=ij; ij=kl; kl=temp
        
        return int(ij * (ij+1) / 2+kl)
    
integral = Integral()
        
class pg:
    def __init__(self, org, pwr, exp, cof):
        self.org=org; self.pwr=pwr
        self.exp=exp; self.cof=cof
        self.norm()

    def ovrlp(self, og):
        return (self.norfac * og.norfac * \
                integral.ovrlp(self.exp, self.pwr, self.org, og.exp, og.pwr, og.org))

    def norm(self):
        l=self.pwr[0]; m=self.pwr[1]; n=self.pwr[2]
        self.norfac = sqrt(pow(2, 2*(l+m+n)+1.5) \
                          * pow(self.exp, l+m+n+1.5) \
                          / integral.factorial2(2*l-1) \
                          / integral.factorial2(2*m-1) \
                          / integral.factorial2(2*n-1) / pow(pi, 1.5))

    def nuc(self, og, cen):
        return (integral.nuc(self.org, self.norfac, self.pwr, \
                self.exp, og.org, og.norfac, og.pwr, og.exp, cen))

    def kinetic(self, og):
        return (self.norfac*og.norfac*integral.kinetic( \
                self.exp, self.pwr, self.org, og.exp, og.pwr, og.org))
        
class cg:
    def __init__(self, org, pwr):
        self.org=org; self.pwr=pwr
        self.pgs = []; self.exps = []; self.cofs = []
        self.norfac = 1.0
        self.primnrm = []

    def addpg(self, exp, cof):        
        self.pgs.append(pg(self.org, self.pwr, exp, cof))
        self.exps.append(exp); self.cofs.append(cof)

    def norm(self):
        self.norfac = 1.0 / sqrt(self.ovrlp(self))        
        for i in range(0, len(self.pgs)):
            self.primnrm.append(self.pgs[i].norfac)

    def ovrlp(self, og):
        sij = 0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                sij += ipg.cof*jpg.cof*ipg.ovrlp(jpg)
        return self.norfac*og.norfac*sij

    def kinetic(self, og):
        tij=0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                tij += ipg.cof*jpg.cof*ipg.kinetic(jpg)
        return self.norfac*og.norfac*tij

    def nuc(self, og, cen):
        vij=0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                vij += ipg.cof*jpg.cof*ipg.nuc(jpg, cen)
        return self.norfac*og.norfac*vij
    
class atmbas:
    def __init__(self, sym, num):
        self.sym=sym; self.num=num;
        self.orbs=[]
                
class basis:
    def __init__(self, mol, bas):
        self.bfs=[]

        for atm in mol.atoms:
            atmbas = curbas[atm[0]]
            for ab in atmbas:
                pwrlst = powers[ab[0]]
                for pwr in pwrlst:
                    c = cg(atm[1], pwr)
                    for orb in ab[1]:
                        c.addpg(orb[0], orb[1])
                    c.norm()
                    self.bfs.append(c)
                
class oneEI:
    def __init__(self, bfs, mol):
        self.H = []; self.S = []        
        nbf = len(bfs)
        for i in range(0, nbf):
            self.H.append([]); self.S.append([])
            for j in range(0, nbf):
                self.H[i].append(0.0);
                self.S[i].append(0.0);
                
        for i in range(0, nbf):
            ibf=bfs[i]
            for j in range(0, nbf):
                jbf=bfs[j]
                self.S[i][j] = ibf.ovrlp(jbf)
                self.H[i][j] = ibf.kinetic(jbf)
                
                for k in range(0, len(mol.atoms)):
                    self.H[i][j] += atmnum[mol.atoms[k][0]] * ibf.nuc(jbf, mol.atoms[k][1])

class twoEI:
    def __init__(self, bfs):
        nbf = len(bfs)
        self.nint = int(nbf*(nbf+1)*(nbf*nbf+nbf+2)/8)
        self.ints = []
        for i in range(0, self.nint): self.ints.append(0.0)

        for i in range(0, nbf):
            bfi = bfs[i]            
            for j in range(0, (i+1)):
                bfj = bfs[j]
                ij = i*(i+1)/2+j                
                for k in range(0, nbf):
                    bfk = bfs[k]                    
                    for l in range(0, (k+1)):
                        bfl = bfs[l]                        
                        kl = k*(k+1)/2+l;
                        if (ij >= kl):
                            ijkl = integral.ijkl2intindex(i, j, k, l)
                            self.ints[ijkl] = integral.coulomb(bfi, bfj, bfk, bfl)

class hf:
    def __init__(self, mol, e1, e2):
        self.e1=e1; self.e2=e2
        ele=0
        for i in range(0, len(mol.atoms)):
            ele += atmnum[mol.atoms[i][0]]
        self.nocc=ele/2

        self.eNuc=0.0
        for i in range(0, len(mol.atoms)):
            ati = mol.atoms[i]
            for j in range(0, i):
                atj = mol.atoms[j]                
                self.eNuc += atmnum[ati[0]] * atmnum[atj[0]] / integral.dist(ati[1], atj[1])
                
    def doRotate(self, a, i, j, k, l, sin, tau):
        g = a[i][j]; h = a[k][l]        
        a[i][j] = g - sin * (h + g*tau)
        a[k][l] = h + sin * (g - h*tau)

    def sortEval(self, eval, evec):
        n = len(evec)
        for i in range(0, n):
            k=i; p=eval[i]
            for j in range(i+1, n):
                if (eval[j] <= p):
                    k=j; p=eval[j]
            if (k!=i):
                eval[k] = eval[i]; eval[i] = p                
                for j in range(0, n):
                    p = evec[j][i]; evec[j][i] = evec[j][k]; evec[j][k] = p
        return eval, evec
                    
    def diagonalize(self, m):
        maxItr=50
        evec=[]; a=[]; n=len(m)
        for i in range(0, n):
            evec.append([]); a.append([])
            for j in range(0, n):
                evec[i].append(0.0)
                a[i].append(m[i][j])
            evec[i][i] = 1.0
        eval=[]; b=[]; z=[]        
        for i in range(0, n):
            eval.append(a[i][i]); b.append(a[i][i]); z.append(0.0)
        for sweeps in range(0, maxItr):
            sum = 0.0
            for i in range(0, n-1):
                for j in range(i+1, n):
                    sum += abs(a[i][j])
            if (sum == 0.0): break
            
            if (sweeps < 3): zeroTol = 0.2*sum/(n*n)
            else: zeroTol = 0.0
            
            for ip in range(0, n-1):
                for iq in range(ip+1, n):
                    g = 100.0 * abs(a[ip][iq])
                    
                    if ((sweeps > 4) and (float(abs(eval[ip])+g)==float(abs(eval[ip]))) \
                        and (float(abs(eval[iq])+g)==float(abs(eval[iq])))): a[ip][iq] = 0.0
                    elif (abs(a[ip][iq]) > zeroTol):
                        h = eval[iq]-eval[ip]                        
                        if (float((abs(h)+g))==float(abs(h))): t = a[ip][iq] / h
                        else:
                            theta = 0.5 * h / a[ip][iq]
                            t = 1.0 / (abs(theta)+sqrt(1.0 + theta*theta))                            
                            if (theta < 0.0): t = -t
                            
                        cos = 1.0 / sqrt(1.0 + t*t)
                        sin = t * cos
                        tau = sin / (1.0 + cos)
                        h  = t * a[ip][iq]                        
                        z[ip] -= h; z[iq] += h
                        eval[ip] -= h; eval[iq] += h                        
                        a[ip][iq] = 0.0
                                                
                        for j in range(0, ip): self.doRotate(a, j, ip, j, iq, sin, tau)
                        for j in range(ip+1, iq): self.doRotate(a, ip, j, j, iq, sin, tau)
                        for j in range(iq+1, n): self.doRotate(a, ip, j, iq, j, sin, tau)
                        for j in range(0, n): self.doRotate(evec, j, ip, j, iq, sin, tau)

            for ip in range(0, n):
                b[ip] += z[ip]; eval[ip] = b[ip]; z[ip] = 0.0
        eval, evec = self.sortEval(eval, evec)
        return eval, self.trans(evec)

    def trans(self, m):
        nm=[]; rc=len(m); cc=len(m[0])
        for i in range(0, cc):
            nm.append([])
            for j in range(0, rc):
                nm[i].append(m[j][i])
        return nm
    
    def mul(self, m1, m2):
        nm=[]; rc=len(m1); brc=len(m2); bcc=len(m2[0])
        for i in range(0, rc):
            nm.append([])
            for j in range(0, bcc):
                cij=0.0
                for k in range(0, brc):
                    cij += m1[i][k]*m2[k][j]
                nm[i].append(cij)
        return nm

    def add(self, m1, m2):
        nm=[]; rc=len(m1);
        for i in range(0, rc):
            nm.append([])
            for j in range(0, rc):
                nm[i].append(m1[i][j]+m2[i][j])
        return nm

    def tr(self, m):
        r=0.0
        for i in range(len(m)): r+=m[i][i]
        return r
    
    def simiTransT(self, m1, m2):
        return self.mul(self.mul(self.trans(m1), m2), m1)
    
    def simiTrans(self, m1, m2):
        return self.mul(self.mul(m1, m2), self.trans(m1))
    
    def symOrth(self, m):
        eval, evec = self.diagonalize(m)
        sHalf=[]; rc=len(m)
        for i in range(0, rc):
            sHalf.append([])
            for j in range(0, rc):
                sHalf[i].append(0.0)
            sHalf[i][i] = 1.0 / sqrt(eval[i])
            
        return (self.simiTransT(evec, sHalf))
        
    def compOrbES(self, h, s):
        x = self.symOrth(s)
        a = self.simiTrans(x, h)
        eval, evec = self.diagonalize(a)
        return eval, self.mul(evec, x)

    def makD(self, mos):
        d = []
        for i in range(0, int(self.nocc)):
            d.append([])
            for j in range(0, len(mos)):
                d[i].append(mos[i][j])
        return self.mul(self.trans(d), d)

    def dot(self, v1, v2):
        res = 0.0
        for i in range(0, len(v1[0])): res += v1[0][i]*v2[0][i]
        return res
    
    def makG(self, mos, D):
        nbf=len(D); G=[]; oneD=[]; oneD.append([]); tmpV=[]; tmpV.append([])
        for i in range(0, nbf):
            G.append([])
            for j in range(0, nbf):
                G[i].append(0.0)
                oneD[0].append(D[i][j])
                tmpV[0].append(0.0)
                
        for i in range(0, nbf):
            for j in range(0, i+1):
                for k in range(0, nbf*nbf): tmpV[0][i]=0.0
                kl = 0                
                for k in range(0, nbf):
                    for l in range(0, nbf):
                        idxJ  = integral.ijkl2intindex(i, j, k, l)
                        idxK1  = integral.ijkl2intindex(i, k, j, l)
                        idxK2  = integral.ijkl2intindex(i, l, k, j)
                        tmpV[0][kl] = 2.0*self.e2.ints[idxJ]-0.5*self.e2.ints[idxK1]-0.5*self.e2.ints[idxK2]
                        kl+=1
                G[i][j] = G[j][i] = self.dot(tmpV, oneD)
                
        return G
    
    def scf(self):                
        maxItr=30; eneTol=1.0e-4
        orbE, mos = self.compOrbES(self.e1.H, self.e1.S)
        oldEnergy=0.0
        for scfIteration in range(0, maxItr):
            D = self.makD(mos);  
            G = self.makG(mos, D)
            F = self.add(self.e1.H, G)
            orbE, mos = self.compOrbES(F, self.e1.S)
            eOne = self.tr(self.mul(D, self.e1.H))
            eTwo = self.tr(self.mul(D, F))
            energy = eOne + eTwo + self.eNuc
            print(scfIteration, energy)
            if (abs(energy - oldEnergy) < eneTol): break            
            oldEnergy = energy
                            
class mobihf:
    def __init__(self):
        #self.filename = "c:\\system\\apps\\python\\my\\h2.inp"
        self.filename = "h2.inp"

    def run(self):
        self.do_run()

    def do_run(self):
        self.read_inp()
        self.chk_inp()
        self.bfs = basis(self.mol, self.basis).bfs
        self.ei1 = oneEI(self.bfs, self.mol)
        self.ei2 = twoEI(self.bfs)
        hf(self.mol, self.ei1, self.ei2).scf()

    def read_inp(self):
        global radii
        """
        xyzf = open(self.filename, "r")
        lines = xyzf.readlines()
        xyzf.close()
        """
        lines = ["2","h2 hf sto-3g", "H 0 0 0", "H 1 0 0"]
        lines = ["3","wat hf sto-3g", "H 0.75 -0.45 0", "O 0 0.11 0", "H -0.75 -0.45 0"]
        
        self.mol = molecule()
        self.mol.noOfAtoms = int(lines[0])
        self.mol.title, self.level, self.basis = lines[1].split()
        self.level = self.level.upper()
        self.basis = self.basis.upper()
        for i in range(2, len(lines)):
            words = lines[i].split()
            self.mol.atoms.append((words[0].upper(), list(map(float, words[1:])), []))
        
        for i in range(0, self.mol.noOfAtoms):
            atm1 = self.mol.atoms[i]
            for j in range(0, i):
                atm2 = self.mol.atoms[j] 
                x = atm1[1][0]-atm2[1][0]
                y = atm1[1][1]-atm2[1][1]
                z = atm1[1][2]-atm2[1][2]
                dist = sqrt(x*x+y*y+z*z)
                radsum = radii[atm1[0]] + radii[atm2[0]]
                if (((radsum - 0.4) < dist) and (dist < (radsum + 0.4))):
                    atm1[2].append(j)
                    atm2[2].append(i)
                    
        lines = None

        print(self.mol.title, self.level, self.basis)        

    def chk_inp(self):
        if (self.level != "HF"):
            print("Unsupported level")
        if (self.basis != "STO-3G"):
            print("Unsupported basis")
            
import time
            
if __name__ == '__main__':
    st = time.time()
    mobihf().run()
    en = time.time()
    print(en-st)
    

