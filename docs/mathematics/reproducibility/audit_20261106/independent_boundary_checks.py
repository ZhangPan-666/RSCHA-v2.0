"""Independent high-precision checks of selected handbook boundary/spectral formulae.

These finite tests supplement the analytic audit; they are not a proof for all
parameters. The formulas under test are evaluated against quadrature,
differentiation or a distinct hypergeometric representation.
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps=45
out=Path(__file__).parent
checks=[]
def check(group,name,lhs,rhs,tol=mp.mpf('1e-32')):
    error=abs(lhs-rhs)/max(1,abs(lhs),abs(rhs))
    assert error<tol,(group,name,mp.nstr(error,8))
    checks.append(dict(group=group,name=name,scaled_error=float(error)))

def P(n,m,t):
    # Negative-order hypergeometric form plus exact integer-order gamma ratio.
    return (-1)**m*mp.rf(n-m+1,2*m)*mp.tan(t/2)**m/mp.factorial(m)*mp.hyp2f1(-n,n+1,m+1,mp.sin(t/2)**2)

for n,m in [(mp.mpf('0.7'),0),(mp.mpf('2.3'),1),(mp.mpf('5.2'),3)]:
    for t in map(mp.mpf,['0.01','.9','2.8']):
        d=mp.diff(lambda x:P(n,m,x),t)
        dd=mp.diff(lambda x:P(n,m,x),t,2)
        check('Legendre','raising',d,m*mp.cot(t)*P(n,m,t)+P(n,m+1,t))
        check('Legendre','ODE',dd+mp.cot(t)*d+(n*(n+1)-m*m/mp.sin(t)**2)*P(n,m,t),0)
        if m:
            check('Legendre','lowering',d,-m*mp.cot(t)*P(n,m,t)-(n+m)*(n-m+1)*P(n,m-1,t))

def K(tau,m,t):
    # A different hypergeometric argument/prefactor from the sine-power form
    # used in the manuscript's explicit normalisation derivation.
    return mp.rf(mp.mpf('.5')+1j*tau-m,2*m)*(-1)**m*mp.tan(t/2)**m/mp.factorial(m)*mp.hyp2f1(mp.mpf('.5')-1j*tau,mp.mpf('.5')+1j*tau,m+1,mp.sin(t/2)**2)

for tau in map(mp.mpf,['.001','.7','7']):
    for m in [0,1,4]:
        for t0 in map(mp.mpf,['.1','1','2.9']):
            integral=mp.quad(lambda t:K(tau,m,t)**2*mp.sin(t),[0,t0])
            val=K(tau,m,t0)
            theta_d=mp.diff(lambda t:K(tau,m,t),t0)
            tau_d=mp.diff(lambda a:K(a,m,t0),tau)
            mixed=mp.diff(lambda a:mp.diff(lambda t:K(a,m,t),t0),tau)
            closed=mp.sin(t0)/(2*tau)*(val*mixed-tau_d*theta_d)
            check('Mehler','parameter-norm',integral,closed)

for m in [1,2,4]:
    for t0 in map(mp.mpf,['.1','1','2.9']):
        u=mp.tan(t0/2)
        direct=mp.quad(lambda t:(mp.gamma(m)*mp.tan(t/2)**m)**2*mp.sin(t),[0,t0])
        lerch=2*u**(2*m)*mp.gamma(m)**2*(m*mp.lerchphi(-u*u,1,m)-1/(1+u*u))
        check('Constant','Lerch-norm',direct,lerch)
        eps=mp.mpf('1e-12')
        check('Constant','regularised-limit',P(-eps,m,t0)/eps,mp.gamma(m)*u**m,mp.mpf('1e-9'))

for s in map(mp.mpf,['.01','.7','5']):
    a=mp.mpf(3); ref=mp.mpf('2.3'); b=a*mp.exp(s)
    def Rp(p,r):
        t=mp.log(r/a)
        return mp.sqrt(2/(s*((2*p*mp.pi)**2+s*s)))*mp.sqrt(a/r)*(2*p*mp.pi*mp.cos(p*mp.pi*t/s)+s*mp.sin(p*mp.pi*t/s))
    for p in [1,2,4]:
        check('Radial','Neumann-a',mp.diff(lambda r:Rp(p,r),a),0)
        check('Radial','Neumann-b',mp.diff(lambda r:Rp(p,r),b),0)
        check('Radial','norm',mp.quad(lambda t:Rp(p,a*mp.exp(t))**2*a*mp.exp(t),[0,s]),a)
        for n in [0,2]:
            Ap=mp.sqrt(2/(s*((2*p*mp.pi)**2+s*s)))
            D=(n+mp.mpf('.5'))**2+(p*mp.pi/s)**2
            lhs=mp.quad(lambda t:(ref/(a*mp.exp(t)))**(n+1)*Rp(p,a*mp.exp(t))*a*mp.exp(t),[0,s])
            rhs=2*p*mp.pi*Ap*(n+1)*(1-(-1)**p*mp.exp(-(n+mp.mpf('.5'))*s))/D*ref**(n+1)/a**n
            check('Radial','internal-overlap',lhs,rhs)
            lhs=mp.quad(lambda t:(a*mp.exp(t)/ref)**n*Rp(p,a*mp.exp(t))*a*mp.exp(t),[0,s])
            rhs=2*p*mp.pi*Ap*n*((-1)**p*mp.exp((n+mp.mpf('.5'))*s)-1)/D*a**(n+1)/ref**n
            check('Radial','external-overlap',lhs,rhs)

# Explicit R-SCHA2D nonuniqueness despite decay at infinity, when only the
# displayed normal trace on r=1 and the north component on its rim are given.
t0=mp.mpf('.7')
def continuation(tau,r,t):
    R=r**(-mp.mpf('.5'))*(mp.cos(tau*mp.log(r))+mp.sin(tau*mp.log(r))/(2*tau))
    return R*K(tau,0,t)/mp.diff(lambda z:K(tau,0,z),t0)
for tau in map(mp.mpf,['.7','1.3']):
    for t in [mp.mpf('.2'),mp.mpf('.5')]:
        check('R-SCHA2D','normal-trace',mp.diff(lambda r:continuation(tau,r,t),mp.mpf(1)),0)
    check('R-SCHA2D','rim-north',mp.diff(lambda t:continuation(tau,mp.mpf(1),t),t0),1)
distinct=abs(continuation(mp.mpf('.7'),mp.mpf(2),mp.mpf('.3'))-continuation(mp.mpf('1.3'),mp.mpf(2),mp.mpf('.3')))
assert distinct>mp.mpf('.01')
# One-radius branch recovery when BOTH surface-potential and normal coefficients
# are provided, contrary to a blanket claim of radial-branch nonuniqueness.
for n in map(mp.mpf,['0','.7','4.2']):
    gi=mp.mpf('1.2');ge=mp.mpf('-.3');c=gi+ge;z=-(n+1)*gi+n*ge
    recovered_i=(n*c-z)/(2*n+1);recovered_e=((n+1)*c+z)/(2*n+1)
    check('R-SCHA2D','two-trace-internal',recovered_i,gi)
    check('R-SCHA2D','two-trace-external',recovered_e,ge)

report=dict(status='PASS',precision=mp.mp.dps,checks=len(checks),maximum_scaled_error=max(c['scaled_error'] for c in checks),
            groups={g:sum(c['group']==g for c in checks) for g in sorted({c['group'] for c in checks})},
            counterexample_distinct_continuations_at_r2_theta03=float(distinct),
            scope='Finite parameter samples; spectral and function-space assertions were audited analytically, not proved by these tests.',
            details=checks)
(out/'independent_boundary_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='details'},indent=2))
