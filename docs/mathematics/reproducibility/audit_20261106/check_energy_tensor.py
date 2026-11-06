"""Independent 45-digit formula audit, with explicit finite parameter scope.

The test constructs regular angular functions from a hypergeometric formula,
normalises by quadrature, and compares differentiated/integrated potentials
with the handbook's closed volume/surface factors. Tensor checks instead
start from Cartesian harmonic polynomials and differentiate in Cartesian
coordinates before rotating into NED. No old verification script is imported.
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps=45
checks=[]

def record(name,actual,expected,tol=mp.mpf('1e-32')):
    error=abs(actual-expected)/(1+abs(actual)+abs(expected))
    checks.append(dict(name=name,scaled_error=float(error)))
    if error>tol:
        raise AssertionError((name,mp.nstr(actual,12),mp.nstr(expected,12),mp.nstr(error,8)))

def regular(nu,m,t):
    return mp.sin(t)**m*mp.hyp2f1(m-nu,m+nu+1,m+1,mp.sin(t/2)**2)

def normalise(f,t0):
    norm=mp.sqrt(mp.quad(lambda t:f(t)**2*mp.sin(t),[0,t0/2,t0]))
    return lambda t:f(t)/norm

def dot(a,b):return sum(x*y for x,y in zip(a,b))

# These independently specified harmonic Cartesian functions exercise every
# spherical mixed derivative and all six NED Hessian components.
potentials=[lambda x,y,z:z,lambda x,y,z:x*y,lambda x,y,z:x*y*z,
            lambda x,y,z:x*x-y*y,lambda x,y,z:z*z-(x*x+y*y)/2,
            lambda x,y,z:1/mp.sqrt(x*x+y*y+z*z),
            lambda x,y,z:z/(x*x+y*y+z*z)**mp.mpf('1.5')]
for pi,vc in enumerate(potentials):
    for r,t,l in [(mp.mpf('1.7'),mp.mpf('.41'),mp.mpf('.73')),
                  (mp.mpf('2.2'),mp.mpf('1.23'),mp.mpf('-1.04'))]:
        xyz=[r*mp.sin(t)*mp.cos(l),r*mp.sin(t)*mp.sin(l),r*mp.cos(t)]
        def spherical(rr,tt,ll):
            return vc(rr*mp.sin(tt)*mp.cos(ll),rr*mp.sin(tt)*mp.sin(ll),rr*mp.cos(tt))
        args=(r,t,l)
        v_r=mp.diff(spherical,args,(1,0,0));v_t=mp.diff(spherical,args,(0,1,0));v_l=mp.diff(spherical,args,(0,0,1))
        v_rr=mp.diff(spherical,args,(2,0,0));v_tt=mp.diff(spherical,args,(0,2,0));v_ll=mp.diff(spherical,args,(0,0,2))
        v_rt=mp.diff(spherical,args,(1,1,0));v_rl=mp.diff(spherical,args,(1,0,1));v_tl=mp.diff(spherical,args,(0,1,1))
        formula=[[-v_tt/r**2-v_r/r,v_tl/(r*r*mp.sin(t))-mp.cos(t)*v_l/(r*r*mp.sin(t)**2),-v_rt/r+v_t/r**2],
                 [0,-v_ll/(r*r*mp.sin(t)**2)-v_r/r-mp.cot(t)*v_t/r**2,v_rl/(r*mp.sin(t))-v_l/(r*r*mp.sin(t))],
                 [0,0,-v_rr]]
        basis=[[-mp.cos(t)*mp.cos(l),-mp.cos(t)*mp.sin(l),mp.sin(t)],
               [-mp.sin(l),mp.cos(l),0],[-mp.sin(t)*mp.cos(l),-mp.sin(t)*mp.sin(l),-mp.cos(t)]]
        h=[]
        for i in range(3):
            row=[]
            for j in range(3):
                order=[0,0,0];order[i]+=1;order[j]+=1
                row.append(-mp.diff(vc,tuple(xyz),tuple(order)))
            h.append(row)
        for i in range(3):
            for j in range(i,3):
                transformed=sum(basis[i][a]*h[a][b]*basis[j][b] for a in range(3) for b in range(3))
                record(f'Cartesian-to-NED Hessian {pi}/{float(t)}/{i}{j}',formula[i][j],transformed)
        grad=[-mp.diff(vc,tuple(xyz),tuple(int(k==i) for k in range(3))) for i in range(3)]
        b_sph=[v_t/r,-v_l/(r*mp.sin(t)),v_r]
        for i in range(3):record(f'Cartesian-to-NED field {pi}/{float(t)}/{i}',b_sph[i],dot(grad,basis[i]))
        record(f'Harmonic Hessian trace {pi}/{float(t)}',sum(formula[i][i] for i in range(3)),0)

ra=mp.mpf('1.3');rb=mp.mpf('2.4');beta=mp.mpf('1.8');s=mp.log(rb/ra)
for t0 in [mp.pi/4,mp.pi/2]:
    for m in [0,1]:
        # Use an eigenvalue with an exact Dirichlet cap zero when t0=pi/2;
        # otherwise find the first regular Dirichlet degree by bracketing.
        if t0==mp.pi/2:nu=mp.mpf(m+1)
        else:
            step=mp.mpf('.5');lo=mp.mpf(m);fl=regular(lo,m,t0)
            hi=lo+step
            while regular(hi,m,t0)*fl>0:hi+=step
            nu=mp.findroot(lambda n:regular(n,m,t0),(hi-step,hi),solver='anderson')
        P=normalise(lambda t:regular(nu,m,t),t0)
        ell=mp.pi if m else 2*mp.pi
        modes=[]
        for family,U in [('Int',lambda r:(ra/r)**(nu+1)),('Ext',lambda r:(r/ra)**nu)]:
            modes.append((family,U,P,nu*(nu+1)))
        for p in [1,2]:
            tau=p*mp.pi/s;mu=mp.mpf('.25')+tau*tau
            K=normalise(lambda t,tau=tau:mp.re(regular(-mp.mpf('.5')+1j*tau,m,t)),t0)
            factor=mp.sqrt(2/(s*((2*p*mp.pi)**2+s*s)))
            U=lambda r,p=p,factor=factor:(ra/r)**mp.mpf('.5')*factor*(2*p*mp.pi*mp.cos(p*mp.pi*mp.log(r/ra)/s)+s*mp.sin(p*mp.pi*mp.log(r/ra)/s))
            modes.append((f'Meh{p}',U,K,-mu))
        C=normalise(lambda t:mp.tan(t/2)**m,t0)
        R0=mp.sqrt(ra/(rb-ra))
        modes.append(('Cons',lambda r:R0,C,mp.mpf(0)))
        angular={}
        def angular_products(i,j):
            if (i,j) in angular:return angular[i,j]
            A=modes[i][2];B=modes[j][2]
            i0=mp.quad(lambda t:A(t)*B(t)*mp.sin(t),[0,t0/2,t0])
            i1=mp.quad(lambda t:(mp.diff(A,t)*mp.diff(B,t)+(m*m/mp.sin(t)**2*A(t)*B(t) if m else 0))*mp.sin(t),[0,t0/2,t0])
            angular[i,j]=i0,i1
            return i0,i1
        for i,(family,U,A,lam) in enumerate(modes):
            i0,i1=angular_products(i,i)
            jr=mp.quad(lambda r:r*r*mp.diff(U,r)**2,[ra,rb])
            j0=mp.quad(lambda r:U(r)**2,[ra,rb])
            measured=ra*ra*ell*(jr*i0+j0*i1)
            if family=='Int':expected=ra**3*ell*(nu+1)*(1-mp.exp(-(2*nu+1)*s))
            elif family=='Ext':expected=ra**3*ell*nu*(mp.exp((2*nu+1)*s)-1)
            else:expected=ra**3*ell*mp.sin(t0)*A(t0)*mp.diff(A,t0)
            record(f'Volume diagonal {float(t0)}/{m}/{family}',measured,expected)
            if m==0:
                mean=mp.quad(lambda t:A(t)*mp.sin(t),[0,t0])
                flux=ra**3*mp.diff(U,ra)*mean-ra*rb*rb*mp.diff(U,rb)*mean-ra*mp.sin(t0)*mp.diff(A,t0)*mp.quad(U,[ra,rb])
                record(f'Closed-shell flux {float(t0)}/{family}',2*mp.pi*flux,0)
        for i in range(len(modes)):
            for j in range(i+1,len(modes)):
                fi,U,A,li=modes[i];fj,W,B,lj=modes[j]
                i0,i1=angular_products(i,j)
                measured_volume=ra*ra*ell*(mp.quad(lambda r:r*r*mp.diff(U,r)*mp.diff(W,r),[ra,rb])*i0+mp.quad(lambda r:U(r)*W(r),[ra,rb])*i1)
                record(f'Volume off-diagonal {float(t0)}/{m}/{fi}/{fj}',measured_volume,0)
                if fi in ['Int','Ext'] and fj.startswith('Meh'):angular_expected=lj*i0
                elif fi in ['Int','Ext'] and fj=='Cons':angular_expected=0
                elif fi.startswith('Meh'):
                    angular_expected=lj*i0+mp.sin(t0)*A(t0)*mp.diff(B,t0)
                else:angular_expected=i1
                actual=ra*ra*ell*(beta*beta*mp.diff(U,beta)*mp.diff(W,beta)*i0+U(beta)*W(beta)*i1)
                expected=ra*ra*ell*(beta*beta*mp.diff(U,beta)*mp.diff(W,beta)*i0+U(beta)*W(beta)*angular_expected)
                record(f'Fixed-surface cross block {float(t0)}/{m}/{fi}/{fj}',actual,expected)

# H^1 does not control the boundary full-gradient norm. For L2-unit spherical
# harmonics Y_n on the full sphere and V_N=sum n^(-3/2) r^n Y_n, the volume
# energies sum to sum n^(-2), while the surface gradient energies sum to
# sum (2/n+1/n^2). This is an exact orthogonal-series counterexample.
trace_example=[]
for N in [10,100,1000,10000]:
    trace_example.append(dict(N=N,volume_energy=float(mp.fsum(mp.mpf(1)/(n*n) for n in range(1,N+1))),
                              outer_surface_energy=float(mp.fsum(2/mp.mpf(n)+1/mp.mpf(n*n) for n in range(1,N+1)))))
report=dict(status='PASS',precision_digits=mp.mp.dps,comparisons=len(checks),
            maximum_scaled_error=max(c['scaled_error'] for c in checks),checks=checks,
            h1_without_l2_gradient_trace_example=trace_example,
            scope='Seven Cartesian harmonic functions at two off-axis points; finite energy/flux blocks at two cap angles and two azimuthal orders. Not an arbitrary-parameter proof.')
path=Path(__file__).with_name('energy_tensor_checks.json')
path.write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
