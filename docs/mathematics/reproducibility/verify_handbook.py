"""Independent equation checks for the October 2026 handbook revision.

Run from this directory or any working directory: python verify_handbook.py.
Requires numpy, scipy and mpmath. This is a finite numerical cross-check,
not a proof of every theorem or a high-degree stress test of the C++ runtime.
"""
from pathlib import Path
import sys, json, math
from collections import defaultdict
import numpy as np
import mpmath as mp
from scipy.special import lpmv, roots_legendre

ROOT=Path(__file__).resolve().parents[1]
if not (ROOT/'figures'/'code'/'basis.py').exists():
    ROOT=Path(r'E:\OneDrive\学习\参考文献\修正球冠谐分析\Mathematical Fundations of RSCHA')
sys.path.insert(0,str(ROOT/'figures'/'code'))
import basis as plotted
mp.mp.dps=45
results=defaultdict(list)
def check(group, actual, expected, tol=1e-11):
    a=np.asarray(actual,dtype=complex); b=np.asarray(expected,dtype=complex)
    error=float(np.max(np.abs(a-b))/max(1.,float(np.max(np.abs(b)))))
    results[group].append(error)
    if not np.isfinite(error) or error>tol:
        raise AssertionError((group,error,tol,actual,expected))

def P(nu,m,theta):
    # Independent arbitrary-precision Ferrers representation, Gamma ratio
    # expressed as a finite product to include integer-degree zero cases.
    product=mp.fprod([nu-m+1+j for j in range(2*m)])
    return (-1)**m*product/(2**m*mp.factorial(m))*mp.sin(theta)**m*mp.hyp2f1(m-nu,nu+m+1,m+1,mp.sin(theta/2)**2)

def sch(nu,m):
    return mp.mpf(1) if m==0 else (-1)**m*mp.sqrt(2*mp.gamma(nu-m+1)/mp.gamma(nu+m+1))

print('Checking Legendre identities...',flush=True)
for m in range(5):
    for offset in ['2.37','5.61']:
        nu=m+mp.mpf(offset)
        for theta in map(mp.mpf,['0.11','0.37','0.81']):
            x=mp.cos(theta); y=P(nu,m,theta)
            dy=mp.diff(lambda th:P(nu,m,th),theta)
            ddy=mp.diff(lambda th:P(nu,m,th),theta,2)
            degree=(nu-m+1)*P(nu+1,m,theta)-(2*nu+1)*x*y+(nu+m)*P(nu-1,m,theta)
            check('Legendre degree recurrence',degree,0)
            check('Legendre first and second derivatives',m*mp.cot(theta)*y+P(nu,m+1,theta),dy)
            d2=(-m/mp.sin(theta)**2+m*m*mp.cot(theta)**2-(nu-m)*(nu+m+1))*y-mp.cot(theta)*P(nu,m+1,theta)
            check('Legendre first and second derivatives',d2,ddy)
            yn=lambda v:sch(v,m)*P(v,m,theta)
            normrec=(nu-m+1)*sch(nu,m)/sch(nu+1,m)*yn(nu+1)-(2*nu+1)*x*yn(nu)+(nu+m)*sch(nu,m)/sch(nu-1,m)*yn(nu-1)
            check('Schmidt-normalised recurrences',normrec,0)
            if m:
                order=P(nu,m+1,theta)+2*m*x/mp.sin(theta)*y+(nu+m)*(nu-m+1)*P(nu,m-1,theta)
                check('Legendre order recurrence (m>=1)',order,0)
                nr=sch(nu,m)/sch(nu,m+1)*(sch(nu,m+1)*P(nu,m+1,theta))+2*m*x/mp.sin(theta)*yn(nu)+(nu+m)*(nu-m+1)*sch(nu,m)/sch(nu,m-1)*(sch(nu,m-1)*P(nu,m-1,theta))
                check('Schmidt-normalised recurrences',nr,0)
                check('Legendre lowering derivative (m>=1)',-m*mp.cot(theta)*y-(nu+m)*(nu-m+1)*P(nu,m-1,theta),dy)
            tangent=(-1)**m*mp.fprod([nu-m+1+j for j in range(2*m)])/mp.factorial(m)*mp.tan(theta/2)**m*mp.hyp2f1(-nu,nu+1,m+1,mp.sin(theta/2)**2)
            check('Equivalent Ferrers hypergeometric forms',tangent,y)

t0=mp.pi/6
print('Checking cap norms and angular Green identities...',flush=True)
for m in range(5):
    nu=mp.findroot(lambda v:P(v,m,t0), (plotted.degree(m,m)-.01,plotted.degree(m,m)+.01))
    x0=mp.cos(t0)
    pnu=mp.diff(lambda v:P(v,m,t0),nu)
    px=-mp.diff(lambda th:P(nu,m,th),t0)/mp.sin(t0)
    norm=mp.quad(lambda th:P(nu,m,th)**2*mp.sin(th),[0,t0])
    check('Legendre root parameter norm',-(1-x0*x0)/(2*nu+1)*pnu*px,norm)
    for n in [m,m+2,m+5]:
        integ=mp.quad(lambda th:P(n,m,th)*P(nu,m,th)*mp.sin(th),[0,t0])
        boundary=(1-x0*x0)*P(n,m,t0)*px/(nu*(nu+1)-n*(n+1))
        check('SH to cap angular overlap',boundary,integ)
    for tau in map(mp.mpf,['1.7','4.1']):
        K=lambda ta,th:P(-mp.mpf('.5')+1j*ta,m,th)
        y=K(tau,t0); yt=mp.diff(lambda th:K(tau,th),t0)
        ytau=mp.diff(lambda ta:K(ta,t0),tau)
        ytaut=mp.diff(lambda ta:mp.diff(lambda th:K(ta,th),t0),tau)
        norm=mp.quad(lambda th:K(tau,th)**2*mp.sin(th),[0,t0])
        check('Mehler parameter norm',mp.sin(t0)/(2*tau)*(y*ytaut-ytau*yt),norm)
        F=lambda order,ta:mp.hyp2f1(order+mp.mpf('.5')-1j*ta,order+mp.mpf('.5')+1j*ta,order+1,mp.sin(t0/2)**2)
        H=F(m,tau)*F(m+1,tau)*2*tau/(tau*tau+(m+mp.mpf('.5'))**2)+F(m,tau)*mp.diff(lambda ta:F(m+1,ta),tau)-F(m+1,tau)*mp.diff(lambda ta:F(m,ta),tau)
        pref=(-1)**m*mp.sin(t0)**m/(2**m*mp.factorial(m))*mp.fprod([-mp.mpf('.5')+1j*tau-m+1+j for j in range(2*m)])
        explicit=(mp.sin(t0)*pref)**2*(tau*tau+(m+mp.mpf('.5'))**2)*H/(4*tau*(m+1))
        check('Mehler explicit H norm',explicit,norm)
        alpha=tau+mp.mpf('1.3')
        cross=mp.quad(lambda th:K(tau,th)*K(alpha,th)*mp.sin(th),[0,t0])
        crossb=mp.sin(t0)/(tau*tau-alpha*alpha)*(K(alpha,t0)*yt-y*mp.diff(lambda th:K(alpha,th),t0))
        check('Mehler bilinear boundary integral',crossb,cross)

for m in range(1,5):
    for aperture in [t0,mp.mpf('2.1')]:
        u=mp.tan(aperture/2)
        norm=mp.quad(lambda th:(mp.gamma(m)*mp.tan(th/2)**m)**2*mp.sin(th),[0,aperture])
        lerch=2*u**(2*m)*mp.gamma(m)**2*(m*mp.lerchphi(-u*u,1,m)-1/(1+u*u))
        check('Constant Lerch norm including large cap',lerch,norm)
    theta=mp.mpf('.43');eps=mp.mpf('1e-20')
    check('Constant regularised degree-zero branch',P(-eps,m,theta)/eps,mp.gamma(m)*mp.tan(theta/2)**m)

def radial(p,r,ra,rb):
    s=mp.log(rb/ra);t=mp.log(r/ra)
    if p==0:return mp.sqrt(ra/(rb-ra))
    ap=mp.sqrt(2/(s*((2*p*mp.pi)**2+s*s)))
    return ap*mp.sqrt(ra/r)*(2*p*mp.pi*mp.cos(p*mp.pi*t/s)+s*mp.sin(p*mp.pi*t/s))

print('Checking radial norms, derivatives and mapping...',flush=True)
for ra,rb,ref in [(mp.mpf('1.3'),mp.mpf('2.6'),mp.mpf('1.9')),(mp.mpf('2.7'),mp.mpf('3.4'),mp.mpf('2.3'))]:
    s=mp.log(rb/ra)
    for p in range(5):
        fun=lambda r:radial(p,r,ra,rb)
        check('Radial Neumann endpoint conditions',[mp.diff(fun,ra),mp.diff(fun,rb)],[0,0])
        check('Radial norms',mp.quad(lambda r:fun(r)**2,[ra,rb]),ra)
        if p:
            mu=mp.mpf('.25')+(p*mp.pi/s)**2
            r=ra+mp.mpf('.37')*(rb-ra)
            check('Radial ODE and second derivative',mp.diff(fun,r,2)+2/r*mp.diff(fun,r)+mu/r**2*fun(r),0)
            check('Radial zero mean',mp.quad(fun,[ra,rb]),0)
        for n in [0,1,4]:
            iint=mp.quad(lambda r:(ref/r)**(n+1)*fun(r),[ra,rb])
            eint=mp.quad(lambda r:(r/ref)**n*fun(r),[ra,rb])
            if p:
                ap=mp.sqrt(2/(s*((2*p*mp.pi)**2+s*s)));den=(n+mp.mpf('.5'))**2+(p*mp.pi/s)**2
                ic=2*p*mp.pi*ap*(n+1)*(1-(-1)**p*mp.exp(-(n+mp.mpf('.5'))*s))/den*ref**(n+1)/ra**n
                ec=2*p*mp.pi*ap*n*((-1)**p*mp.exp((n+mp.mpf('.5'))*s)-1)/den*ra**(n+1)/ref**n
            else:
                r0=fun(ra)
                ic=r0*ref*s if n==0 else r0*ref*(ref/ra)**n*(1-mp.exp(-n*s))/n
                ec=r0*ra*(ra/ref)**n*(mp.exp((n+1)*s)-1)/(n+1)
            check('Closed SH radial overlaps',ic,iint)
            check('Closed SH radial overlaps',ec,eint)
    for p in range(5):
        for q in range(p):
            check('Radial cross orthogonality',mp.quad(lambda r:radial(p,r,ra,rb)*radial(q,r,ra,rb),[ra,rb]),0)
    mode2d=lambda r:mp.sqrt(ra/r)*(mp.log(r/ra)+2)/2
    check('One-surface repeated-root condition',mp.diff(mode2d,ra),0)
    r=ra+mp.mpf('.4')*(rb-ra)
    check('One-surface repeated-root ODE',r*r*mp.diff(mode2d,r,2)+2*r*mp.diff(mode2d,r)+mode2d(r)/4,0)
    for nu in map(mp.mpf,['3.7','9.1']):
        a,b=mp.mpf('.7'),mp.mpf('-.2');alpha=(ra/rb)**(nu+2);beta=(rb/ra)**(nu-1)
        ba=(nu+1)*a-nu*b;bb=-(nu+1)*alpha*a+nu*beta*b
        check('Two-boundary coefficient solution',[(beta*ba+bb)/((nu+1)*(beta-alpha)),(alpha*ba+bb)/(nu*(beta-alpha))],[a,b])

print('Checking Cartesian field, Hessian and rotations...',flush=True)
def xyzpotential(x,y,z):
    return x**3-3*x*y*y+mp.mpf('.3')*(2*z*z-x*x-y*y)+mp.mpf('.2')*z*(x*x-y*y)+mp.mpf('.1')/mp.sqrt(x*x+y*y+z*z)
def spherical(r,t,l):return xyzpotential(r*mp.sin(t)*mp.cos(l),r*mp.sin(t)*mp.sin(l),r*mp.cos(t))
for r,t,l in [(mp.mpf('1.4'),mp.mpf('.2'),mp.mpf('.7')),(mp.mpf('2.1'),mp.mpf('.8'),mp.mpf('1.3')),(mp.mpf('1.7'),mp.mpf('1.2'),mp.mpf('2.4'))]:
    st,ct,sl,cl=map(float,[mp.sin(t),mp.cos(t),mp.sin(l),mp.cos(l)])
    Q=np.array([[st*cl,ct*cl,-sl],[st*sl,ct*sl,cl],[ct,-st,0.]])
    coords=(r*mp.sin(t)*mp.cos(l),r*mp.sin(t)*mp.sin(l),r*mp.cos(t))
    Bxyz=-np.array([complex(mp.diff(xyzpotential,coords,tuple(1 if i==j else 0 for i in range(3)))) for j in range(3)])
    Hxyz=np.array([[complex(mp.diff(xyzpotential,coords,tuple((1 if k==i else 0)+(1 if k==j else 0) for k in range(3)))) for j in range(3)] for i in range(3)])
    vr=mp.diff(spherical,(r,t,l),(1,0,0));vt=mp.diff(spherical,(r,t,l),(0,1,0));vl=mp.diff(spherical,(r,t,l),(0,0,1))
    check('Spherical field vs Cartesian derivatives',[-vr,-vt/r,-vl/(r*mp.sin(t))],Q.T@Bxyz)
    hrr=mp.diff(spherical,(r,t,l),(2,0,0))
    hrt=mp.diff(spherical,(r,t,l),(1,1,0))/r-vt/r**2
    hrl=mp.diff(spherical,(r,t,l),(1,0,1))/(r*mp.sin(t))-vl/(r*r*mp.sin(t))
    htt=vr/r+mp.diff(spherical,(r,t,l),(0,2,0))/r**2
    htl=mp.diff(spherical,(r,t,l),(0,1,1))/(r*r*mp.sin(t))-mp.cos(t)*vl/(r*r*mp.sin(t)**2)
    hll=vr/r+mp.cot(t)*vt/r**2+mp.diff(spherical,(r,t,l),(0,0,2))/(r*r*mp.sin(t)**2)
    hs=np.array([[complex(hrr),complex(hrt),complex(hrl)],[complex(hrt),complex(htt),complex(htl)],[complex(hrl),complex(htl),complex(hll)]])
    check('Spherical Hessian vs Cartesian derivatives',hs,Q.T@Hxyz@Q)
    check('Harmonic tensor trace',np.trace(hs),0)
    S=np.array([[0.,-1,0],[0,0,1],[-1,0,0]])
    check('NED rotation orthogonality',S@Q.T@Q@S.T,np.eye(3))
    T=-Hxyz.real;R=S@Q.T
    pairs=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    M=np.array([[R[i,k]*R[j,k] if k==z else R[i,k]*R[j,z]+R[i,z]*R[j,k] for k,z in pairs] for i,j in pairs])
    check('Six-component tensor rotation',M@np.array([T[i,j] for i,j in pairs]),np.array([(R@T@R.T)[i,j] for i,j in pairs]))

def rz(a):c,s=np.cos(a),np.sin(a);return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def ry(a):c,s=np.cos(a),np.sin(a);return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Y(ell,m,v):
    x,y,z=v;lon=np.arctan2(y,x);k=abs(m)
    yv=np.sqrt((2*ell+1)/(4*np.pi)*math.factorial(ell-k)/math.factorial(ell+k))*lpmv(k,ell,z)*np.exp(1j*k*lon)
    return yv if m>=0 else (-1)**k*np.conj(yv)
a,b,g=.31,.67,-.24;A=rz(a)@ry(b)@rz(g)
for ell in [1,2,4]:
    D=np.zeros((2*ell+1,2*ell+1),complex)
    for i,m in enumerate(range(-ell,ell+1)):
        for j,q in enumerate(range(-ell,ell+1)):
            pref=math.sqrt(math.factorial(ell+m)*math.factorial(ell-m)*math.factorial(ell+q)*math.factorial(ell-q))
            d=sum((-1)**(t-m+q)*np.cos(b/2)**(2*ell+m-q-2*t)*np.sin(b/2)**(q-m+2*t)/(math.factorial(ell+m-t)*math.factorial(t)*math.factorial(q-m+t)*math.factorial(ell-q-t)) for t in range(max(0,m-q),min(ell+m,ell-q)+1))*pref
            D[i,j]=np.exp(-1j*m*g)*d*np.exp(-1j*q*a)
    for v in [np.array([.2,.3,.8]),np.array([.7,-.4,.1])]:
        v=v/np.linalg.norm(v)
        check('Wigner defining rotation identity',D@np.array([Y(ell,m,v) for m in range(-ell,ell+1)]),np.array([Y(ell,m,A.T@v) for m in range(-ell,ell+1)]))
    check('Wigner unitarity',D@D.conj().T,np.eye(2*ell+1))

print('Checking full volume/surface energy blocks...',flush=True)
z,w=roots_legendre(160);th=(z+1)*float(t0)/2;wt=w*float(t0)/2*np.sin(th)
rr=(z+1)/2+1;wr=w/2;ra,rb,beta=1.,2.,1.37;s=np.log(2.)
for m in range(3):
    modes=[('Int',m),('Int',m+1),('Ext',m),('Ext',m+1),('Meh',1),('Meh',2),('Cons',0)]
    angular=[];boundary=[];rv=[];rd=[];rbv=[];rbd=[];nus=[]
    for fam,i in modes:
        family={'Int':'legendre','Ext':'legendre','Meh':'mehler','Cons':'constant'}[fam]
        ar,at,_=plotted.basis(family,i,m,th)
        norm=np.sqrt(np.dot(wt,ar*ar));ar/=norm;at/=norm
        ay,ady,_=plotted.basis(family,i,m,np.array([float(t0)]))
        boundary.append((ay[0]/norm,ady[0]/norm));angular.append((ar,at))
        if fam in ['Int','Ext']:
            nu=plotted.degree(i,m);power=-nu-1 if fam=='Int' else nu
            u=rr**power;d=power*rr**(power-1);ub=beta**power;db=power*beta**(power-1)
        else:
            nu=None;u,d,_=plotted.radial(i,rr);ub,db,_=[v.item() for v in plotted.radial(i,np.array(beta))]
        nus.append(nu);rv.append(u);rd.append(d);rbv.append(ub);rbd.append(db)
    angular=np.array(angular);rv=np.array(rv);rd=np.array(rd);rbv=np.array(rbv);rbd=np.array(rbd);boundary=np.array(boundary)
    I0=(angular[:,0]*wt)@angular[:,0].T
    I1=(angular[:,1]*wt)@angular[:,1].T+(angular[:,0]*(wt*m*m/np.sin(th)**2))@angular[:,0].T
    Jr=(rd*(wr*rr**2))@rd.T;J0=(rv*wr)@rv.T
    direct_volume=ra*ra*(Jr*I0+J0*I1)
    direct_surface=ra*ra*(beta*beta*np.outer(rbd,rbd)*I0+np.outer(rbv,rbv)*I1)
    analytic_volume=np.zeros((7,7));analytic_surface=np.zeros((7,7))
    for i,(fam,k) in enumerate(modes):
        if fam=='Int':e=(nus[i]+1)*(1-np.exp(-(2*nus[i]+1)*s))
        elif fam=='Ext':e=nus[i]*(np.exp((2*nus[i]+1)*s)-1)
        else:e=np.sin(float(t0))*boundary[i,0]*boundary[i,1]
        analytic_volume[i,i]=e
        for j,(other,q) in enumerate(modes):
            if fam in ['Int','Ext'] and other in ['Int','Ext']:
                val=0. if i!=j else (nus[i]+1 if fam=='Int' else nus[i])*(2*nus[i]+1)*rbv[i]**2
            elif (fam in ['Int','Ext'] and other=='Cons') or (other in ['Int','Ext'] and fam=='Cons'):val=0.
            elif fam in ['Int','Ext'] and other=='Meh':
                mu=.25+(q*np.pi/s)**2;val=(beta*beta*rbd[i]*rbd[j]-mu*rbv[i]*rbv[j])*I0[i,j]
            elif other in ['Int','Ext'] and fam=='Meh':
                mu=.25+(k*np.pi/s)**2;val=(beta*beta*rbd[i]*rbd[j]-mu*rbv[i]*rbv[j])*I0[i,j]
            else:
                mu=0. if other=='Cons' else .25+(q*np.pi/s)**2
                h=-mu*I0[i,j]+np.sin(float(t0))*boundary[i,0]*boundary[j,1]
                val=beta*beta*rbd[i]*rbd[j]*I0[i,j]+rbv[i]*rbv[j]*h
            analytic_surface[i,j]=val
    check('All volume energy blocks by quadrature',direct_volume,analytic_volume)
    check('All fixed-surface energy blocks by quadrature',direct_surface,analytic_surface)
    check('Surface Gram symmetry',analytic_surface,analytic_surface.T)
    eig=np.linalg.eigvalsh(analytic_surface)
    assert eig.min()>-1e-9*max(1.,eig.max())
    results['Surface Gram non-negative spectrum'].append(max(0.,-float(eig.min()))/max(1.,float(eig.max())))

# Check the omitted mode against every retained p in a representative finite
# triangle. The proof in the text uses exact SL/Fourier orthogonality for all p.
R1=plotted.radial(1,rr)[0]
for p in [0,2,3,4,5,6]:
    check('Triangular missing (1,2) radial counterexample',np.dot(wr,R1*plotted.radial(p,rr)[0]),0)

report={'status':'PASS','precision_decimal_digits':mp.mp.dps,
        'comparison_count':sum(map(len,results.values())),
        'checks':{k:{'count':len(v),'max_scaled_error':max(v)} for k,v in results.items()},
        'scope':'Finite representative numerical identities and matrix checks; see handbook audit report for analytical proof review and limits.'}
dest=Path(__file__).with_name('handbook_math_validation.json')
dest.write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
