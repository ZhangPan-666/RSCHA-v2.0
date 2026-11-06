from pathlib import Path
import json, math
import numpy as np
from scipy.special import sph_harm_y

def rz(a):
    c,s=math.cos(a),math.sin(a)
    return np.array([[c,-s,0],[s,c,0],[0,0,1.]])
def ry(a):
    c,s=math.cos(a),math.sin(a)
    return np.array([[c,0,s],[0,1.,0],[-s,0,c]])

# Implement only the finite d-sum printed in the handbook. Compare it with
# direct coordinate rotation of independently evaluated SciPy harmonics.
def d_element(l,m,q,beta):
    pref=math.sqrt(math.factorial(l+m)*math.factorial(l-m)*math.factorial(l+q)*math.factorial(l-q))
    c,s=math.cos(beta/2),math.sin(beta/2)
    return pref*sum((-1)**(t-m+q)*c**(2*l+m-q-2*t)*s**(q-m+2*t)/
        (math.factorial(l+m-t)*math.factorial(t)*math.factorial(q-m+t)*math.factorial(l-q-t))
        for t in range(max(0,m-q),min(l+m,l-q)+1))

def matrix(l,angles):
    a,b,g=angles
    return np.array([[np.exp(-1j*m*g)*d_element(l,m,q,b)*np.exp(-1j*q*a)
        for q in range(-l,l+1)] for m in range(-l,l+1)])

def harmonics(l,v):
    theta=math.acos(float(np.clip(v[2],-1,1)))
    phi=math.atan2(v[1],v[0])
    return np.array([sph_harm_y(l,m,theta,phi) for m in range(-l,l+1)])

maxerr=0.;unitarity=0.;coefficient_error=0.;count=0
for angles in [(.37,.83,-.51),(.9,1.8,2.1),(-1.2,1e-8,.6)]:
    a,b,g=angles;rot=rz(a)@ry(b)@rz(g)
    for l in range(7):
        D=matrix(l,angles)
        unitarity=max(unitarity,float(np.max(np.abs(D@D.conj().T-np.eye(2*l+1)))))
        coef=np.array([complex(math.sin(m+.3),math.cos(m+.7)) for m in range(-l,l+1)])
        rotated_coef=D.T@coef
        for theta,phi in [(.17,.4),(.7,-.6),(1.2,2.4),(2.8,-2.)]:
            v=np.array([math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)])
            local=harmonics(l,v);global_y=harmonics(l,rot.T@v)
            maxerr=max(maxerr,float(np.max(np.abs(global_y-D@local))))
            coefficient_error=max(coefficient_error,float(abs(coef@global_y-rotated_coef@local)))
            count+=2*l+1
assert maxerr<2e-13 and coefficient_error<5e-13 and unitarity<2e-13

# Unequal radial-branch truncation: hemispherical Dirichlet eigenmode nu=3,
# ra=1, rb=2; target V=r^3 P_3(cos theta) is purely external.
nu=3.;ra=1.;rb=2.
A=np.array([[nu+1,-nu],[-(nu+1)*(ra/rb)**(nu+2),nu*(rb/ra)**(nu-1)]])
boundary=np.array([-3.,12.]);paired=np.linalg.solve(A,boundary)
unpaired=paired.copy();unpaired[1]=0
residual=A@unpaired-boundary
assert np.max(np.abs(A@paired-boundary))<1e-13
assert np.linalg.norm(residual)>1

# A smooth harmonic potential u=z has different normal derivatives on
# adjacent spherical/lateral faces. A scalar jump is not globally H^{1/2}
# on the boundary. The cross-edge seminorm contains the divergent integral
# int_eps^1 int_eps^1 (s+t)^(-2) ds dt.
edge_theta=math.pi/3
face_normal_values=[-math.cos(edge_theta),-math.sin(edge_theta)]
eps=[1e-3,1e-6,1e-9]
cross_edge_integrals=[math.log((1+e)**2/(4*e)) for e in eps]
assert face_normal_values[0]!=face_normal_values[1]
assert cross_edge_integrals[0]<cross_edge_integrals[1]<cross_edge_integrals[2]

result=dict(status='PASS',wigner_defining_identity_checks=count,
    max_wigner_identity_absolute_error=maxerr,max_wigner_unitarity_error=unitarity,
    max_coefficient_rotation_error=coefficient_error,
    paired_boundary_example=dict(A=A.tolist(),target=boundary.tolist(),paired_coefficients=paired.tolist(),
        unpaired_boundary_residual=residual.tolist()),
    global_normal_trace_counterexample=dict(potential='u=z',edge_theta_rad=edge_theta,
        adjacent_normal_derivatives=face_normal_values,cutoffs=eps,seminorm_lower_bound_integrals=cross_edge_integrals),
    scope='Finite-degree numerical consistency checks and explicit analytical counterexamples; not a proof of all function-space statements.')
out=Path(__file__).parent/'independent_checks.json'
out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
