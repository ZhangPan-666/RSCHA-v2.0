from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent / 'deps'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from matplotlib.transforms import Bbox

OUT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).parent / 'images'
OUT.mkdir(exist_ok=True)
GRID = True
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix',
    'font.size': 14, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'savefig.facecolor': 'white', 'axes.unicode_minus': False,
})
INK = '#263B4B'
MUTED = '#73838E'
TEAL = '#347F87'
BLUE = '#5B7EA7'
GOLD = '#B88942'
LIGHT = '#D5DFE5'

fig = plt.figure(figsize=(12.8, 12.8), facecolor='white')
left = fig.add_axes([0.035, 0.555, 0.443, 0.385])
right = fig.add_axes([0.5365, 0.555, 0.414, 0.385])
compression = fig.add_axes([.043, .055, .429, .385])
coupled = fig.add_axes([.531, .055, .429, .385])
for ax in [left, right, compression, coupled]:
    ax.set_axis_off()
fig.text(.035, .965, '(a)  Spherical-conical geometry', color=INK, fontsize=18.5, fontweight='bold')
fig.text(.522, .965, '(b)  Family-specific stable evaluation', color=INK, fontsize=18.5, fontweight='bold')
fig.text(.035, .472, '(c)  Residual compression and recovery', color=INK, fontsize=18.5, fontweight='bold')
fig.text(.522, .472, '(d)  Legendre roots and cap normalization', color=INK, fontsize=18.5, fontweight='bold')

# Orthographic projection: cap axis tilted relative to the geographic axis.
elev = np.deg2rad(19)
EX = np.array([1., 0., 0.])
EY = np.array([0., np.sin(elev), np.cos(elev)])
VIEW = np.array([0., -np.cos(elev), np.sin(elev)])
tilt = np.deg2rad(34)
AXIS = np.array([np.cos(tilt), 0., np.sin(tilt)])
U = np.array([-np.sin(tilt), 0., np.cos(tilt)])
W = np.array([0., 1., 0.])
# Clockwise azimuth rotation as seen from +z_g, keeping the geographic axis fixed.
azimuth = np.deg2rad(-12.)
RZ = np.array([[np.cos(azimuth),-np.sin(azimuth),0.],
               [np.sin(azimuth), np.cos(azimuth),0.], [0.,0.,1.]])
AXIS, U, W = (RZ @ vector for vector in [AXIS, U, W])
theta0 = np.deg2rad(32)
r1, r2, re = 1.00, 1.38, .94

def project(points):
    return np.asarray(points) @ np.column_stack([EX, EY])

def direction(theta, phi):
    theta, phi = np.broadcast_arrays(theta, phi)
    return (np.cos(theta)[..., None] * AXIS
            + np.sin(theta)[..., None] *
            (np.cos(phi)[..., None] * U + np.sin(phi)[..., None] * W))

def line3(points, **kw):
    p = project(points)
    return left.plot(p[:, 0], p[:, 1], **kw)

def label3(point, text, offset=(0, 0), **kw):
    x, y = project(point)
    kw.setdefault('zorder', 10)
    return left.text(x+offset[0], y+offset[1], text, color=kw.pop('color', INK), **kw)

def radial(theta, phi, low, high, **kw):
    d = direction(theta, phi)
    return line3(np.array([low*d, high*d]), **kw)

left.set_xlim(-1.02, 1.40)
left.set_ylim(-.98, 1.24)
left.set_aspect('equal', adjustable='box')

# Pale globe with deliberately sparse graticule.
left.add_patch(Circle((0, 0), re, fc='#F0F4F6', ec='#CCD7DE', lw=.9, zorder=0))
t = np.linspace(0, 2*np.pi, 401)
equator = re*np.array([np.cos(t), np.sin(t), np.zeros_like(t)]).T
near = equator @ VIEW >= 0
for mask, style, alpha in [(~near, (0, (3, 4)), .6), (near, '-', .85)]:
    p = project(equator)
    q = np.where(mask[:, None], p, np.nan)
    left.plot(q[:, 0], q[:, 1], color='#A7B7C1', lw=.85, ls=style, alpha=alpha, zorder=1)
for longitude in [-55, 45]:
    a = np.linspace(-np.pi/2, np.pi/2, 201)
    lon = np.deg2rad(longitude)
    pts = re*np.array([np.cos(a)*np.cos(lon), np.cos(a)*np.sin(lon), np.sin(a)]).T
    line3(pts, color='#D7E0E5', lw=.55, zorder=1)

# Smooth vector silhouettes avoid the distracting grid of a triangulated surface.
def hull(points):
    pts = sorted(set(map(tuple, points)))
    def cross(o,a,b):
        return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    low=[]
    for p in pts:
        while len(low)>=2 and cross(low[-2],low[-1],p)<=0: low.pop()
        low.append(p)
    high=[]
    for p in reversed(pts):
        while len(high)>=2 and cross(high[-2],high[-1],p)<=0: high.pop()
        high.append(p)
    return np.array(low[:-1]+high[:-1])

theta = np.linspace(0, theta0, 61)
phi = np.linspace(0, 2*np.pi, 241)
th, ph = np.meshgrid(theta, phi, indexing='ij')
cap = direction(th,ph).reshape(-1,3)
left.add_patch(Polygon(hull(project(r2*cap)),closed=True,fc='#CFDDED',ec='none',alpha=.68,zorder=3))
front_phi=np.linspace(np.pi,2*np.pi,301)
wall=np.vstack([project(r2*direction(theta0,front_phi)),project(r1*direction(theta0,front_phi[::-1]))])
left.add_patch(Polygon(wall,closed=True,fc='#E8CF9E',ec='none',alpha=.45,zorder=3.1))
left.add_patch(Polygon(hull(project(r1*cap)),closed=True,fc='#AED0CD',ec='none',alpha=.60,zorder=3.2))

# A sparse white graticule on S2 only.
# This preserves the distinct outline of S1 and keeps all labels above the grid.
if GRID:
    def white_grid(points):
        front = np.asarray(points) @ VIEW >= 0
        xy = project(points)
        xy = np.where(front[:,None],xy,np.nan)
        left.plot(xy[:,0],xy[:,1],color='white',alpha=.80,lw=1.20,
                  ls=(0,(2.4,2.5)),zorder=3.8,dash_capstyle='round')
    for polar_fraction in [.36,.72]:
        white_grid(r2*direction(theta0*polar_fraction,np.linspace(0,2*np.pi,801)))
    for azimuth_deg in [-135,-75,-15,45]:
        white_grid(r2*direction(np.linspace(0,theta0,401),np.deg2rad(azimuth_deg)))

# Thin surface rims; hidden halves remain visible as unobtrusive dashed lines.
phi_dense = np.linspace(0, 2*np.pi, 801)
for radius, color in [(r1, TEAL), (r2, BLUE)]:
    rim = radius*direction(theta0, phi_dense)
    d = rim @ VIEW
    for mask, style, alpha, width in [(d < 0, (0, (3, 3)), .46, .9), (d >= 0, '-', .95, 1.45)]:
        xy = project(rim)
        xy = np.where(mask[:,None], xy, np.nan)
        left.plot(xy[:,0], xy[:,1], color=color, ls=style, alpha=alpha, lw=width, zorder=4)
# Meridian boundaries of the cone.
for angle in [0., np.pi]:
    radial(theta0, angle, r1, r2, color=GOLD, lw=1.45, zorder=5)

# Geographic axis and tilted cap axis.
line3(np.array([[0.,0.,-1.0], [0.,0.,1.2]]), color='#A1B0BA', lw=.95, ls=(0,(5,4)), zorder=2)
label3(np.array([0.,0.,1.2]), '$z_g$', (.02, .02), fontsize=17)
cap_end = project(1.62*AXIS)
left.add_patch(FancyArrowPatch((-.09*project(AXIS)), cap_end, arrowstyle='-|>', mutation_scale=11,
                               color=INK, lw=1.05, ls=(0,(5,4)), zorder=5))
left.text(cap_end[0]-.03, cap_end[1]+.065, 'Cap axis', color=INK, fontsize=13, ha='right')

# Radii follow a common meridian and use the manuscript's r1/r2 notation.
radial(theta0, 0., 0, r2, color=MUTED, lw=.95, ls=(0,(4,4)), zorder=5)
radial(theta0, np.pi, 0, r1, color=MUTED, lw=.95, ls=(0,(4,4)), zorder=5)
label3(.63*r2*direction(theta0,0), '$r_2$', (-.105,.015), fontsize=18)
label3(.55*r1*direction(theta0,np.pi), '$r_1$', (.0,-.12), fontsize=18)
left.scatter([0],[0], s=14, color=INK, zorder=7)
left.text(-.075,-.13,'$O$', fontsize=17, color=INK)

# Cap half-angle, measured from the cap axis to the conical boundary.
arc_t = np.linspace(0,theta0,101)
line3(.45*direction(arc_t,0), color=GOLD, lw=1.5, zorder=7)
label3(.50*direction(theta0*.53,0), r'$\theta_0$', (.025,.015), fontsize=18, color=GOLD)

# A field point and its colatitude in cap coordinates.
ptheta, pphi = np.deg2rad(22), np.deg2rad(185)
p = 1.19*direction(ptheta,pphi)
radial(ptheta,pphi,0,1.19,color=INK,lw=.85,ls=(0,(2,3)),zorder=6)
left.scatter(*project(p),s=23,color=INK,zorder=8)
xy = project(p)
left.annotate(r'$P(r,\theta,\varphi)$', xy=xy, xytext=(.72,-.34),
              fontsize=15,color=INK,arrowprops=dict(arrowstyle='-',color=MUTED,lw=.75),zorder=8)
line3(.40*direction(np.linspace(0,ptheta,100),pphi),color=INK,lw=1.15,zorder=7)
left.text(.205,.024,r'$\theta$',fontsize=17,color=INK,zorder=10)

# Boundary labels are kept outside the colored surfaces.
label3(r1*direction(np.deg2rad(23),np.deg2rad(150)), '$S_1$', (-.10,-.02), fontsize=18, color=TEAL)
left.text(1.225,.64,'$S_2$',fontsize=18,color=BLUE,zorder=10)
omega = 1.18*direction(theta0,np.deg2rad(135))
label3(omega,r'$\Omega$',(.015,-.005),fontsize=22,color=GOLD)

# Geometric scope below the panel; angular symbols match the manuscript.
fig.text(.2565,.515, r'$\Omega:\quad r_1\leq r\leq r_2,\quad 0\leq\theta\leq\theta_0,\quad 0\leq\varphi<2\pi$',
         ha='center',color=INK,fontsize=14)

# Match the flowchart's visible footprint to the projected geometry.
for ax in [right,compression,coupled]:
    ax.set_xlim(0,1); ax.set_ylim(-.025,1.)

def box(ax,x,y,w,h,title,detail=None,color=INK,fill='#F7F9FA', title_size=15, detail_size=13.4):
    patch = FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.012',
                           fc=fill,ec=color,lw=1.05)
    ax.add_patch(patch)
    if detail is None:
        ax.text(x+w/2,y+h/2,title,ha='center',va='center',fontsize=title_size,color=color)
    else:
        title_y = .75 if '\n' in detail else .69
        detail_y = .32 if '\n' in detail else .28
        ax.text(x+w/2,y+h*title_y,title,ha='center',va='center',fontsize=title_size,
                   color=color,fontweight='bold')
        ax.text(x+w/2,y+h*detail_y,detail,ha='center',va='center',fontsize=detail_size,color=INK,
                   linespacing=1.35)
    return patch

def arrow(ax,start,end,color=MUTED):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=10,
                                    lw=1.,color=color,shrinkA=1.6,shrinkB=2.2))

# (b): family-specific evaluation, with explicitly different public scales.
box(right,.055,.855,.89,.125,'Validate input and select the regime',
    r'$m,\theta,\theta_0,\nu\ \mathrm{or}\ \tau$; regular axis / endpoint limits',
    title_size=15,detail_size=13.1)
right.plot([.5,.5],[.847,.803],color=MUTED,lw=1.)
right.plot([.165,.835],[.803,.803],color=MUTED,lw=1.)
for cx in [.165,.5,.835]: arrow(right,(cx,.803),(cx,.755))
families=[(.018,TEAL,'#F0F7F6','Legendre',r'$\lambda=\nu(\nu+1)$',
           'Scaled forms\n+ residual atlas\nAnalytic '+r'$\partial_\nu$'),
          (.353,BLUE,'#F1F5FA','Mehler',r'$\lambda=-(\tau^2+\frac{1}{4})$',
           'Positive series\nor residual atlas\nRiccati cap norm'),
          (.688,GOLD,'#FAF6EE','Constant',r'$\lambda=0$',
           'Analytic shape\nPositive integral\n2-D norm atlas')]
for x,c,f,title,spectrum,detail in families:
    box(right,x,.425,.294,.32,'',color=c,fill=f)
    right.text(x+.147,.697,title,ha='center',va='center',fontsize=15.5,fontweight='bold',color=c)
    right.text(x+.147,.640,spectrum,ha='center',va='center',fontsize=14,color=INK)
    right.plot([x+.045,x+.249],[.604,.604],color=c,lw=.6,alpha=.35)
    right.text(x+.147,.516,detail,ha='center',va='center',fontsize=13.3,color=INK,linespacing=1.40)
for cx in [.165,.5,.835]: right.plot([cx,cx],[.417,.376],color=MUTED,lw=1.)
right.plot([.165,.835],[.376,.376],color=MUTED,lw=1.)
arrow(right,(.5,.376),(.5,.333))
box(right,.055,.205,.89,.120,'Preserve scale and floating-point range',
    'Component-wise signs, log magnitudes and status',title_size=14.7,detail_size=13.1)
arrow(right,(.5,.197),(.5,.149))
box(right,.055,.010,.89,.131,'Values and angular derivatives',
    r'$L,L_\theta,L_{\theta\theta};\qquad\widehat K,\widehat C\ \mathrm{and\ derivatives}$',
    color=TEAL,fill='#F0F7F6',title_size=15,detail_size=14)
fig.text(.7435,.515,'Residual recovery: (c).  Legendre roots and normalization: (d).',
         ha='center',color=MUTED,fontsize=12.6)

# (c): offline compression is separate from the online admissibility decision.
compression.text(.005,.992,'OFFLINE',ha='left',va='bottom',fontsize=11.8,color=MUTED)
box(compression,.055,.835,.89,.120,'Generate carrier-removed source fields',
    'Direct kernels; Riccati and positive-series norms',title_size=14.6,detail_size=13.1)
arrow(compression,(.5,.827),(.5,.797))
box(compression,.055,.665,.89,.124,'Compress with Chebyshev polynomials',
    'Cell eligibility checks; checksummed atlas',title_size=14.6,detail_size=13.3)
compression.plot([.01,.99],[.612,.612],color=LIGHT,lw=.75)
compression.text(.005,.586,'ONLINE',ha='left',va='center',fontsize=11.8,color=MUTED)
arrow(compression,(.5,.657),(.5,.556))
diamond=np.array([[.5,.548],[.723,.455],[.5,.362],[.277,.455]])
compression.add_patch(Polygon(diamond,closed=True,fc='#F7F9FA',ec=INK,lw=1.05))
compression.text(.5,.455,'Admissible\nlookup?',ha='center',va='center',color=INK,fontsize=14.6,linespacing=1.15)
arrow(compression,(.277,.455),(.245,.323))
arrow(compression,(.723,.455),(.755,.323))
compression.text(.176,.383,'Yes',ha='center',color=TEAL,fontsize=13.1)
compression.text(.824,.383,'No',ha='center',color=GOLD,fontsize=13.1)
box(compression,.017,.165,.456,.150,'Nested Clenshaw',
    'Restore carrier\nand normalization',color=BLUE,fill='#F1F5FA',title_size=14.7,detail_size=13.5)
box(compression,.527,.165,.456,.150,'Non-table route',
    'Family-specific\nanalytic / direct path',color=GOLD,fill='#FAF6EE',title_size=14.7,detail_size=13.5)
for cx in [.245,.755]: compression.plot([cx,cx],[.157,.123],color=MUTED,lw=1.)
compression.plot([.245,.755],[.123,.123],color=MUTED,lw=1.)
arrow(compression,(.5,.123),(.5,.106))
box(compression,.055,-.008,.89,.106,'Consistent signed-log output',
    'Report failure if no admissible route exists',title_size=14.6,detail_size=13.1)
fig.text(.2565,.024,'Smooth residuals and log norms are tabulated.',ha='center',color=MUTED,fontsize=12.6)

# (d): preserve the coupling among degree sensitivity, roots and cap norms.
box(coupled,.055,.855,.89,.125,'Analytic degree sensitivity',
    r'$(L,L_\theta,S,S_\theta),\quad S=\partial_\nu L$',title_size=15.1,detail_size=14.3)
coupled.plot([.5,.5],[.847,.805],color=MUTED,lw=1.)
coupled.plot([.25,.75],[.805,.805],color=MUTED,lw=1.)
for cx in [.25,.75]: arrow(coupled,(cx,.805),(cx,.773))
box(coupled,.025,.350,.45,.415,'',color=TEAL,fill='#F0F7F6')
box(coupled,.525,.350,.45,.415,'',color=BLUE,fill='#F1F5FA')
coupled.text(.25,.724,'Boundary roots',ha='center',color=TEAL,fontsize=15.4,fontweight='bold')
coupled.text(.75,.724,'Cap squared norm',ha='center',color=BLUE,fontsize=15.4,fontweight='bold')
for x in [.25,.75]: coupled.plot([x-.17,x+.17],[.699,.699],color=LIGHT,lw=.8)
for y,text in [(.655,'Action-based predictor'),(.579,'Computed sign bracket'),(.503,'Safeguarded refinement')]:
    coupled.text(.25,y,text,ha='center',va='center',color=INK,fontsize=13.1)
    arrow(coupled,(.25,y-.025),(.25,y-.047))
coupled.text(.25,.413,'Public-evaluator residual\nand ordering checks',ha='center',va='center',color=INK,fontsize=12.9,linespacing=1.25)
coupled.text(.75,.645,r'$I=\dfrac{\sin\theta_0}{2\nu+1}$',ha='center',va='center',color=INK,fontsize=15.5)
coupled.text(.75,.550,r'$\times\,(L_\theta S-LS_\theta)_{\theta_0}$',ha='center',va='center',color=INK,fontsize=14.3)
coupled.text(.75,.463,'Integral alternative',ha='center',va='center',color=INK,fontsize=13.3)
coupled.text(.75,.422,'when required',ha='center',va='center',color=INK,fontsize=13.1)
coupled.text(.75,.375,'Positivity / scale checks',ha='center',va='center',color=INK,fontsize=12.8)
for cx in [.25,.75]: coupled.plot([cx,cx],[.342,.307],color=MUTED,lw=1.)
coupled.plot([.25,.75],[.307,.307],color=MUTED,lw=1.)
arrow(coupled,(.5,.307),(.5,.280))
box(coupled,.055,.155,.89,.117,'Cap-unit Legendre basis',
    r'$\widehat L=L/\sqrt{I},\qquad\widehat L_\theta,\ \widehat L_{\theta\theta}$',title_size=14.9,detail_size=14.5)
arrow(coupled,(.5,.147),(.5,.109))
box(coupled,.055,-.008,.89,.109,'R-SCHA field kernels',
    r'$V\quad\longrightarrow\quad\mathbf{B}\quad\longrightarrow\quad\mathbf{T}$',
    color=TEAL,fill='#F0F7F6',title_size=14.9,detail_size=14.2)
fig.text(.7435,.024,'Spectral sensitivities link values, roots and cap norms.',ha='center',color=MUTED,fontsize=12.6)

for lo,hi in [(.545,.94),(.055,.445)]:
    fig.add_artist(plt.Line2D([.5,.5],[lo,hi],transform=fig.transFigure,color='#E0E6EA',lw=.7))
fig.add_artist(plt.Line2D([.035,.965],[.493,.493],transform=fig.transFigure,color='#E0E6EA',lw=.7))

base = OUT / 'RSCHA_Geometry_Algorithm'
for ext in ['pdf','svg']:
    fig.savefig(base.with_suffix('.'+ext),bbox_inches=None,pad_inches=0,
                metadata={'Title':'R-SCHA geometry, stable evaluation, residual compression and spectral coupling'})
fig.savefig(base.with_suffix('.png'),dpi=240,bbox_inches=None,pad_inches=0)
fig.savefig(OUT/(base.name+'_geometry.png'),dpi=220,
            bbox_inches=Bbox.from_bounds(.025*12.8,.500*12.8,.460*12.8,.495*12.8),pad_inches=0)
plt.close(fig)
print(str(base))
