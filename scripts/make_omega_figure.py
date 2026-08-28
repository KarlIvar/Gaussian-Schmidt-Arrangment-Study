"""Figure: the atoms of Omega = {X in SL(2,Z[i]) : X(H) subset H}.

Left:  the atom disks -- the Apollonian gasket in the strip 0 < Im z < 1 generated
       by the line Im z = 0 (= R-hat), the line Im z = 1 and the Ford circles of
       radius 1/2 at the integers.  Shaded: the Ford disks (alpha = 1, the atoms
       associate to z -> z+i).
Right: three successive zooms of the nested chain D_1 > D_2 > ... produced by
       scripts/apollonian_chain_demo.py for z = 31/101 + 7/53 i.  Inside every
       atom disk sits a congruent copy of the whole gasket -- that self-similarity
       is what makes the factorisation algorithm run.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from fractions import Fraction
from omega import gasket, apollonian_chain, QC, ID, _pull_back

GK = [D for D in gasket(700) if D.A != 0]

fig = plt.figure(figsize=(15, 4.9))
gs = fig.add_gridspec(1, 4, width_ratios=[1.75, 1, 1, 1], wspace=0.28)
ax1 = fig.add_subplot(gs[0, 0])

# ---------------- left: the atoms ----------------
ax1.set_aspect('equal')
for D in GK:
    c, r = D.centre(), float(D.radius())
    if not (-0.15 <= float(c.real) <= 2.15):
        continue
    ford = (D.B.im == 1)
    ax1.add_patch(Circle((float(c.real), float(c.imag)), r,
                         facecolor=('#cfe3f5' if ford else 'none'),
                         edgecolor=('#2b6cb0' if ford else '#b04a2b'), lw=0.5, zorder=2))
ax1.add_patch(Rectangle((-0.15, 1.0), 2.35, 0.32, facecolor='#dddddd',
                        edgecolor='none', zorder=1))
ax1.axhline(0, color='k', lw=1.6, zorder=3)
ax1.axhline(1, color='#444', lw=1.0, zorder=3)
ax1.text(2.08, 1.14, r'$\{\operatorname{Im}z>1\}=T_i(\mathbb{H})$', ha='right', fontsize=10)
ax1.text(2.08, -0.14, r'$\hat{\mathbb{R}}$', ha='right', fontsize=12)
ax1.set_xlim(-0.15, 2.15); ax1.set_ylim(-0.22, 1.34)
ax1.set_title(r'the atom disks of $\Omega$: an Apollonian gasket'
              '\n' r'(shaded: the Ford disks, $\alpha=1$)', fontsize=11)
ax1.set_xticks([0, 1, 2]); ax1.set_yticks([0, 0.5, 1])

# ---------------- right: successive zooms ----------------
z = QC(Fraction(31, 101), Fraction(7, 53))
chain, _ = apollonian_chain(z, 6, nmax=4000)
zx, zy = float(z.real), float(z.imag)
Y = ID
mats = []
for X, D in chain:
    Y = Y * X
    mats.append((Y, D))

for panel, k in enumerate([0, 1, 2]):
    ax = fig.add_subplot(gs[0, panel + 1])
    ax.set_aspect('equal')
    Yk, Dk = mats[k]
    Dnext = mats[k + 1][1]
    for Gd in GK:                                  # the gasket inside D_k is Y_k(gasket)
        im = _pull_back(Gd, Yk.inv())
        if im.A == 0:
            continue
        c, r = im.centre(), float(im.radius())
        if r < 3e-4 * float(Dk.radius() if Dk.A else 1):
            continue
        ax.add_patch(Circle((float(c.real), float(c.imag)), r, facecolor='none',
                            edgecolor='#c9c9c9', lw=0.5, zorder=1))
    for D, col, lab in ((Dk, '#2b6cb0', f'$D_{{{k+1}}}$'), (Dnext, '#c23b22', f'$D_{{{k+2}}}$')):
        if D.A == 0:
            ax.axhline(float(D.height()), color=col, lw=1.8, zorder=3, label=lab)
        else:
            c, r = D.centre(), float(D.radius())
            ax.add_patch(Circle((float(c.real), float(c.imag)), r, facecolor='none',
                                edgecolor=col, lw=1.8, zorder=3, label=lab))
    ax.plot([zx], [zy], marker='*', color='#111111', ms=7, zorder=6, lw=0)
    h = 1.35 * float(Dk.radius()) if Dk.A else 1.2
    cx = float(Dk.centre().real) if Dk.A else zx
    cy = float(Dk.centre().imag) if Dk.A else 1.0
    ax.set_xlim(cx - h, cx + h); ax.set_ylim(cy - h, cy + h)
    ax.set_title(f'curvatures {-Dk.A}, {-Dnext.A}', fontsize=10)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.set_xticks([round(cx - h/2, 6), round(cx + h/2, 6)])
    ax.set_yticks([round(cy - h/2, 6), round(cy + h/2, 6)])

fig.text(0.68, 0.015, r'$z=\frac{31}{101}+\frac{7}{53}i$ (star): successive zooms into the nested atom disks '
         r'$\mathbb{H}\supset D_1\supset D_2\supset D_3\supset D_4$; grey = the gasket inside $D_k$',
         ha='center', fontsize=10)
fig.subplots_adjust(left=0.045, right=0.99, top=0.84, bottom=0.13)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures', 'omega-atoms.png')
plt.savefig(out, dpi=140)
print("wrote", os.path.normpath(out))
