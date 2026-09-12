"""Figure: the elementary inverse (mirror) and composition (magnification)
recipes for Schmidt circles, illustrated at level n = 9 (D = -80).

Left:  inverse = mirror image (in the imaginary axis, or in the unit circle).
Right: composition (3,4) * (4,0) = (12,16): the composed circle magnifies
       by 4 onto a translate of the first factor and by 3 onto a translate
       of the second; the three centers are collinear with the origin.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

n = 9

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.6),
                               gridspec_kw={'width_ratios': [1, 1.35]})

# ---------------- left: inverse = mirror ----------------
ax1.set_aspect('equal')
w = dict(q=3, x=4)                       # class (3,2,7), center (4+9i)/6
cx, cy, r = w['x'] / 6, n / 6, 1 / 6
ax1.add_patch(Circle((cx, cy), r, fc=(0.2, 0.4, 0.8, 0.35), ec='C0', lw=1.6))
ax1.add_patch(Circle((-cx, cy), r, fc=(0.9, 0.5, 0.1, 0.35), ec='C1', lw=1.6))
ax1.annotate(r'$\omega$   [$f$]', (cx + r + 0.05, cy), fontsize=12)
ax1.annotate(r'$[f]^{-1}$', (-cx - r - 0.42, cy), fontsize=12)
# mirror line
ax1.plot([0, 0], [0, 2.1], color='0.3', lw=1.2, ls='--')
ax1.annotate('mirror', (0.03, 2.0), fontsize=9, color='0.3')
# unit circle mirror + inversion image (m=8): center (4+9i)/16, r=1/16
th = [3.141592653589793 * t / 300 for t in range(301)]
import math
ax1.plot([math.cos(t) for t in th], [math.sin(t) for t in th],
         color='0.55', lw=1.0, ls=':')
ax1.add_patch(Circle((4 / 16, 9 / 16), 1 / 16,
                     fc=(0.9, 0.5, 0.1, 0.35), ec='C1', lw=1.4))
ax1.annotate(r'$[f]^{-1}$', (4 / 16 + 0.1, 9 / 16 - 0.02), fontsize=10)
ax1.annotate('unit-circle mirror', (0.75, 0.78), fontsize=8, color='0.45',
             rotation=-38)
ax1.plot([-1.6, 1.6], [0, 0], color='0.6', lw=1.0)
ax1.set_xlim(-1.55, 1.55)
ax1.set_ylim(-0.06, 2.15)
ax1.set_title('inverse of a circle = its mirror image', fontsize=11)

# ---------------- right: composition by magnification ----------------
ax2.set_aspect('equal')
# omega1 = (3, 4), omega2 = (4, 0), omega3 = (12, 16)
c1 = (4 / 6, n / 6);   r1 = 1 / 6
c2 = (0 / 8, n / 8);   r2 = 1 / 8
c3 = (16 / 24, n / 24); r3 = 1 / 24
g1 = (16 / 6, n / 6)               # 4 * omega3 = omega1 + 2
g2 = (16 / 8, n / 8)               # 3 * omega3 = omega2 + 2

ax2.add_patch(Circle(c1, r1, fc=(0.2, 0.4, 0.8, 0.35), ec='C0', lw=1.6))
ax2.add_patch(Circle(c2, r2, fc=(0.1, 0.6, 0.3, 0.35), ec='C2', lw=1.6))
ax2.add_patch(Circle(c3, r3, fc=(0.7, 0.2, 0.5, 0.45), ec='C3', lw=1.6))
ax2.add_patch(Circle(g1, r1, fc='none', ec='C0', lw=1.4, ls='--'))
ax2.add_patch(Circle(g2, r2, fc='none', ec='C2', lw=1.4, ls='--'))

ax2.annotate(r'$\omega_1$  ($q=3$)', (c1[0] - 0.30, c1[1] + r1 + 0.06), fontsize=11)
ax2.annotate(r'$\omega_2$  ($q=4$)', (c2[0] - 0.14, c2[1] + r2 + 0.06), fontsize=11)
ax2.annotate(r'$\omega_1 * \omega_2$  ($q=12$)', (c3[0] + 0.08, c3[1] - 0.03),
             fontsize=11)
ax2.annotate(r'$4\,\omega_3 = \omega_1 + 2$', (g1[0] - 0.33, g1[1] + r1 + 0.06),
             fontsize=10, color='C0')
ax2.annotate(r'$3\,\omega_3 = \omega_2 + 2$', (g2[0] - 0.30, g2[1] + r2 + 0.06),
             fontsize=10, color='C2')

# the magnification ray from the origin through the three collinear centers
tmax = 1.12
ax2.plot([0, g1[0] * tmax], [0, g1[1] * tmax], color='0.4', lw=1.0, ls=':')
ax2.annotate(r'magnification ray', (g1[0] * 0.52, g1[1] * 0.52 - 0.13),
             fontsize=9, color='0.35', rotation=29)
ax2.annotate(r'$\times 3$', (g2[0] * 0.82, g2[1] * 0.82 - 0.11), fontsize=10,
             color='0.35')
ax2.annotate(r'$\times 4$', (g1[0] * 0.93, g1[1] * 0.93 - 0.12), fontsize=10,
             color='0.35')

ax2.plot([-0.3, 3.2], [0, 0], color='0.6', lw=1.0)
ax2.set_xlim(-0.28, 3.15)
ax2.set_ylim(-0.06, 1.95)
ax2.set_title(r'composition: $x_3 \equiv x_1\ (2q_1)$, '
              r'$x_3 \equiv x_2\ (2q_2)$, $q_3 = q_1 q_2$', fontsize=11)

fig.suptitle(f'The class group acting on Schmidt circles at level '
             r'$\alpha = 9$  (discriminant $-80$)', fontsize=12)
fig.tight_layout()
fig.savefig('figures/composition-n9.png', dpi=150, bbox_inches='tight')
print('saved figures/composition-n9.png')
