"""Rebuild the author-selected cylinder schematic as a vector PDF.

One-shot, manuscript-specific artwork source (2026-11-17), intentionally not a
general research instrument. Existing build_artifacts.py:make_geometry creates
an unrelated two-panel diagram; the selected WPS PDF contains only a raster
image, so neither provides editable vector artwork for this correction.

Run with Python + reportlab, from any directory. Output is Figure 1 only.
All drawing lengths and the illustrative Snell pair are schematic, not inputs
to a thermal calculation. Physical +y is upward, +z enters the front face,
and +x slopes downward/right, matching the selected artwork. The central ray
has x=0 and y=z*tan(theta_t), with 0 < theta_t < theta. Visual verification:
render the PDF, check the ray reaches the rear face above its center, and
check no labels or heat-loss arrows obscure the ray or coordinate axes.
"""

from math import atan2, cos, degrees, hypot, pi, radians, sin, tan, asin
from pathlib import Path

from reportlab.lib.colors import Color, HexColor, black
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'manuscript/figures/fig_physical_model_schematic_final.pdf'
W, H = 750, 430
ORIGIN = (240., 205.)
DZ = (255., 69.)
EX = (45., -38.)
EY = (0., 95.)
Z_SCALE = 180.
THETA = radians(50.)  # illustrative geometry, no nominal input is changed
THETA_T = asin(sin(THETA) / 1.45)
GRAY = HexColor('#8b9295')
BLUE = HexColor('#75a8cf')
RED = HexColor('#b84943')


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(a, s):
    return a[0] * s, a[1] * s


def ellipse(center, t, factor=1.):
    return add(center, add(scale(EX, factor*cos(t)), scale(EY, factor*sin(t))))


def path(points, fill=None, stroke=black, width=.9, closed=False):
    p = C.beginPath()
    p.moveTo(*points[0])
    for pt in points[1:]:
        p.lineTo(*pt)
    if closed:
        p.close()
    if fill is not None:
        C.setFillColor(fill)
    C.setStrokeColor(stroke)
    C.setLineWidth(width)
    C.drawPath(p, stroke=1, fill=int(fill is not None))


def line(a, b, color=black, width=.9, dash=None):
    C.saveState()
    C.setStrokeColor(color)
    C.setLineWidth(width)
    if dash:
        C.setDash(*dash)
    C.line(*a, *b)
    C.restoreState()


def arrow(a, b, color=black, width=.9, head=6):
    line(a, b, color, width)
    d = (b[0]-a[0], b[1]-a[1])
    n = hypot(*d)
    u, v = (d[0]/n, d[1]/n), (-d[1]/n, d[0]/n)
    path([b, add(add(b, scale(u, -head)), scale(v, .35*head)),
          add(b, scale(u, -.72*head)),
          add(add(b, scale(u, -head)), scale(v, -.35*head))],
         fill=color, stroke=color, width=.3, closed=True)


def text(x, y, s, size=14, font='Times-Roman', color=black, angle=0, align='left'):
    C.saveState()
    C.translate(x, y)
    C.rotate(angle)
    C.setFillColor(color)
    C.setFont({'Times-Roman':'FigureTimes', 'Times-Italic':'FigureTimesItalic',
               'Times-Bold':'FigureTimesBold'}.get(font,font), size)
    {'left': C.drawString, 'center': C.drawCentredString,
     'right': C.drawRightString}[align](0, 0, s)
    C.restoreState()


def theta_label(x, y, subscript=False):
    text(x, y, '\u03b8', size=15, font='Times-Italic')
    if subscript:
        text(x+8, y-3, 't', size=10, font='Times-Italic')


def arc(center, radius, a0, a1, color=GRAY):
    path([add(center, (radius*cos(a0+(a1-a0)*i/45),
                       radius*sin(a0+(a1-a0)*i/45))) for i in range(46)],
         stroke=color, width=.7)


def main():
    global C
    # Embed fonts, including theta: no dependence on a PDF viewer's Symbol font.
    for name, filename in [('FigureTimes','times.ttf'),
                           ('FigureTimesItalic','timesi.ttf'),
                           ('FigureTimesBold','timesbd.ttf')]:
        pdfmetrics.registerFont(TTFont(name, str(Path('C:/Windows/Fonts')/filename)))
    C = canvas.Canvas(str(OUTPUT), pagesize=(W, H), pdfVersion=(1, 4),
                      pageCompression=1, invariant=1,
                      initialFontName='FigureTimes')
    C.setTitle('Finite cylinder: positive-y oblique heating and refracted ray')
    C.setAuthor('Oblique thermal OPD manuscript')
    rear = add(ORIGIN, DZ)
    # Tangency points of the projected cylindrical sidewall.
    t = atan2(-(EY[0]*DZ[1]-EY[1]*DZ[0]),
              -(EX[0]*DZ[1]-EX[1]*DZ[0]))
    tangents = [ellipse(ORIGIN, t), ellipse(ORIGIN, t+pi)]
    top, bottom = sorted(tangents, key=lambda p:p[1], reverse=True)
    path([top, add(top,DZ), add(bottom,DZ), bottom],
         fill=HexColor('#f4f5f6'), stroke=black, closed=True)
    # Rear boundary: dashed hidden half, solid visible half.
    rear_points = [ellipse(rear, 2*pi*i/240) for i in range(241)]
    path(rear_points, fill=HexColor('#e8ebeb'), stroke=GRAY, width=.65, closed=True)
    C.saveState()
    C.setDash(1.8, 1.5)
    path([ellipse(rear, a) for a in [pi/2+pi*i/120 for i in range(121)]],
         stroke=GRAY, width=.6)
    C.restoreState()
    path([ellipse(rear, a) for a in [-pi/2+pi*i/120 for i in range(121)]],
         stroke=black, width=.9)
    # Coated front face and elliptical Gaussian footprint; vector-only bands.
    path([ellipse(ORIGIN, 2*pi*i/240) for i in range(241)],
         fill=HexColor('#e5e7e7'), stroke=black, width=1., closed=True)
    for j in range(30, 0, -1):
        f = .78*j/30
        weight = (1-j/31)*.38
        tint = Color(.91-.14*weight, .77-.40*weight, .76-.40*weight)
        path([ellipse(ORIGIN, 2*pi*i/180, f) for i in range(181)],
             fill=tint, stroke=tint, width=0, closed=True)
    # Gray reference axes retain the original projection and front-face origin.
    normal_angle = atan2(DZ[1], DZ[0])
    line(add(ORIGIN, scale(DZ, -.49)), add(rear, scale(DZ,.13)), GRAY, .75, (2,2))
    arrow(ORIGIN, (ORIGIN[0],388), black, .8, 6)
    text(256,374,'y',16,'Times-Italic')
    arrow(ORIGIN, add(ORIGIN, scale(EX,2.0)), GRAY, .8, 6)
    text(327,145,'x',16,'Times-Italic')
    arrow(rear, (630,310), black, .85, 6)
    text(640,303,'z',16,'Times-Italic')
    # The incident and refracted directions are projections of (0,tan(theta),1).
    incident_slope = DZ[1]/DZ[0] + EY[1]/Z_SCALE*tan(THETA)
    ray_slope = DZ[1]/DZ[0] + EY[1]/Z_SCALE*tan(THETA_T)
    start = (75., ORIGIN[1] - 165*incident_slope)
    band = (0.,12.)
    path([add(start,band), add(ORIGIN,band),
          add(ORIGIN,scale(band,-1)), add(start,scale(band,-1))],
         fill=HexColor('#efd2cf'), stroke=RED, width=1.2, closed=True)
    arrow(start, ORIGIN, RED, 1.6, 8)
    angle_i = atan2(incident_slope,1.)
    text(85,35,'Collimated Gaussian Beam',13,'Times-Roman',black,
         degrees(angle_i))
    # Refracted ray lies at x=0 and exits the rear plane with positive y.
    exit_y = EY[1] * (DZ[0]/Z_SCALE) * tan(THETA_T)
    exit_point = add(rear,(0.,exit_y))
    arrow(ORIGIN,exit_point,RED,1.9,8)
    text(328,276,'Refracted reference ray',13,'Times-Roman',RED,
         degrees(atan2(ray_slope,1.)))
    text(209,274,'Absorbed Gaussian',10,'Times-Roman',RED,-76)
    text(222,274,'heat flux',10,'Times-Roman',RED,-76)
    arc(ORIGIN,58,normal_angle+pi,angle_i+pi)
    theta_label(169,167)
    arc(ORIGIN,70,normal_angle,atan2(ray_slope,1.))
    theta_label(316,226,True)
    # Loss arrows on front, sidewall, and rear surfaces.
    for f in (.17,.45,.72):
        a=add(top,scale(DZ,f)); arrow(a,add(a,(-22.,28.)),BLUE,.85,6)
        a=add(bottom,scale(DZ,f)); arrow(a,add(a,(22.,-28.)),BLUE,.85,6)
    for p in [(534,326),(545,272),(537,221)]:
        arrow(p,add(p,(35.,10.)),BLUE,.85,6)
    for p in [(210,295),(205,240)]:
        arrow(p,add(p,(-27.,7.)),BLUE,.85,6)
    text(349,58,'Convective and radiative heat loss',14,'Times-Roman',BLUE,15)
    # Radius and face labels preserve the original semantic labels.
    radius_end = ellipse(rear,-pi/2+.45)
    arrow(rear,radius_end,GRAY,.65,5)
    text(510,223,'R',15,'Times-Italic',GRAY)
    line((221,85),(233,111),black,.75)
    text(231,64,'Front Face',14,'Times-Bold',align='center')
    text(231,44,'(Coated)',14,'Times-Italic',align='center')
    text(231,23,'z = 0',13,'Times-Italic',align='center')
    line((580,259),(543,268),GRAY,.7)
    text(600,267,'Rear Face',14,'Times-Bold')
    text(600,247,'(Uncoated)',14,'Times-Italic')
    text(632,225,'z = L',13,'Times-Italic',align='center')
    C.showPage()
    C.save()
    print(OUTPUT)
    print('Schematic checks: +y exit displacement = %.3f drawing units; '
          '0 < theta_t < theta; internal ray reaches z=L inside rear face.' % exit_y)
    assert 0 < THETA_T < THETA
    assert 0 < exit_y < EY[1]


if __name__ == '__main__':
    main()
