#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 10-17-2018"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "v1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"

from bokeh.layouts import column
from bokeh.models.widgets import RadioButtonGroup
from math import *
from bokeh.models import LinearAxis
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, Range1d, Label, TapTool, OpenURL, CustomJS, CrosshairTool
from bokeh.embed import components
from bokeh.plotting import *
from kapteyn import wcs
import numpy as np
import os
import sys
reload(sys)
sys.setdefaultencoding('utf8')

##########################################


def v3k2vls(el, b, v3k):

    alpha = pi / 180.
    cosb = cos(b*alpha)
    sinb = sin(b*alpha)
    cosl = cos(el*alpha)
    sinl = sin(el*alpha)

    Vh = v3k + 25.2*cosl*cosb+245.7*sinl*cosb-276.8*sinb

    # vgsrm=Vh+11.1*cosl*cosb+251.*sinl*cosb+7.25*sinb

    return Vh2Vls(el, b, Vh)
##########################################


def Vh2Vls(el, b, Vh):

    alpha = pi / 180.
    cosb = cos(b*alpha)
    sinb = sin(b*alpha)
    cosl = cos(el*alpha)
    sinl = sin(el*alpha)

    vls = float(Vh)-26.*cosl*cosb+317.*sinl*cosb-8.*sinb

    return vls
##########################################


def lineaR(x, x1, y1, x2, y2):
    y = (x-x1)*(y2-y1)/(x2-x1)+y1
    return y
##########################################


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
##########################################


def XYZ(sgl, sgb, dist):

    cl1 = cos(sgl*pi/180.)
    sl1 = sin(sgl*pi/180.)
    cb1 = cos(sgb*pi/180.)
    sb1 = sin(sgb*pi/180.)

    X = dist * cl1 * cb1
    Y = dist * sl1 * cb1
    Z = dist * sb1

    return X, Y, Z
##########################################


def loadGrid(npz_file):

    grid = np.load(npz_file)
    X = grid['X']
    Y = grid['Y']
    Z = grid['Z']
    Vx = grid['Vx']
    Vy = grid['Vy']
    Vz = grid['Vz']
    grid.close()

    return X, Y, Z, Vx, Vy, Vz


###############################################################
def read_dat():

    datfile = '../secure/cf3_model_yehuda/CF3_EDD.csv'
    test = np.genfromtxt(datfile, delimiter=',',
                         filling_values=None, names=True, dtype=None)

    C = 299792.458  # speed of light km/s
    rtod = 180./np.pi

    pgc = test['PGC']
    nest = test['Nest']
    pgc1 = test['PGC1']
    dist = test['Dist']
    ed = test['eDM']*0.5
    Vls_observe = test['Vls']      # test['Vls_observe']
    Vls_model = test['Vls']
    sgl_catalog = test['SGL']
    sgb_catalog = test['SGB']
    name = test['Name']

    n = len(pgc)

    ind, = np.where(dist < 250)  # only those galaxies within 250 Mpc

    return nest[ind], pgc1[ind], dist[ind], Vls_observe[ind], Vls_model[ind], sgl_catalog[ind], sgb_catalog[ind], name[ind]


###############################################################

def interpgrid(l, b, cone=20, coordinates='supergalactic', veldist=-1, VD=10):

    if (coordinates == 'blank'):

        TOOLS = ['reset']
        p = figure(tools=TOOLS, toolbar_location="right",
                   plot_width=650, plot_height=550, title='Please Wait ...')
        p.title.text_font_size = '13pt'
        p.title.text_color = "red"
        p.grid.grid_line_color = "gainsboro"
        p.line(np.arange(250), np.arange(250)*75, line_width=1,
               color="black", line_dash='dashed', legend='H0=75')
        p.xaxis.axis_label = 'D_input [Mpc]'
        p.yaxis.axis_label = 'Vls [km/s]'

        p.xaxis.axis_label_text_font_size = "14pt"
        p.yaxis.axis_label_text_font_size = "14pt"
        p.yaxis.major_label_text_font_size = "12pt"
        p.xaxis.major_label_text_font_size = "12pt"

        p.x_range = Range1d(0, 250)
        p.y_range = Range1d(-200, 15000)

        p.legend.location = "top_left"
        p.legend.label_text_font_size = "12pt"
        p.legend.label_text_font = "times"
        p.legend.label_text_color = "black"

        # Setting the second y axis range name and range
        p.extra_y_ranges = {"foo": p.y_range}
        p.extra_x_ranges = {"joo": p.x_range}

        # Adding the second axis to the plot.
        p.add_layout(LinearAxis(y_range_name="foo"), 'right')
        p.add_layout(LinearAxis(x_range_name="joo"), 'above')

        cross = CrosshairTool()
        # cross.dimensions='height'
        cross.line_alpha = 0.3
        cross.line_color = 'green'
        p.add_tools(cross)

        print('<tr><td>')
        script, div = components(p)
        script = '\n'.join(['' + line for line in script.split('\n')])
        print(div)
        print(script)
        print('</td></tr>')
        return l, b

    DDD = False
    VVV = False
    if veldist == 0 and VD <= 250 and VD >= 0:
        DDD = True
    if veldist == 1:
        VVV = True

    if coordinates == 'equatorial':
        with HiddenPrints():
            tr = wcs.Transformation("equatorial j2000 j2000", "supergalactic")
        sgl, sgb = tr((l, b))
    elif coordinates == 'galactic':
        with HiddenPrints():
            tr = wcs.Transformation("galactic j2000 j2000", "supergalactic")
        sgl, sgb = tr((l, b))
    else:
        sgl = l
        sgb = b

    # These are any kind of input parameters
    L = l
    B = b
    SGL = sgl
    SGB = sgb

    if coordinates == 'equatorial':
        with HiddenPrints():
            tr = wcs.Transformation("equatorial j2000 j2000", "galactic")
        l, b = tr((l, b))
    elif coordinates == 'supergalactic':
        with HiddenPrints():
            tr = wcs.Transformation("supergalactic j2000 j2000", "galactic")
        l, b = tr((sgl, sgb))

    #################################
    h0 = 0.76

    # x = (np.linspace(0, 1000, 512, endpoint=False)-500)#/h0
    # y = (np.linspace(0, 1000, 512, endpoint=False)-500)#/h0
    # z = (np.linspace(0, 1000, 512, endpoint=False)-500)#/h0
    #xyz = np.meshgrid(x,y,z)
    #X = xyz[0].flatten()
    #Y = xyz[1].flatten()
    #Z = xyz[2].flatten()
    Vx = np.load('../secure/cf3_model_yehuda/VX.npy')
    Vy = np.load('../secure/cf3_model_yehuda/VY.npy')
    Vz = np.load('../secure/cf3_model_yehuda/VZ.npy')

    #p = 0
    # for i in range(512):
    # for j in range(512):
    # for k in range(512):
    #X[p] = x[i]
    #Y[p] = y[j]
    #Z[p] = z[k]
    # p+=1

    #np.save('../secure/cf3_model_yehuda/X.npy', X)
    #np.save('../secure/cf3_model_yehuda/Y.npy', Y)
    #np.save('../secure/cf3_model_yehuda/Z.npy', Z)

    X = np.load('../secure/cf3_model_yehuda/X.npy')/h0
    Y = np.load('../secure/cf3_model_yehuda/Y.npy')/h0
    Z = np.load('../secure/cf3_model_yehuda/Z.npy')/h0

    #################################

    Dn = np.arange(1, 250, 2)
    dtor = np.pi/180.

    xhat = np.cos(sgl*dtor)*np.cos(sgb*dtor)
    yhat = np.sin(sgl*dtor)*np.cos(sgb*dtor)
    zhat = np.sin(sgb*dtor)

    x = xhat*Dn
    y = yhat*Dn
    z = zhat*Dn

    ndis = len(Dn)
    cz = np.zeros(ndis)

    if DDD:
        x0 = xhat*VD
        y0 = yhat*VD
        z0 = zhat*VD

        delta = 1000./512./0.76
        ii = np.floor((x0+500/0.76)/delta)
        jj = np.floor((y0+500/0.76)/delta)
        kk = np.floor((z0+500/0.76)/delta)

        inds = np.zeros(27).astype('int')
        gg = 0
        for pp in [ii-1, ii, ii+1]:
            for qq in [jj-1, jj, jj+1]:
                for rr in [kk-1, kk, kk+1]:
                    inds[gg] = pp*512*512+qq*512+rr
                    gg += 1

        X_ = X[inds]
        Y_ = Y[inds]
        Z_ = Z[inds]
        Vx_ = Vx[inds]
        Vy_ = Vy[inds]
        Vz_ = Vz[inds]

        dis2 = (X_-x0)**2 + (Y_-y0)**2 + (Z_-z0)**2

        # Reduce number of distances to check
        nearby, = np.where(dis2 < 30)

        dis3 = dis2[nearby]

        # Sort distances
        srt = np.argsort(dis3)

        nearest = nearby[srt[0:8]]

        best_X = X_[nearest[0]]
        best_Y = Y_[nearest[0]]
        best_Z = Z_[nearest[0]]
        best_Vx = Vx_[nearest[0]]
        best_Vy = Vy_[nearest[0]]
        best_Vz = Vz_[nearest[0]]
        cz0 = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + 100*h0*VD)
        cz0 = v3k2vls(l, b, cz0)

    for i in range(ndis):

        delta = 1000./512./0.76
        ii = np.floor((x[i]+500/0.76)/delta)
        jj = np.floor((y[i]+500/0.76)/delta)
        kk = np.floor((z[i]+500/0.76)/delta)

        inds = np.zeros(27).astype('int')
        gg = 0
        for pp in [ii-1, ii, ii+1]:
            for qq in [jj-1, jj, jj+1]:
                for rr in [kk-1, kk, kk+1]:
                    inds[gg] = pp*512*512+qq*512+rr
                    gg += 1

        X_ = X[inds]
        Y_ = Y[inds]
        Z_ = Z[inds]
        Vx_ = Vx[inds]
        Vy_ = Vy[inds]
        Vz_ = Vz[inds]

        dis2 = (X_-x[i])**2 + (Y_-y[i])**2 + (Z_-z[i])**2

        if True:  # try:
            # Reduce number of distances to check
            nearby, = np.where(dis2 < 30)

            dis3 = dis2[nearby]

            # Sort distances
            srt = np.argsort(dis3)

            nearest = nearby[srt[0:8]]

            best_X = X_[nearest[0]]
            best_Y = Y_[nearest[0]]
            best_Z = Z_[nearest[0]]
            best_Vx = Vx_[nearest[0]]
            best_Vy = Vy_[nearest[0]]
            best_Vz = Vz_[nearest[0]]

            cz[i] = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + 100*h0*Dn[i])
            cz[i] = v3k2vls(l, b, cz[i])
        # except:
            #print dis2
            # sys.exit()

    if VVV and VD >= np.min(cz) and VD <= np.max(cz):

        dist_list = []
        for jj in range(1, len(cz)):
            if (VD >= cz[jj-1] and VD < cz[jj]) or (VD <= cz[jj-1] and VD > cz[jj]):

                dist_list.append(
                    lineaR(VD, cz[jj-1], Dn[jj-1], cz[jj], Dn[jj]))

        if len(dist_list) == 0:
            VVV = False

    else:
        VVV = False

    nest, pgc1, dist_catalog, Vls_observe, Vls_model, sgl_catalog, sgb_catalog, name_catalog = read_dat()

    N = len(sgl_catalog)
    ind = np.zeros(N)

    ind[np.where(np.abs(sgl_catalog - sgl) < cone/np.cos(sgb*dtor))] += 1
    ind[np.where(np.abs(sgb_catalog - sgb) < cone)] += 1
    ind[np.where(dist_catalog > 5)] += 1

    ind, = np.where(ind == 3)

    for iii in ind:
        if dist_catalog[iii] < 300:

            sgl = sgl_catalog[iii]
            sgb = sgb_catalog[iii]
            xhat = np.cos(sgl*dtor)*np.cos(sgb*dtor)
            yhat = np.sin(sgl*dtor)*np.cos(sgb*dtor)
            zhat = np.sin(sgb*dtor)

            Dist = dist_catalog[iii]
            x = xhat*Dist
            y = yhat*Dist
            z = zhat*Dist

            delta = 1000./512./0.76
            ii = np.floor((x+500/0.76)/delta)
            jj = np.floor((y+500/0.76)/delta)
            kk = np.floor((z+500/0.76)/delta)

            inds = np.zeros(27).astype('int')
            gg = 0
            for pp in [ii-1, ii, ii+1]:
                for qq in [jj-1, jj, jj+1]:
                    for rr in [kk-1, kk, kk+1]:
                        inds[gg] = pp*512*512+qq*512+rr
                        gg += 1

            X_ = X[inds]
            Y_ = Y[inds]
            Z_ = Z[inds]
            Vx_ = Vx[inds]
            Vy_ = Vy[inds]
            Vz_ = Vz[inds]

            dis2 = (X_-x)**2 + (Y_-y)**2 + (Z_-z)**2
            nearby, = np.where(dis2 < 30)
            dis3 = dis2[nearby]
            srt = np.argsort(dis3)
            nearest = nearby[srt[0:8]]
            best_Vx = Vx_[nearest[0]]
            best_Vy = Vy_[nearest[0]]
            best_Vz = Vz_[nearest[0]]
            Vcmb_model = (best_Vx*xhat + best_Vy*yhat +
                          best_Vz*zhat + 100*h0*Dist)
            with HiddenPrints():
                tr = wcs.Transformation(
                    "supergalactic j2000 j2000", "galactic")
            l, b = tr((sgl, sgb))
            Vls_model[iii] = v3k2vls(l, b, Vcmb_model)
    # except:
        #print sgl, sgb
        #print Dist
        #print dis2
        # pass
    # PLOTS
    # if len(ind)>0:
        #plt.plot(d_calc[ind], cz_calc[ind], 'ro')
    #plt.plot(np.arange(40),np.arange(40)*h0*100., 'k--')
    #plt.plot(Dn, cz)
    # plt.show()

    # Bokeh PLOTS

    TOOLS = ['pan', 'tap', 'wheel_zoom', 'box_zoom', 'reset', 'save']

    if coordinates == "galactic":
        title_l = "Glon= "
        title_b = "Glat= "
    elif coordinates == "equatorial":
        title_l = "RA= "
        title_b = "Dec= "
    else:
        title_l = "SGL= "
        title_b = "SGB= "

    p = figure(tools=TOOLS, toolbar_location="right", plot_width=650,
               plot_height=550, title=title_l+str(L)+'   '+title_b+str(B))
    p.title.text_font_size = '13pt'

    p.grid.grid_line_color = "gainsboro"

    source = ColumnDataSource({'D_input': Dn, 'Vls_model': cz})
    curve = p.line('D_input', 'Vls_model', source=source,
                   line_width=2, color="blue", legend='Model')

    ttp = """
    <div>
        <div>
            <span style="font-size: 14px; color: blue;">D_input:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{int}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: blue;">Vls_model:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls_model{int}</span>
        </div>  
    </div>     
     """
    hover = HoverTool(tooltips=ttp, renderers=[curve])

    hover.point_policy = 'snap_to_data'
    hover.line_policy = 'nearest'
    # hover.mode='vline'
    p.add_tools(hover)

    Vls_plot = Vls_observe+0.

    if len(ind) > 0:
        source = ColumnDataSource({'pgc1': pgc1[ind], 'D_input': dist_catalog[ind], 'Vls_model': Vls_model[ind],
                                   'Vls_plot': Vls_plot[ind], 'Vls_observe': Vls_observe[ind], 'name': name_catalog[ind]})
        data_p = p.cross('D_input', 'Vls_plot', source=source, size=7, color="red",
                         alpha=0.8, hover_color="black", hover_alpha=1, legend='Data')

        ttp2 = """
    <div>
        <div>
            <span style="font-size: 14px; color: red;">PGC1:</span>
            <span style="font-size: 14px; font-weight: bold;">@pgc1</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red;">D_input:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{0.2f}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red;">Vls_model:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls_model</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red;">Vls_observe:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls_observe</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red;">Name:</span>
            <span style="font-size: 14px; font-weight: bold;">@name</span>
        </div>  
    </div>     
     """

        hover = HoverTool(tooltips=ttp2, renderers=[data_p])

        hover.point_policy = 'snap_to_data'
        hover.line_policy = 'nearest'
        p.add_tools(hover)

    p.line(np.arange(250), np.arange(250)*75., line_width=1,
           color="black", line_dash='dashed', legend='H0=75')

    if VVV:
        vel_list = []
        p.line([0, np.max(dist_list)], [VD, VD], line_width=2,
               color="maroon", line_dash='dotted')
        for dist in dist_list:
            p.line([dist, dist], [-200, VD], line_width=2,
                   color="maroon", line_dash='dotted')
            vel_list.append(VD)
    if DDD:
        dist_list = [VD]
        vel_list = [cz0]
        p.line([VD, VD], [-200, cz0], line_width=2,
               color="maroon", line_dash='dotted')
        p.line([0, VD], [cz0, cz0], line_width=2,
               color="maroon", line_dash='dotted')

    p.xaxis.axis_label = 'Distance [Mpc]'
    p.yaxis.axis_label = 'Vls [km/s]'

    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = "12pt"
    p.xaxis.major_label_text_font_size = "12pt"

    p.x_range = Range1d(0, 250)
    p.y_range = Range1d(-200, 15000)

    p.legend.location = "top_left"
    p.legend.label_text_font_size = "12pt"
    p.legend.label_text_font = "times"
    p.legend.label_text_color = "black"

    # Setting the second y axis range name and range
    p.extra_y_ranges = {"foo": p.y_range}
    p.extra_x_ranges = {"joo": p.x_range}

    # Adding the second axis to the plot.
    p.add_layout(LinearAxis(y_range_name="foo"), 'right')
    p.add_layout(LinearAxis(x_range_name="joo"), 'above')

    cross = CrosshairTool()
    # cross.dimensions='height'
    cross.line_alpha = 0.3
    cross.line_color = 'green'
    p.add_tools(cross)

    radio_button_group = RadioButtonGroup(
        labels=["Model", "Observed", "Hide"], active=1, callback=CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var Vls_obs = data['Vls_observe'];
    var Vls_mod = data['Vls_model'];
    var Vls_plot = data['Vls_plot'];
    var f = cb_obj.active;
    
    if (f==1) {
        for (var i = 0; i < Vls_plot.length; i++) {
            Vls_plot[i] = Vls_obs[i];
        }    
    }

    if (f==0) {
        for (var i = 0; i < Vls_plot.length; i++) {
            Vls_plot[i] = Vls_mod[i];
        }    
    }
    
    if (f==2) {
        for (var i = 0; i < Vls_plot.length; i++) {
            Vls_plot[i] = NaN;
        }    
    }

    source.change.emit();
"""))

    p = column(p, radio_button_group)
# Online Mode
    script, div = components(p)
    script = '\n'.join(['' + line for line in script.split('\n')])

    print('<tr><td>')
    print(div)
    print(script)
    print('</td></tr>')

    if DDD or VVV:

        print('<tr><td align="center"><table id="calc" border="1">')

        print('<tr><th id="calc"><b>D<sub>input</sub></b><br><i>Mpc</i></th><th id="calc"><b>V<sub>ls</sub> -model</b><br><i>km/s</i></th></tr>')
        for i in range(len(dist_list)):
            print('<tr><td id="calc">'+'%.2f' % dist_list[i]+'</td>')
            print('<td id="calc">'+'%d' % vel_list[i]+'</td></tr>')

        print('</table></td></tr>')


# Online Mode

    # offline mode
    #print 'Showing the Plot ... '
    # show(p)

    return SGL, SGB

###############################################################


def myCoords(sgl, sgb, addCord=True):

    if not addCord:
        s = '<tr><td align="center">'
        s += '<br>'
        s += '<table>'
        s += '<tr><td><p><font color="blue"><b>RA:</b></font></p></td><td></td><td><p><font color="blue"><b>&nbsp;&nbsp;Dec:</b></font></p></td><td></td></tr>'
        s += '<tr><td><p><font color="green"><b>Glon:</b></font></p></td><td></td><td><p><font color="green"><b>&nbsp;&nbsp;Glat:</b></font></p></td><td></td></tr>'
        s += '<tr><td><p><font color="red"><b>SGL:</b></font></p></td><td></td><td><p><font color="red"><b>&nbsp;&nbsp;SGB:</b></font></p></td><td></td></tr>'

        s += '</table>'

        s += '</td></tr>'

        print(s)
        return

    with HiddenPrints():
        tr = wcs.Transformation("supergalactic j2000 j2000", "equatorial")
    ra, dec = tr((sgl, sgb))
    with HiddenPrints():
        tr = wcs.Transformation("supergalactic j2000 j2000", "galactic")
    gl, gb = tr((sgl, sgb))

    s = '<tr><td align="center">'
    s += '<br>'
    s += '<table>'
    s += '<tr><td><p><font color="blue"><b>RA:</b></font></p></td><td>'+'%.5f' % ra + \
        '</td><td><p><font color="blue"><b>&nbsp;&nbsp;Dec:</b></font></p></td><td>' + \
        '%.5f' % dec+'</td></tr>'
    s += '<tr><td><p><font color="green"><b>Glon:</b></font></p></td><td>'+'%.5f' % gl + \
        '</td><td><p><font color="green"><b>&nbsp;&nbsp;Glat:</b></font></p></td><td>' + \
        '%.5f' % gb+'</td></tr>'
    s += '<tr><td><p><font color="red"><b>SGL:</b></font></p></td><td>'+'%.5f' % sgl + \
        '</td><td><p><font color="red"><b>&nbsp;&nbsp;SGB:</b></font></p></td><td>' + \
        '%.5f' % sgb+'</td></tr>'

    s += '</table>'

    s += '</td></tr>'

    print(s)


if __name__ == "__main__":

    #fname = 'grid_ascii.dat'
    # grid2npz(fname)

    print('<table>')

    try:
        sgl = float(sys.argv[1])
        sgb = float(sys.argv[2])
        cone = float(sys.argv[3])
        coordinates = str(sys.argv[4])

        while sgl < 0:
            sgl += 360.
        while sgl > 360:
            sgl -= 360.

        if sgb > 90 or sgb < -90:
            sgb = float('exception')

        if (coordinates != 'blank'):

            if len(sys.argv) > 5:
                veldist = float(sys.argv[5])
                VD = float(sys.argv[6])
                sgl, sgb = interpgrid(
                    sgl, sgb, cone=cone, coordinates=coordinates, veldist=veldist, VD=VD)
            else:
                sgl, sgb = interpgrid(
                    sgl, sgb, cone=cone, coordinates=coordinates)
            myCoords(sgl, sgb, addCord=True)
        else:
            interpgrid(sgl, sgb, cone=cone, coordinates=coordinates)
            myCoords(0, 0, addCord=False)

    except:
        pass
        print('<tr><td>')
        print('<p><font color="red">&nbsp;&nbsp;Something went wrong !&nbsp;&nbsp;</font></p>')
        print('<p><font color="green">&nbsp;&nbsp;Check the input parameters.&nbsp;&nbsp;</font></p>')
        print('</td></tr>')

    print('</table>')


###############################################################
