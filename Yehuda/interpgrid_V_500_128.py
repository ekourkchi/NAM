#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 10-17-2018"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "v1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"

from bokeh.transform import linear_cmap
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


def deAdjust(cz):

    Cs = 299792.458
    zz = cz/Cs
    Omega_m = 0.27
    Omega_lambda = 0.73
    j0 = 1
    q0 = 0.5*(Omega_m-2.*Omega_lambda)
    fz = 1.+0.5*(1-q0)*zz-(1-q0-3*q0**2+j0)*(zz**2)/6.

    return cz/fz

##########################################

##########################################


def distAverage(X, D):
    return np.sum(X/D**2)/np.sum(1./D**2)
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

    datfile = '../secure/cf3_BGc_500_128/table3_linearDV.csv'
    test = np.genfromtxt(datfile, delimiter=',',
                         filling_values=None, names=True, dtype=None)

    pgc1 = test['PGC1']
    nest = test['Nest']
    dist = test['d']
    et = test['et']
    eb = test['eb']
    Vls = test['Vls']      # test['Vls']
    Vlsm = test['Vlsm']
    sgl_catalog = test['sgl']
    sgb_catalog = test['sgb']
    name1 = test['Name1']
    name2 = test['Name2']
    M12 = test['M12']

    for i in range(len(pgc1)):
        name1[i] = name1[i].strip()+'('+name2[i].strip()+')'

    return nest, pgc1, dist, et, eb, Vls, Vlsm, sgl_catalog, sgb_catalog, name1, M12

###############################################################


def interpgrid(l, b, cone=20, coordinates='supergalactic', veldist=-1, VD=10):

    if (coordinates == 'blank'):

        TOOLS = ['reset']
        p = figure(tools=TOOLS, toolbar_location="right",
                   plot_width=650, plot_height=450, title='Please Wait ...')
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
    h0 = 0.75

    Vx = np.load('../secure/cf3_BGc_500_128/VX_500_128.npy')
    Vy = np.load('../secure/cf3_BGc_500_128/VY_500_128.npy')
    Vz = np.load('../secure/cf3_BGc_500_128/VZ_500_128.npy')

    X = 500.*np.load('../secure/cf3_BGc_500_128/X_500_128.npy') / \
        h0/128.-500.*64/h0/128.
    Y = 500.*np.load('../secure/cf3_BGc_500_128/Y_500_128.npy') / \
        h0/128.-500.*64/h0/128.
    Z = 500.*np.load('../secure/cf3_BGc_500_128/Z_500_128.npy') / \
        h0/128.-500.*64/h0/128.

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

        delta = 500./128./h0
        ii = np.floor((x0+250/h0)/delta)
        jj = np.floor((y0+250/h0)/delta)
        kk = np.floor((z0+250/h0)/delta)

        inds = np.zeros(27).astype('int')
        gg = 0
        for pp in [ii-1, ii, ii+1]:
            for qq in [jj-1, jj, jj+1]:
                for rr in [kk-1, kk, kk+1]:
                    inds[gg] = pp*128*128+qq*128+rr
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

        if len(srt) < 8:
            ppp = len(srt)
        else:
            ppp = 8

        nearest = nearby[srt[0:8]]

        Dw = dis3[srt[0:ppp]]

        # averaging over the nearest points
        best_Vx = distAverage(Vx_[nearest[0:ppp]], Dw)
        best_Vy = distAverage(Vy_[nearest[0:ppp]], Dw)
        best_Vz = distAverage(Vz_[nearest[0:ppp]], Dw)
        cz0 = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + 100*h0*VD)
        cz0 = v3k2vls(l, b, cz0)

    for i in range(ndis):

        delta = 500./128./h0
        ii = np.floor((x[i]+250/h0)/delta)
        jj = np.floor((y[i]+250/h0)/delta)
        kk = np.floor((z[i]+250/h0)/delta)

        inds = np.zeros(27).astype('int')
        gg = 0
        for pp in [ii-1, ii, ii+1]:
            for qq in [jj-1, jj, jj+1]:
                for rr in [kk-1, kk, kk+1]:
                    inds[gg] = pp*128*128+qq*128+rr
                    gg += 1

        X_ = X[inds]
        Y_ = Y[inds]
        Z_ = Z[inds]
        Vx_ = Vx[inds]
        Vy_ = Vy[inds]
        Vz_ = Vz[inds]

        dis2 = (X_-x[i])**2 + (Y_-y[i])**2 + (Z_-z[i])**2

        nearby, = np.where(dis2 < 30)

        dis3 = dis2[nearby]

        # Sort distances
        srt = np.argsort(dis3)

        if len(srt) < 8:
            ppp = len(srt)
        else:
            ppp = 8

        nearest = nearby[srt[0:8]]

        Dw = dis3[srt[0:ppp]]

        # averaging over the nearest points
        best_Vx = distAverage(Vx_[nearest[0:ppp]], Dw)
        best_Vy = distAverage(Vy_[nearest[0:ppp]], Dw)
        best_Vz = distAverage(Vz_[nearest[0:ppp]], Dw)
        cz[i] = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + 100*h0*Dn[i])
        cz[i] = v3k2vls(l, b, cz[i])

    cz_d = deAdjust(cz)

    if VVV and VD >= np.min(cz_d) and VD <= np.max(cz_d):

        dist_list1 = []
        for jj in range(1, len(cz_d)):
            if (VD >= cz_d[jj-1] and VD < cz_d[jj]) or (VD <= cz_d[jj-1] and VD > cz_d[jj]):

                dist_list1.append(
                    lineaR(VD, cz_d[jj-1], Dn[jj-1], cz_d[jj], Dn[jj]))

        if len(dist_list1) == 0:
            VVV = False

    else:
        VVV = False

    nest, pgc1, dist_catalog, et, eb, Vls, Vlsm, sgl_catalog, sgb_catalog, name_catalog, M12 = read_dat()

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

            delta = 500./128./h0
            ii = np.floor((x+250/h0)/delta)
            jj = np.floor((y+250/h0)/delta)
            kk = np.floor((z+250/h0)/delta)

            inds = np.zeros(27).astype('int')
            gg = 0
            for pp in [ii-1, ii, ii+1]:
                for qq in [jj-1, jj, jj+1]:
                    for rr in [kk-1, kk, kk+1]:
                        inds[gg] = pp*128*128+qq*128+rr
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

            # Sort distances
            srt = np.argsort(dis3)

            if len(srt) < 8:
                ppp = len(srt)
            else:
                ppp = 8

            nearest = nearby[srt[0:8]]

            Dw = dis3[srt[0:ppp]]

            # averaging over the nearest points
            best_Vx = distAverage(Vx_[nearest[0:ppp]], Dw)
            best_Vy = distAverage(Vy_[nearest[0:ppp]], Dw)
            best_Vz = distAverage(Vz_[nearest[0:ppp]], Dw)
            Vcmb_model = (best_Vx*xhat + best_Vy*yhat +
                          best_Vz*zhat + 100*h0*Dist)

            with HiddenPrints():
                tr = wcs.Transformation(
                    "supergalactic j2000 j2000", "galactic")
            l, b = tr((sgl, sgb))

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

    p = figure(tools=TOOLS, toolbar_location="right",
               plot_width=650, plot_height=450)
    p.title.text_font_size = '13pt'

    p.grid.grid_line_color = "gainsboro"

    Vls_plot = Vls+0.

    if len(ind) > 0:

        M12_ind = M12[ind]
        indx = np.argsort(M12_ind)
        ind = ind[indx]

        dist_ind = dist_catalog[ind]
        indx, = np.where(dist_ind < 250)
        ind = ind[indx]
        ind = ind[::-1]

        M12_ind = M12[ind]
        indx = np.where(M12_ind < 10)
        ind_ = ind[indx]

        palette = ['black', 'chocolate', 'red', 'green']

        n_ind = len(ind)
        dye = np.zeros(n_ind)
        size = np.zeros(n_ind)
        dist_ind = dist_catalog[ind]
        Vls_plot_ind = Vls_plot[ind]
        et_ind = et[ind]
        eb_ind = eb[ind]
        M12_ind = M12[ind]

        # create the coordinates for the errorbars
        err_xs = []
        err_ys = []

        for iind in range(n_ind):
            if M12_ind[iind] < 10:
                size[iind] = 1
                dye[iind] = 3
            elif M12_ind[iind] < 100:
                size[iind] = 9
                dye[iind] = 2
            else:
                size[iind] = 20
                dye[iind] = 1

            err_xs.append((dist_ind[iind] - eb_ind[iind],
                           dist_ind[iind] + et_ind[iind]))
            err_ys.append((Vls_plot_ind[iind], Vls_plot_ind[iind]))

        # errorbars
        if len(err_xs) > 0:
            p.multi_line(err_xs, err_ys, color='gray', alpha=0.5)

        mapper = linear_cmap(field_name='DYE', palette=palette, low=0, high=3)

        source = ColumnDataSource({'pgc1': pgc1[ind],
                                   'nest': nest[ind],
                                   'D_input': dist_catalog[ind],
                                   'Vlsm': Vlsm[ind],
                                   'Vls_plot': Vls_plot_ind,
                                   'Vls': Vls[ind],
                                   'name': name_catalog[ind],
                                   'M12': M12_ind,
                                   'size': size,
                                   'DYE': dye
                                   })

        data_p = p.cross('D_input', 'Vls_plot', source=source,
                         size=7,
                         color="black",
                         alpha=0.4)

        _data = p.circle('D_input', 'Vls_plot', size="size", source=source,
                         line_color='black', fill_color=mapper, line_width=1)

        _data_p = p.cross(dist_catalog[ind_], Vls[ind_],
                          size=7,
                          color="green",
                          alpha=0.4,
                          legend='Data')

        ttp2 = """
    <div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">PGC1:</span>
            <span style="font-size: 14px; font-weight: bold;">@pgc1</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">Nest:</span>
            <span style="font-size: 14px; font-weight: bold;">@nest</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">Dist:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{0.0}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">V<sup>c</sup><sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vlsm</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">V<sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">M<sub>12</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@M12{int}</span>
        </div>          
    </div>     
     """

        hover = HoverTool(tooltips=ttp2, renderers=[data_p])

        hover.point_policy = 'snap_to_data'
        hover.line_policy = 'nearest'
        p.add_tools(hover)

    source = ColumnDataSource({'D_input': Dn, 'Vls': cz_d})
    curve = p.line('D_input', 'Vls', source=source,
                   line_width=2, color="blue", legend='Model')

    ttp = """
    <div>
        <div>
            <span style="font-size: 14px; color: blue;">D_input:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{int}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: blue;">V<sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls{int}</span>
        </div>  
    </div>     
     """
    hover = HoverTool(tooltips=ttp, renderers=[curve])

    hover.point_policy = 'snap_to_data'
    hover.line_policy = 'nearest'
    # hover.mode='vline'
    p.add_tools(hover)

    if VVV:
        vel_list1 = []
        p.line([0, np.max(dist_list1)], [VD, VD], line_width=2,
               color="maroon", line_dash='dotted')
        for dist in dist_list1:
            p.line([dist, dist], [-200, VD], line_width=2,
                   color="maroon", line_dash='dotted')
            vel_list1.append(VD)
    if DDD:
        dist_list1 = [VD]
        vel_list1 = [cz0]
        p.line([VD, VD], [-200, cz0], line_width=2,
               color="maroon", line_dash='dotted')
        p.line([0, VD], [cz0, cz0], line_width=2,
               color="maroon", line_dash='dotted')

    p.xaxis.axis_label = 'Distance [Mpc]'
    p.yaxis.axis_label = 'V_ls [km/s]'

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

    VgivesD = np.genfromtxt('../secure/cf3_BGc_500_128/VgivesD.csv',
                            delimiter=',', filling_values=None, names=True, dtype=None)
    Hub_d = VgivesD['d']
    Hub_Vls = VgivesD['Vls']
    Hub_Vlsm = VgivesD['Vlsm']

    # straigt Hubble line
    p.line(Hub_d, Hub_Vlsm, line_width=1,
           color="black", line_dash='dotted', alpha=0.8)

    # curved line
    p.line(Hub_d, Hub_Vls, line_width=1.5, color="red",
           line_dash='dashed', alpha=0.8, legend='H0=75')


# 3
    js_move = '''
        if(cb_obj.x >= fig.x_range.start && cb_obj.x <= fig.x_range.end &&
           cb_obj.y >= fig.y_range.start && cb_obj.y <= fig.y_range.end)
        {
            cross.spans.height.computed_location = cb_obj.sx;
            cross.spans.width.computed_location = cb_obj.sy;
        }
        else
        {
            cross.spans.height.computed_location = null;
            cross.spans.width.computed_location = null;
        }
     '''
    js_leave = '''
            cross.spans.height.computed_location = null;
            cross.spans.width.computed_location = null;
     '''
# 3


# 3
    if DDD:
        x0 = xhat*VD
        y0 = yhat*VD
        z0 = zhat*VD
        dis2 = (X-x0)**2 + (Y-y0)**2 + (Z-z0)**2

        # Reduce number of distances to check
        nearby, = np.where(dis2 < 120.)

        dis3 = dis2[nearby]

        # Sort distances
        srt = np.argsort(dis3)

        if len(srt) < 8:
            ppp = len(srt)
        else:
            ppp = 8

        nearest = nearby[srt[0:8]]

        Dw = dis3[srt[0:ppp]]

        # averaging over the nearest points
        best_Vx = distAverage(Vx[nearest[0:ppp]], Dw)
        best_Vy = distAverage(Vy[nearest[0:ppp]], Dw)
        best_Vz = distAverage(Vz[nearest[0:ppp]], Dw)
        cz0 = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + 100*h0*VD)
        cz0 = v3k2vls(l, b, cz0)

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
# 3
# 3

    p2 = figure(tools=TOOLS, toolbar_location="right",
                plot_width=650, plot_height=450)

    p2.grid.grid_line_color = "gainsboro"

    if len(ind) > 0:
        M12_ind = M12[ind]
        indx = np.argsort(M12_ind)
        ind = ind[indx]

        dist_ind = dist_catalog[ind]
        indx, = np.where(dist_ind < 250)
        ind = ind[indx]
        ind = ind[::-1]

        M12_ind = M12[ind]
        indx = np.where(M12_ind < 10)
        ind_ = ind[indx]

        palette = ['black', 'chocolate', 'red', 'green']

        n_ind = len(ind)
        dye = np.zeros(n_ind)
        size = np.zeros(n_ind)
        dist_ind = dist_catalog[ind]
        Vls_plot_ind = Vls_plot[ind]
        et_ind = et[ind]
        eb_ind = eb[ind]
        M12_ind = M12[ind]

        # create the coordinates for the errorbars
        err_xs = []
        err_ys = []

        for iind in range(n_ind):
            if M12_ind[iind] < 10:
                size[iind] = 1
                dye[iind] = 3
            elif M12_ind[iind] < 100:
                size[iind] = 9
                dye[iind] = 2
            else:
                size[iind] = 20
                dye[iind] = 1

            err_xs.append((dist_ind[iind] - eb_ind[iind],
                           dist_ind[iind] + et_ind[iind]))
            err_ys.append((Vls_plot_ind[iind], Vls_plot_ind[iind]))

        # errorbars
        # if len(err_xs)>0:
            #p2.multi_line(err_xs, err_ys, color='gray', alpha=0.5)

        mapper = linear_cmap(field_name='DYE', palette=palette, low=0, high=3)

        source = ColumnDataSource({'pgc1': pgc1[ind],
                                   'nest': nest[ind],
                                   'D_input': dist_catalog[ind],
                                   'Vlsm': Vlsm[ind],
                                   'Vls_plot': Vls_plot_ind,
                                   'Vls': Vls[ind],
                                   'name': name_catalog[ind],
                                   'M12': M12_ind,
                                   'size': size,
                                   'DYE': dye
                                   })

        data_p = p2.cross('D_input', 'Vlsm', source=source,
                          size=7,
                          color="black",
                          alpha=0.4)

        _data = p2.circle('D_input', 'Vlsm', size="size", source=source,
                          line_color='black', fill_color=mapper, line_width=1)

        _data_p = p2.cross(dist_catalog[ind_], Vlsm[ind_],
                           size=7,
                           color="black",
                           alpha=0.4,
                           legend='Data')

        ttp2 = """
    <div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">PGC1:</span>
            <span style="font-size: 14px; font-weight: bold;">@pgc1</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">Nest:</span>
            <span style="font-size: 14px; font-weight: bold;">@nest</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">Dist:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{0.0}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">V<sup>c</sup><sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vlsm</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">V<sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vls</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red; font-weight: bold;">M<sub>12</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@M12{int}</span>
        </div>          
    </div>     
      """
        hover = HoverTool(tooltips=ttp2, renderers=[data_p])

        hover.point_policy = 'snap_to_data'
        hover.line_policy = 'nearest'
        p2.add_tools(hover)

    source = ColumnDataSource({'D_input': Dn, 'Vlsm': cz})
    curve = p2.line('D_input', 'Vlsm', source=source,
                    line_width=2, color="blue", legend='Model')

    ttp = """
    <div>
        <div>
            <span style="font-size: 14px; color: blue;">D_input:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{int}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: blue;">V<sup>c</sup><sub>ls</sub>:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vlsm{int}</span>
        </div>  
    </div>     
     """
    hover = HoverTool(tooltips=ttp, renderers=[curve])

    hover.point_policy = 'snap_to_data'
    hover.line_policy = 'nearest'
    # hover.mode='vline'
    p2.add_tools(hover)

    if VVV:
        vel_list = []
        p2.line([0, np.max(dist_list)], [VD, VD], line_width=2,
                color="maroon", line_dash='dotted')
        for dist in dist_list:
            p2.line([dist, dist], [-200, VD], line_width=2,
                    color="maroon", line_dash='dotted')
            vel_list.append(VD)
    if DDD:
        dist_list = [VD]
        vel_list = [cz0]
        p2.line([VD, VD], [-200, cz0], line_width=2,
                color="maroon", line_dash='dotted')
        p2.line([0, VD], [cz0, cz0], line_width=2,
                color="maroon", line_dash='dotted')

    p2.xaxis.axis_label = 'Distance [Mpc]'
    p2.yaxis.axis_label = 'V^c_ls [km/s]'

    p2.xaxis.axis_label_text_font_size = "14pt"
    p2.yaxis.axis_label_text_font_size = "14pt"
    p2.yaxis.major_label_text_font_size = "12pt"
    p2.xaxis.major_label_text_font_size = "12pt"

    p2.x_range = p.x_range  # Range1d(0, 250)
    p2.y_range = p.y_range  # Range1d(-200, 15000)

    p2.legend.location = "top_left"
    p2.legend.label_text_font_size = "12pt"
    p2.legend.label_text_font = "times"
    p2.legend.label_text_color = "black"

    # Setting the second y axis range name and range
    p2.extra_y_ranges = {"foo": p2.y_range}
    p2.extra_x_ranges = {"joo": p2.x_range}

    # Adding the second axis to the plot.
    p2.add_layout(LinearAxis(y_range_name="foo"), 'right')
    p2.add_layout(LinearAxis(x_range_name="joo"), 'above')

    cross2 = CrosshairTool()
    # cross.dimensions='height'
    cross2.line_alpha = 0.3
    cross2.line_color = 'green'
    p2.add_tools(cross2)

    p2.line(np.arange(250), np.arange(250)*0.75*100., line_width=1.5,
            color="red", line_dash='dashed', legend='H0=75')
# 3
    args = {'cross': cross2, 'fig': p}
    p.js_on_event('mousemove', CustomJS(args=args, code=js_move))
    p.js_on_event('mouseleave', CustomJS(args=args, code=js_leave))
    args = {'cross': cross, 'fig': p2}
    p2.js_on_event('mousemove', CustomJS(args=args, code=js_move))
    p2.js_on_event('mouseleave', CustomJS(args=args, code=js_leave))
# 3

    lbl = Label(x=300, y=20, x_units='screen', y_units='screen',
                text='Observed', render_mode='css',
                border_line_color=None, border_line_alpha=1.0,
                background_fill_color='white', background_fill_alpha=1.0,
                text_font_size='14pt', text_color='red')
    p.add_layout(lbl)

    lbl = Label(x=200, y=20, x_units='screen', y_units='screen',
                text='Cosmologically Adjusted', render_mode='css',
                border_line_color=None, border_line_alpha=1.0,
                background_fill_color='white', background_fill_alpha=1.0,
                text_font_size='14pt', text_color='green')
    p2.add_layout(lbl)

    p = column(p, p2)


# Online Mode
    script, div = components(p)
    script = '\n'.join(['' + line for line in script.split('\n')])

    title = '<b>&nbsp;&nbsp;'+title_l+'</b>' + \
        str(L)+'&nbsp;&nbsp;&nbsp;<b>'+title_b+'</b>'+str(B)
    print('<tr><td>')
    print(title)
    print('</td></tr>')

    print('<tr><td>')
    print(div)
    print(script)
    print('</td></tr>')

    if DDD or VVV:

        print('<tr><td align="center"><table id="calc" border="1">')
        print('<tr><th id="calc"><b>D<sub>input</sub></b><br><i>Mpc</i></th><th id="calc"><span style="font-size: 14px; color: red; font-weight: bold;">V<sub>ls</sub> - Observed</span><br><i>km/s</i></th></tr>')
        for i in range(len(dist_list1)):
            print('<tr><td id="calc">'+'%.2f' % dist_list1[i]+'</td>')
            print('<td id="calc">'+'%d' % vel_list1[i]+'</td></tr>')
        print('</table></td></tr>')

        print('<tr><td></td></tr>')

        print('<tr><td align="center"><table id="calc" border="1">')
        print('<tr><th id="calc"><b>D<sub>input</sub></b><br><i>Mpc</i></th><th id="calc"><span style="font-size: 14px; color: green; font-weight: bold;">V<sup>c</sup><sub>ls</sub> - Adjusted</span><br><i>km/s</i></th></tr>')
        for i in range(len(dist_list)):
            print('<tr><td id="calc">'+'%.2f' % dist_list[i]+'</td>')
            print('<td id="calc">'+'%d' % vel_list[i]+'</td></tr>')
        print('</table></td></tr>')


# Online Mode

    # offline mode
    ##print 'Showing the Plot ... '
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
