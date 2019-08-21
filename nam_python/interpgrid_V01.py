#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 10-17-2018"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "v1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"

import os, sys
reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np
from kapteyn import wcs
from bokeh.plotting import *
from bokeh.embed import components
from bokeh.models import ColumnDataSource, LabelSet, HoverTool, Range1d, Label, TapTool, OpenURL, CustomJS, CrosshairTool
from bokeh.models import LinearAxis
from bokeh.transform import linear_cmap
##########################################
def lineaR(x,x1,y1,x2,y2):
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
# Reading the input file, e.g. grid_ascii.dat
# 6 columns (x,y,z,Vx,Vy,Vz)
def grid2npz(fname, n_skip=0):
     
     name = fname.split('.dat')[0]
     dtype1 = np.dtype([('x', 'f8'), ('y', 'f8'), ('z', 'f8'), ('Vx', 'f8'), ('Vy', 'f8'), ('Vz', 'f8')])
     a = np.loadtxt(fname, dtype=dtype1, skiprows=n_skip)
     
     X = a['x']
     Y = a['y']
     Z = a['z']     
     Vx = a['Vx']
     Vy = a['Vy']
     Vz = a['Vz']   

     np.savez_compressed(name+'.npz', X=X,Y=Y,Z=Z,Vx=Vx,Vy=Vy,Vz=Vz)  
     
     return X,Y,Z,Vx,Vy,Vz
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
     
     return X,Y,Z,Vx,Vy,Vz



###############################################################
def read_dat(version='V6'):
    

    datfile = '../secure/nam_python/test'+version+'_merged.csv'
    dtype1 = np.dtype([('nest', 'i'), ('pgc1', 'i'), ('dist', 'f8'), ('ed', 'f8'), ('Vgsr_observe', 'i'), ('Vgsr_model', 'i'), ('sgl', 'f8'), ('sgb', 'f8'), ('name', 'U30')])
    test = np.genfromtxt(datfile, dtype=dtype1, delimiter='\t', usecols=(0,1,2,4,6,7, 15, 16, 22), filling_values=None)      
    
    nest = test['nest']
    pgc1 = test['pgc1']
    dist = test['dist']
    ed   = test['ed']
    Vgsr_observe = test['Vgsr_observe']
    Vgsr_model   = test['Vgsr_model']
    sgl_catalog  = test['sgl'] 
    sgb_catalog  = test['sgb']
    name = test['name']
    
    ind = np.where(dist>0)
    
  
    return nest[ind],pgc1[ind],dist[ind],Vgsr_observe[ind],Vgsr_model[ind],sgl_catalog[ind],sgb_catalog[ind], name[ind], ed[ind]
    

###############################################################

def interpgrid(l, b, cone=20, coordinates='supergalactic', veldist=-1, VD=10):
     
     if (coordinates=='blank'):
         
         TOOLS = ['reset']
         p = figure(tools=TOOLS, toolbar_location="below", plot_width=650, plot_height=550, title='Please Wait ...')
         p.title.text_font_size = '13pt'
         p.title.text_color = "red"
         p.grid.grid_line_color="gainsboro"
         p.line(np.arange(50),np.arange(50)*75, line_width=1, color="black", line_dash='dashed', legend='H0=75')
         p.xaxis.axis_label = 'D_input [Mpc]'
         p.yaxis.axis_label = 'Vgsr_model [km/s]'
        
         p.xaxis.axis_label_text_font_size = "14pt"
         p.yaxis.axis_label_text_font_size = "14pt"
         p.yaxis.major_label_text_font_size = "12pt"
         p.xaxis.major_label_text_font_size = "12pt"
        
         p.x_range = Range1d(0, 42)
         p.y_range = Range1d(-200, 3250)
        
         p.legend.location = "top_left"
         p.legend.label_text_font_size = "12pt"
         p.legend.label_text_font = "times"
         p.legend.label_text_color = "black"
   
         # Setting the second y axis range name and range
         p.extra_y_ranges = {"foo": p.y_range}
         p.extra_x_ranges = {"joo": p.x_range}
        
         #Adding the second axis to the plot.  
         p.add_layout(LinearAxis(y_range_name="foo"), 'right')
         p.add_layout(LinearAxis(x_range_name="joo"), 'above')
               
         cross = CrosshairTool()
         #cross.dimensions='height'
         cross.line_alpha = 0.3
         cross.line_color = 'green'
         p.add_tools(cross) 
         
         print '<tr><td>'
         script, div = components(p)
         script = '\n'.join(['' + line for line in script.split('\n')])
         print div
         print script 
         print '</td></tr>'
         return l, b
     
     DDD = False
     VVV = False
     if veldist==0 and VD<=38 and VD>=0: DDD=True
     if veldist==1: VVV=True
 
     if coordinates=='equatorial':
         with HiddenPrints():
            tr = wcs.Transformation("equatorial j2000 j2000", "supergalactic")     
         sgl, sgb = tr((l, b))
     elif coordinates=='galactic':
         with HiddenPrints():
            tr = wcs.Transformation("galactic j2000 j2000", "supergalactic")
         sgl, sgb = tr((l, b))
     else:
         sgl = l
         sgb = b 
 
     
     X,Y,Z,Vx,Vy,Vz = loadGrid('../secure/nam_python/grid_ascii.npz')
     
     Dn = np.arange(38)+1
     dtor = np.pi/180.
     h0 = 0.75
     vsun = 240.
     
     xhat= np.cos(sgl*dtor)*np.cos(sgb*dtor)
     yhat= np.sin(sgl*dtor)*np.cos(sgb*dtor)
     zhat= np.sin(sgb*dtor)
     
     x = xhat*Dn
     y = yhat*Dn
     z = zhat*Dn   
     
     ndis = len(Dn)
     cz = np.zeros(ndis)
     
     if DDD:
        x0 = xhat*VD
        y0 = yhat*VD
        z0 = zhat*VD          
        dis2 = (X-x0)**2 + (Y-y0)**2 + (Z-z0)**2
         
        # Reduce number of distances to check
        nearby, = np.where(dis2<10.)
         
        dis3 = dis2[nearby]
         
        # Sort distances
        srt = np.argsort(dis3)
         
        nearest = nearby[srt[0:8]]

        best_X = X[nearest[0]]
        best_Y = Y[nearest[0]]
        best_Z = Z[nearest[0]]
        best_Vx = Vx[nearest[0]]
        best_Vy = Vy[nearest[0]]
        best_Vz = Vz[nearest[0]]
        cz0 = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + h0*VD)*100.
         
     
     for i in range(ndis):
         
         dis2 = (X-x[i])**2 + (Y-y[i])**2 + (Z-z[i])**2
         
         # Reduce number of distances to check
         nearby, = np.where(dis2<10.)
         
         dis3 = dis2[nearby]
         
         # Sort distances
         srt = np.argsort(dis3)
         
         nearest = nearby[srt[0:8]]

         best_X = X[nearest[0]]
         best_Y = Y[nearest[0]]
         best_Z = Z[nearest[0]]
         best_Vx = Vx[nearest[0]]
         best_Vy = Vy[nearest[0]]
         best_Vz = Vz[nearest[0]]
         cz[i] = (best_Vx*xhat + best_Vy*yhat + best_Vz*zhat + h0*Dn[i])*100.
         
     if VVV and VD>=np.min(cz) and VD<=np.max(cz):
         
         dist_list = []
         for jj in range(1,len(cz)):
             if (VD>=cz[jj-1] and VD<cz[jj]) or (VD<=cz[jj-1] and VD>cz[jj]):
                 
                 dist_list.append(lineaR(VD,cz[jj-1],Dn[jj-1], cz[jj],Dn[jj]))
        
         if len(dist_list)==0: VVV=False
         
     else: VVV=False
     
     
     nest, pgc1, dist_catalog, Vgsr_observe, Vgsr_model, sgl_catalog, sgb_catalog, name_catalog, ed_catalog = read_dat(version='V6')
     

     N = len(sgl_catalog)    
     ind = np.zeros(N)
     
     ind[np.where(np.abs(sgl_catalog - sgl) < cone/np.cos(sgb*dtor))]+=1
     ind[np.where(np.abs(sgb_catalog - sgb) < cone)]+=1
     ind[np.where(dist_catalog>5)]+=1
     
     ind, = np.where(ind==3)
     
     
     ##################### PLOTS
     #if len(ind)>0:
         #plt.plot(d_calc[ind], cz_calc[ind], 'ro')
     #plt.plot(np.arange(40),np.arange(40)*h0*100., 'k--')
     #plt.plot(Dn, cz)
     #plt.show()
     
     
     ##################### Bokeh PLOTS

     
     TOOLS = ['pan', 'tap', 'wheel_zoom', 'box_zoom', 'reset', 'save']
     
     if coordinates=="galactic": 
        title_l = "Glon= "
        title_b = "Glat= "
     elif coordinates=="equatorial":
        title_l = "RA= "
        title_b = "Dec= "
     else:
        title_l = "SGL= "
        title_b = "SGB= "         
     
     p = figure(tools=TOOLS, toolbar_location="below", plot_width=650, plot_height=550, title=title_l+str(l)+'   '+title_b+str(b))
     p.title.text_font_size = '13pt'
     
     p.grid.grid_line_color="gainsboro"
     

     if len(ind)>0:
         
            palette = ['black', 'chocolate', 'red', 'green']

            n_ind = len(ind)
            dye  = np.zeros(n_ind)
            size = np.zeros(n_ind)
            dist_ind = dist_catalog[ind]
            ed_ind = ed_catalog[ind]
            pgc1_ind = pgc1[ind]
            error = dist_ind*ed_ind
            Vgsr_ind = Vgsr_model[ind]
            # create the coordinates for the errorbars
            err_xs = []
            err_ys = []

            
            for iind in range(n_ind):
                if ed_ind[iind]<=0.05:
                    size[iind]=11
                    dye[iind] = 2
                elif ed_ind[iind]<=0.15:
                    size[iind]=6
                    dye[iind] = 3
                else:
                    size[iind]= 0
                if pgc1_ind[iind]==41361 or pgc1_ind[iind]==12651:
                    size[iind]=17
                    dye[iind] = 1
                
                if ed_ind[iind]<=0.15:
                    err_xs.append((dist_ind[iind]- error[iind], dist_ind[iind]+ error[iind]))
                    err_ys.append((Vgsr_ind[iind] , Vgsr_ind[iind]))
            
            if len(err_xs)>0:
                p.multi_line(err_xs, err_ys, color='gray')    
         
            mapper = linear_cmap(field_name='DYE', palette=palette, low=0, high=3)
            
            source = ColumnDataSource({'pgc1': pgc1_ind,
                                       'D_input': dist_ind, 
                                       'Vgsr_model': Vgsr_ind, 'Vgsr_observe':Vgsr_observe[ind], 'name':name_catalog[ind],
                                       'ed': ed_ind,
                                       'error': error,
                                       'size': size, 
                                       'DYE': dye
                                       })
            
            data_p = p.cross('D_input', 'Vgsr_model', 
                             source=source, size=7, color="black",                             legend='Data')
            _data = p.circle('D_input', 'Vgsr_model', size="size", source=source,
                             line_color="black", fill_color=mapper, line_width=1)


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
            <span style="font-size: 14px; color: red;">Vgsr_model:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vgsr_model</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red;">Vgsr_observe:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vgsr_observe</span>
        </div>  
        <div>
            <span style="font-size: 14px; color: red;">Name:</span>
            <span style="font-size: 14px; font-weight: bold;">@name</span>
        </div>  
    </div>     
     """     
            
            hover = HoverTool(tooltips=ttp2 , renderers=[data_p])
            
            #hover = HoverTool(tooltips=[ 
               #("PGC1", "@pgc1"),
               #("D_input", "@D_input{0.2f}"),
               #("Vgsr_model", "@Vgsr_model"),    
               #("Vgsr_observe", "@Vgsr_observe"),
               #("Name", "@name"),
               #], renderers=[data_p])
            hover.point_policy='snap_to_data'
            hover.line_policy='nearest'
            p.add_tools(hover)
     

     source = ColumnDataSource({'D_input': Dn, 'Vgsr_model': cz})
     curve = p.line('D_input', 'Vgsr_model', source=source, line_width=2, color="blue", legend='Model')
     
     
     ttp = """
    <div>
        <div>
            <span style="font-size: 14px; color: blue;">D_input:</span>
            <span style="font-size: 14px; font-weight: bold;">@D_input{int}</span>
        </div>
        <div>
            <span style="font-size: 14px; color: blue;">Vgsr_model:</span>
            <span style="font-size: 14px; font-weight: bold;">@Vgsr_model{int}</span>
        </div>  
    </div>     
     """
     hover = HoverTool(tooltips=ttp, renderers=[curve])
     
     hover.point_policy='snap_to_data'
     hover.line_policy='nearest'
     #hover.mode='vline'
     p.add_tools(hover) 
 
 
     p.line(np.arange(50),np.arange(50)*h0*100., line_width=1, color="black", line_dash='dashed', legend='H0=75')
     
     if VVV:
         vel_list = []
         p.line([0,np.max(dist_list)],[VD, VD], line_width=2, color="maroon", line_dash='dotted')
         for dist in dist_list:
             p.line([dist,dist],[-200, VD], line_width=2, color="maroon", line_dash='dotted')
             vel_list.append(VD) 
     if DDD:
         dist_list = [VD]
         vel_list = [cz0]
         p.line([VD,VD],[-200, cz0], line_width=2, color="maroon", line_dash='dotted')
         p.line([0,VD],[cz0, cz0], line_width=2, color="maroon", line_dash='dotted')
                
                
     p.xaxis.axis_label = 'D_input [Mpc]'
     p.yaxis.axis_label = 'Vgsr_model [km/s]'
     
     p.xaxis.axis_label_text_font_size = "14pt"
     p.yaxis.axis_label_text_font_size = "14pt"
     p.yaxis.major_label_text_font_size = "12pt"
     p.xaxis.major_label_text_font_size = "12pt"
     
     p.x_range = Range1d(0, 42)
     p.y_range = Range1d(-200, 3250)
     
     p.legend.location = "top_left"
     p.legend.label_text_font_size = "12pt"
     p.legend.label_text_font = "times"
     p.legend.label_text_color = "black"
     
     # Setting the second y axis range name and range
     p.extra_y_ranges = {"foo": p.y_range}
     p.extra_x_ranges = {"joo": p.x_range}
     
     #Adding the second axis to the plot.  
     p.add_layout(LinearAxis(y_range_name="foo"), 'right')
     p.add_layout(LinearAxis(x_range_name="joo"), 'above')
     
     
     cross = CrosshairTool()
     #cross.dimensions='height'
     cross.line_alpha = 0.3
     cross.line_color = 'green'
     p.add_tools(cross) 

######################### Online Mode
     script, div = components(p)
     script = '\n'.join(['' + line for line in script.split('\n')])
     
     print '<tr><td>'
     print div
     print script 
     print '</td></tr>'
     
     if DDD or VVV:
         
         print '<tr><td align="center"><table id="calc" border="1">'
         
         print '<tr><th id="calc"><b>D<sub>input</sub></b><br><i>Mpc</i></th><th id="calc"><b>V<sub>gsr</sub> -model</b><br><i>km/s</i></th></tr>'
         for i in range(len(dist_list)):
             print '<tr><td id="calc">'+'%.2f'%dist_list[i]+'</td>'
             print '<td id="calc">'+'%d'%vel_list[i]+'</td></tr>'
             
         
         print '</table></td></tr>'         
         

######################### Online Mode
    
     ## offline mode
     #show(p)
     
     return sgl, sgb

###############################################################
def myCoords(sgl,sgb, addCord = True):
    
    if not addCord:
        s= '<tr><td align="center">'
        s+= '<br>'
        s+= '<table>'
        s+= '<tr><td><p><font color="blue"><b>RA:</b></font></p></td><td></td><td><p><font color="blue"><b>&nbsp;&nbsp;Dec:</b></font></p></td><td></td></tr>'
        s+= '<tr><td><p><font color="green"><b>Glon:</b></font></p></td><td></td><td><p><font color="green"><b>&nbsp;&nbsp;Glat:</b></font></p></td><td></td></tr>'
        s+= '<tr><td><p><font color="red"><b>SGL:</b></font></p></td><td></td><td><p><font color="red"><b>&nbsp;&nbsp;SGB:</b></font></p></td><td></td></tr>'

        s+= '</table>'
        
        s+= '</td></tr>'
        
        print s  
        return
        
    
    with HiddenPrints():
        tr = wcs.Transformation("supergalactic j2000 j2000", "equatorial")     
    ra, dec = tr((sgl, sgb))
    with HiddenPrints():
        tr = wcs.Transformation("supergalactic j2000 j2000", "galactic")     
    gl, gb = tr((sgl, sgb))   
    
    

    s= '<tr><td align="center">'
    s+= '<br>'
    s+= '<table>'
    s+= '<tr><td><p><font color="blue"><b>RA:</b></font></p></td><td>'+'%.5f'%ra+'</td><td><p><font color="blue"><b>&nbsp;&nbsp;Dec:</b></font></p></td><td>'+'%.5f'%dec+'</td></tr>'
    s+= '<tr><td><p><font color="green"><b>Glon:</b></font></p></td><td>'+'%.5f'%gl+'</td><td><p><font color="green"><b>&nbsp;&nbsp;Glat:</b></font></p></td><td>'+'%.5f'%gb+'</td></tr>'
    s+= '<tr><td><p><font color="red"><b>SGL:</b></font></p></td><td>'+'%.5f'%sgl+'</td><td><p><font color="red"><b>&nbsp;&nbsp;SGB:</b></font></p></td><td>'+'%.5f'%sgb+'</td></tr>'

    s+= '</table>'
    
    s+= '</td></tr>'
    
    print s 


if __name__ == "__main__":

    #fname = 'grid_ascii.dat'
    #grid2npz(fname)    
    
    print '<table>'
    
    try:
        sgl = float(sys.argv[1])
        sgb = float(sys.argv[2])
        cone = float(sys.argv[3])
        coordinates = str(sys.argv[4])

        while sgl<0: sgl+=360.
        while sgl>360: sgl-=360.
        
        if sgb>90 or sgb<-90: sgb=float('exception')
        
        if (coordinates!='blank'):
            
            if len(sys.argv)>5:
                veldist = float(sys.argv[5])
                VD = float(sys.argv[6])
                sgl, sgb = interpgrid(sgl, sgb, cone=cone, coordinates=coordinates, veldist=veldist, VD=VD)
            else:
                 sgl, sgb = interpgrid(sgl, sgb, cone=cone, coordinates=coordinates)
            myCoords(sgl,sgb, addCord=True)
        else:
            interpgrid(sgl, sgb, cone=cone, coordinates=coordinates)
            myCoords(0,0, addCord=False)

        
    except:
        print '<tr><td>'
        print '<p><font color="red">&nbsp;&nbsp;Something went wrong !&nbsp;&nbsp;</font></p>'
        print '<p><font color="green">&nbsp;&nbsp;Check the input parameters.&nbsp;&nbsp;</font></p>'
        print '</td></tr>'
    

    print '</table>'


    













###############################################################
