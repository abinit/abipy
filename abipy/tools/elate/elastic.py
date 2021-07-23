# -*- coding: utf-8 -*-

import json
import math
import os
import platform
import random
import re
import sys
import time

from collections import OrderedDict
from io import StringIO
import requests

import numpy as np
from scipy import optimize


__author__ = "Romain Gaillac and François-Xavier Coudert"
__version__ = "2020.07.26"
__license__ = "MIT"


def removeHTMLTags(s):
    """Remove HTML tags, notably for use as page title"""
    return re.sub('<[^<]+?>', '', s)


def finishWebPage(outbuffer):
    """ Write the footer and finish the page """

    print('<div id="footer" class="content">')
    print('Code version: ' + __version__ + ' (running on Python ' + platform.python_version() + ')<br/>')
    print('<script type="text/javascript">var endTime = %g;' % time.perf_counter())
    print('document.write("Execution time: " + (endTime-startTime).toFixed(3) + " seconds<br/>");')
    print('if(typeof isOrtho !== \'undefined\') document.write("Specific (faster) code for orthorhombic case was used.");')
    print('</script></div>')
    print('</div>')
    print('</body></html>')
    return outbuffer.getvalue()


def writeHeader(outbuffer, title="Elastic Tensor Analysis"):
    """ Write the header of the HTML page """

    print("""
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html>
    <head>
        <title>%s</title>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
        <!--
        <link rel="stylesheet" type="text/css" href="/default.css" />
        -->
        %s
        <link rel="stylesheet" type="text/css" href="/assets/css/elate.css" />
        <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/jsxgraph@1.1.0/distrib/jsxgraph.css" />
        <script src="https://cdn.jsdelivr.net/npm/jsxgraph@1.1.0/distrib/jsxgraphcore.js"></script>
        <script src="http://cdn.plot.ly/plotly-latest.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
        </head>
    """ % (title, INTERNAL_CSS))


# printTitle writes the introduction of Elate
def printTitle(outbuffer, title="Elastic Tensor Analysis"):
    writeHeader(outbuffer, title)
    print("""
        <body>

        <div class="content">
        <h1><a href="/elate">ELATE: Elastic tensor analysis</a></h1>

        <p>Welcome to ELATE, the online tool for analysis of elastic tensors, developed by <b>Romain Gaillac</b> and <b><a
        href="http://coudert.name">François-Xavier Coudert</a></b> at <a href="http://www.chimie-paristech.fr/molsim/">CNRS / Chimie
        ParisTech</a>. <br/> If you use the software in published results (paper, conference, etc.), please cite the <a
        href="http://dx.doi.org/10.1088/0953-8984/28/27/275201">corresponding paper</a> (<em>J. Phys. Condens. Matter</em>, 2016, 28, 275201) and give the
        website URL.</p>

        <p>ELATE is <a href="https://github.com/fxcoudert/elate">open source software</a>. Any queries or comments are welcome at
    <script type="text/javascript">
    //<![CDATA[
    var c_="";for(var o5=0;o5<411;o5++)c_+=String.fromCharCode(("s%oz65j5>oJ.~~vs!Kt00}.~|}{\\"$s~%}!s0Kv#\\"wv<s!~tjjK{j5wo#zH}<j5s!z~qo6s~=u=i:00ikk>97a6!#|w<u!t{}vQ!o}Qsr?6F8G9:B8D9>@?7>a9!#|w<u!t{}vQ!o}QsrB67Dj59}qr$!s8#vq{wsw~;!oAA\\"wA#qsj5v!<~sozsq=6=A:u00970i0<ikk>a9!#|w<u!t{}vQ!o}QsrA69DDD>:E\\'7@<7s!z~qo6sjj==8:uN070j59j5jj.0|}}{\\"$}s#$0Kv#\\"wv<s!Ktj5jjj5jjL0\\'t14>O>>DBqI$}sr#!14>>>>BDqIwvw{sO~;!o\\"ws#vq14>>B>ID!t=JLo<j5s!z~qo6sO=u=0:705<!s~zoqs6=6<76<7=u:02@2?07<\\"$p\\"#!6?77".charCodeAt(o5)-(14)+0x3f)%(2*6+83)+64-32);document.write(eval(c_))
    //]]>
    </script>
        </p>
    """)


# 3D plot functions
################################################################################################

def write3DPlotData(dataX, dataY, dataZ, dataR, n, opacity=1.0):

    showcont = "true"
    if (opacity != 1.0):
        showcont = "false"
    if (n == 1):
        js = OrderedDict([
            ("x", dataX),
            ("y", dataY),
            ("z", dataZ),
            ("text", dataR),
            ("showscale", "false"),
            ("colorscale", "[[\'0\',\'rgb(22,136,51)\'],[\'0.125\',\'rgb(61,153,85)\'],[\'0.25\',\'rgb(121,178,136)\'],[\'0.375\',\'rgb(181,204,187)\'],[\'0.5\',\'rgb(195,230,200)\'],[\'0.625\',\'rgb(181,204,187)\'],[\'0.75\',\'rgb(121,178,136)\'],[\'0.875\',\'rgb(61,153,85)\'],[\'1\',\'rgb(22,136,51)\']]"),
            ("zsmooth", "'fast'"),
            ("type", "'surface'"),
            ("hoverinfo", "'text'"),
            ("opacity", opacity),
            ("contours", "{x :{ show:"+showcont+", color: 'rgb(192,192,192)'},y :{ show:"+showcont+", color: 'rgb(192,192,192)'},z :{ show:"+showcont+", color: 'rgb(192,192,192)'}}")
        ])

    if (n == 2):
        js = OrderedDict([
            ("x", dataX),
            ("y", dataY),
            ("z", dataZ),
            ("text", dataR),
            ("showscale", "false"),
            ("colorscale", "[[\'0\',\'rgb(180,4,38)\'],[\'0.125\',\'rgb(222,96,77)\'],[\'0.25\',\'rgb(244,154,123)\'],[\'0.375\',\'rgb(245,196,173)\'],[\'0.5\',\'rgb(246,216,201)\'],[\'0.625\',\'rgb(245,196,173)\'],[\'0.75\',\'rgb(244,154,123)\'],[\'0.875\',\'rgb(222,96,77)\'],[\'1\',\'rgb(180,4,38)\']]"),
            ("zsmooth", "'fast'"),
            ("type", "'surface'"),
            ("hoverinfo", "'text'"),
            ("opacity", opacity),
            ("contours", "{x :{ show:"+showcont+", color: 'rgb(192,192,192)'},y :{ show:"+showcont+", color: 'rgb(192,192,192)'},z :{ show:"+showcont+", color: 'rgb(192,192,192)'}}")
        ])

    if (n == 3):
        js = OrderedDict([
            ("x", dataX),
            ("y", dataY),
            ("z", dataZ),
            ("text", dataR),
            ("showscale", "false"),
            ("colorscale", "[[\'0\',\'rgb(59,76,192)\'],[\'0.125\',\'rgb(98,130,234)\'],[\'0.25\',\'rgb(141,176,254)\'],[\'0.375\',\'rgb(184,208,249)\'],[\'0.5\',\'rgb(207,223,250)\'],[\'0.625\',\'rgb(184,208,249)\'],[\'0.75\',\'rgb(141,176,254)\'],[\'0.875\',\'rgb(98,130,234)\'],[\'1\',\'rgb(59,76,192)\']]"),
            ("zsmooth", "'fast'"),
            ("type", "'surface'"),
            ("hoverinfo", "'text'"),
            ("opacity", opacity),
            ("contours", "{x :{ show:"+showcont+", color: 'rgb(192,192,192)'},y :{ show:"+showcont+", color: 'rgb(192,192,192)'},z :{ show:"+showcont+", color: 'rgb(192,192,192)'}}")
        ])

    print(json.dumps(js, indent=3).replace('\"', '') + ";")


def make3DPlot(func, legend='', width=600, height=600, npoints=200):

    str1 = legend.split("\'")[0]
    str2 = legend.split("\'")[1]

    u = np.linspace(0, np.pi, npoints)
    v = np.linspace(0, 2*np.pi, 2*npoints)
    r = np.zeros(len(u)*len(v))

    dataX = [[0.0 for i in range(len(v))] for j in range(len(u))]
    dataY = [[0.0 for i in range(len(v))] for j in range(len(u))]
    dataZ = [[0.0 for i in range(len(v))] for j in range(len(u))]
    dataR = [["0.0" for i in range(len(v))] for j in range(len(u))]

    count = 0
    for cu in range(len(u)):
        for cv in range(len(v)):
            r_tmp = func(u[cu], v[cv])
            z = r_tmp * np.cos(u[cu])
            x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
            y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
            dataX[cu][cv] = x
            dataY[cu][cv] = y
            dataZ[cu][cv] = z
            dataR[cu][cv] = "'E = "+str(float(int(10*r_tmp))/10.0)+" GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
            count = count+1

    i = random.randint(0, 100000)
    print('<div class="plot3D">')
    print('<div id="box%d" style="width: %dpx; height: %dpx; display:block;"></div>' % (i, width, height))
    print('</div>')
    print('<script type="text/javascript">')
    print("var trace =")
    write3DPlotData(dataX, dataY, dataZ, dataR, 1)
    print("var data = [trace]")
    print("var layout =")
    layout = {"title": "\'"+str1+"\\"+"\'"+str2+"\'", "width": "650", "height": "700", "autosize": "false", "autorange": "true", "margin": "{l: 65, r: 50, b: 65, t: 90}"}
    print(json.dumps(layout, indent=3).replace('\\\\', '\\').replace('\"', '') + ";")
    print("Plotly.newPlot('box%d',data,layout);" % (i))
    print('</script>')


def make3DPlotPosNeg(func, legend='', width=600, height=600, npoints=200):

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, 2*np.pi, 2*npoints)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  for cu in range(len(u)):
    for cv in range(len(v)):
      r_tmp = max(0, func(u[cu], v[cv]))
      z = r_tmp * np.cos(u[cu])
      x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
      y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
      dataX1[cu][cv] = x
      dataY1[cu][cv] = y
      dataZ1[cu][cv] = z
      dataR1[cu][cv] = "'"+"\u03B2 = "+str(float(int(10*r_tmp))/10.0)+" TPa'"+"+'-1'.sup()+"+"', \u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  for cu in range(len(u)):
    for cv in range(len(v)):
      r_tmp = max(0, -func(u[cu], v[cv]))
      z = r_tmp * np.cos(u[cu])
      x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
      y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
      dataX2[cu][cv] = x
      dataY2[cu][cv] = y
      dataZ2[cu][cv] = z
      dataR2[cu][cv] = "'"+"\u03B2 = -"+str(float(int(10*r_tmp))/10.0)+" TPa'"+"+'-1'.sup()+"+"', \u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  i = random.randint(0, 100000)
  print('<div class="plot3D">')
  print('<div id="box%d" style="width: %dpx; height: %dpx; display:block;"></div>' % (i, width, height))
  print('</div>')
  print('<script type="text/javascript">')
  print("var trace1 =")
  write3DPlotData(dataX1, dataY1, dataZ1, dataR1, 1)
  print("var trace2 =")
  write3DPlotData(dataX2, dataY2, dataZ2, dataR2, 2)
  print("var data = [trace1, trace2]")
  print("var layout =")
  layout = {"title": "\'"+legend+"\'", "width": "650", "height": "700", "autosize": "false", "autorange": "true", "margin": "{l: 65, r: 50, b: 65, t: 90}"}
  print(json.dumps(layout, indent=3).replace('\\\\', '\\').replace('\"', '') + ";")
  print("Plotly.newPlot('box%d',data,layout);" % (i))
  print('</script>')


def make3DPlot2(func, legend='', width=600, height=600, npoints=50):

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, np.pi, npoints)
  w = [v[i]+np.pi for i in range(1,len(v))]
  v = np.append(v, w)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  r = [0.0,0.0,np.pi/2.0,np.pi/2.0]
  for cu in range(len(u)):
    for cv in range(len(v)):

      r = func(u[cu],v[cv],r[2],r[3])
      z = np.cos(u[cu])
      x = np.sin(u[cu]) * np.cos(v[cv])
      y = np.sin(u[cu]) * np.sin(v[cv])

      r1_tmp = r[0]
      z1 = r1_tmp * z
      x1 = r1_tmp * x
      y1 = r1_tmp * y
      dataX1[cu][cv] = x1
      dataY1[cu][cv] = y1
      dataZ1[cu][cv] = z1
      dataR1[cu][cv] = "'"+"G'"+"+'min'.sub()+"+"' = "+str(float(int(10*r1_tmp))/10.0)+"GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"

      r2_tmp = r[1]
      z2 = r2_tmp * z
      x2 = r2_tmp * x
      y2 = r2_tmp * y
      dataX2[cu][cv] = x2
      dataY2[cu][cv] = y2
      dataZ2[cu][cv] = z2
      dataR2[cu][cv] = "'"+"G'"+"+'max'.sub()+"+"' = "+str(float(int(10*r1_tmp))/10.0)+"GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  i = random.randint(0, 100000)
  print('<div class="plot3D">')
  print('<div id="box%d" style="width: %dpx; height: %dpx; display:block;"></div>' % (i, width, height))
  print('</div>')
  print('<script type="text/javascript">')
  print("var trace1 =")
  write3DPlotData(dataX1, dataY1, dataZ1, dataR1, 1)
  print("var trace2 =")
  write3DPlotData(dataX2, dataY2, dataZ2, dataR2, 3, 0.5)
  print("var data = [trace1, trace2]")
  print("var layout =")
  layout = {"title": "\'"+legend+"\'", "width":"650", "height":"700" , "autosize":"false", "autorange":"true", "margin": "{l: 65, r: 50, b: 65, t: 90}"}
  print(json.dumps(layout, indent=3).replace('\\\\','\\').replace('\"','') + ";")
  print("Plotly.newPlot('box%d',data,layout);" % (i))
  print('</script>')


def make3DPlot3(func, legend='', width=600, height=600, npoints=50):

  str1 = legend.split("\'")[0]
  str2 = legend.split("\'")[1]

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, np.pi, npoints)
  w = [v[i]+np.pi for i in range(1,len(v))]
  v = np.append(v, w)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR3 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  r = [0.0, 0.0, 0.0, np.pi/2.0, np.pi/2.0]
  ruv = [[r for i in range(len(u))] for j in range(len(v))]
  for cu in range(len(u)):
    for cv in range(len(v)):
       ruv[cv][cu] = func(u[cu],v[cv],r[3],r[4])

  for cu in range(len(u)):
    for cv in range(len(v)):

      z = np.cos(u[cu])
      x = np.sin(u[cu]) * np.cos(v[cv])
      y = np.sin(u[cu]) * np.sin(v[cv])

      r = ruv[cv][cu]
      r1_tmp = r[0]
      dataX1[cu][cv] = r1_tmp * x
      dataY1[cu][cv] = r1_tmp * y
      dataZ1[cu][cv] = r1_tmp * z
      dataR1[cu][cv] = "'"+"\u03BD'"+"+'min'.sub()+"+"' = "+str(float(int(100*r1_tmp))/100.0)+", "+"\u03B8 = "+str(float(int(100*u[cu]*180/np.pi))/100.0)+"\u00B0, "+"\u03c6 = "+str(float(int(100*v[cv]*180/np.pi))/100.0)+"\u00B0'"

      r2_tmp = r[1]
      dataX2[cu][cv] = r2_tmp * x
      dataY2[cu][cv] = r2_tmp * y
      dataZ2[cu][cv] = r2_tmp * z
      dataR2[cu][cv] = "'"+"\u03BD'"+"+'min'.sub()+"+"' = "+str(float(int(100*r2_tmp))/100.0)+", "+"\u03B8 = "+str(float(int(100*u[cu]*180/np.pi))/100.0)+"\u00B0, "+"\u03c6 = "+str(float(int(100*v[cv]*180/np.pi))/100.0)+"\u00B0'"

      r3_tmp = r[2]
      dataX3[cu][cv] = r3_tmp * x
      dataY3[cu][cv] = r3_tmp * y
      dataZ3[cu][cv] = r3_tmp * z
      dataR3[cu][cv] = "'"+"\u03BD'"+"+'max'.sub()+"+"' = "+str(float(int(100*r3_tmp))/100.0)+", "+"\u03B8 = "+str(float(int(100*u[cu]*180/np.pi))/100.0)+"\u00B0, "+"\u03c6 = "+str(float(int(100*v[cv]*180/np.pi))/100.0)+"\u00B0'"
      count = count+1

  i = random.randint(0, 100000)
  print('<div class="plot3D">')
  print('<div id="box%d" style="width: %dpx; height: %dpx; display:block;"></div>' % (i, width, height))
  print('</div>')
  print('<script type="text/javascript">')
  print("var trace1 =")
  write3DPlotData(dataX1, dataY1, dataZ1, dataR1, 2, 0.5)
  print("var trace2 =")
  write3DPlotData(dataX2, dataY2, dataZ2, dataR2, 1, 1.0)
  print("var trace3 =")
  write3DPlotData(dataX3, dataY3, dataZ3, dataR3, 3, 0.5)
  print("var data = [trace1, trace2, trace3]")
  print("var layout =")
  layout = {"title": "\'"+str1+"\\"+"\'"+str2+"\'", "width":"650", "height":"700" , "autosize":"false", "autorange":"true", "margin": "{l: 65, r: 50, b: 65, t: 90}"}
  print(json.dumps(layout, indent=3).replace('\\\\','\\').replace('\"','') + ";")
  print("Plotly.newPlot('box%d',data,layout);" % (i))
  print('</script>')




# Polar plot functions
################################################################################################

def writePolarPlotData(dataX, dataY, suffix):
    """Write data for a polar plot, taking care of the center of inversion"""

    print("var dataX" + suffix + " = [")
    print((len(dataX) * "%.5f,") % tuple(dataX))
    print(((len(dataX)-1) * "%.5f," + "%.5f") % tuple(-dataX))
    print("];")
    print("var dataY" + suffix + " = [")
    print((len(dataX) * "%.5f,") % tuple(dataY))
    print(((len(dataX)-1) * "%.5f," + "%.5f") % tuple(-dataY))
    print("];")



def makePolarPlot(func, maxrad, legend='', p='xy', width=300, height=300, npoints=90, color='#009010', linewidth=2):

    i = random.randint(0, 100000)
    print('<div class="plot">')
    print('<div id="box%d" class="jxgbox" style="width: %dpx; height: %dpx; display:inline-block;"></div>' % (i, width, height))
    print('<br />%s</div>' % legend)
    print('<script type="text/javascript">')
    print('var b = JXG.JSXGraph.initBoard(\'box%d\', {boundingbox: [-%f, %f, %f, -%f], axis:true, showcopyright: 0});'
          % (i, maxrad, maxrad, maxrad, maxrad))

    u = np.linspace(0, np.pi, npoints)
    r = list(map(func, u))
    if (p=="xy"):
        x = r * np.cos(u)
        y = r * np.sin(u)
    else:
        y = r * np.cos(u)
        x = r * np.sin(u)

    writePolarPlotData (x, y, "")
    print("b.create('curve', [dataX,dataY], {strokeColor:'%s', strokeWidth: %d});" % (color, linewidth))
    print('</script>')

def makePolarPlotPosNeg(func, maxrad, legend='', p='xy', width=300, height=300, npoints=90, linewidth=2):
    i = random.randint(0, 100000)
    print('<div class="plot">')
    print('<div id="box%d" class="jxgbox" style="width: %dpx; height: %dpx; display:inline-block;"></div>' % (i, width, height))
    print('<br />%s</div>' % legend)
    print('<script type="text/javascript">')
    print('var b = JXG.JSXGraph.initBoard(\'box%d\', {boundingbox: [-%f, %f, %f, -%f], axis:true, showcopyright: 0});'
          % (i, maxrad, maxrad, maxrad, maxrad))

    u = np.linspace(0, np.pi, npoints)
    r = list(map(lambda x: max(0, func(x)), u))
    if (p=="xy"):
        x1 = r * np.cos(u)
        y1 = r * np.sin(u)
    else:
        y1 = r * np.cos(u)
        x1 = r * np.sin(u)
    r = list(map(lambda x: max(0, -func(x)), u))
    if (p=="xy"):
        x2 = r * np.cos(u)
        y2 = r * np.sin(u)
    else:
        y2 = r * np.cos(u)
        x2 = r * np.sin(u)

    writePolarPlotData (x1, y1, "1")
    writePolarPlotData (x2, y2, "2")
    print("b.create('curve', [dataX1,dataY1], {strokeColor:'green', strokeWidth: %d});" % (linewidth))
    print("b.create('curve', [dataX2,dataY2], {strokeColor:'red', strokeWidth: %d});" % (linewidth))
    print('</script>')

def makePolarPlot2(func, maxrad, legend='', p='xy', width=300, height=300, npoints=61, linewidth=2):
    i = random.randint(0, 100000)
    print('<div class="plot">')
    print('<div id="box%d" class="jxgbox" style="width: %dpx; height: %dpx; display:inline-block;"></div>' % (i, width, height))
    print('<br />%s</div>' % legend)
    print('<script type="text/javascript">')
    print('var b = JXG.JSXGraph.initBoard(\'box%d\', {boundingbox: [-%f, %f, %f, -%f], axis:true, showcopyright: 0});'
          % (i, maxrad, maxrad, maxrad, maxrad))

    u = np.linspace(0, np.pi, npoints)
    r = list(map(func, u))

    if (p=="xy"):
        x1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
        y1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
        x2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
        y2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
    else:
        y1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
        x1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
        y2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
        x2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])

    writePolarPlotData (x1, y1, "1")
    writePolarPlotData (x2, y2, "2")
    print("b.create('curve', [dataX1,dataY1], {strokeColor:'green', strokeWidth: %d});" % (linewidth))
    print("b.create('curve', [dataX2,dataY2], {strokeColor:'blue', strokeWidth: %d});" % (linewidth))
    print('</script>')

def makePolarPlot3(func, maxrad, legend='', p='xy', width=300, height=300, npoints=61, linewidth=2):
    i = random.randint(0, 100000)
    print('<div class="plot">')
    print('<div id="box%d" class="jxgbox" style="width: %dpx; height: %dpx; display:inline-block;"></div>' % (i, width, height))
    print('<br />%s</div>' % legend)
    print('<script type="text/javascript">')
    print('var b = JXG.JSXGraph.initBoard(\'box%d\', {boundingbox: [-%f, %f, %f, -%f], axis:true, showcopyright: 0});'
          % (i, maxrad, maxrad, maxrad, maxrad))

    u = np.linspace(0, np.pi, npoints)
    r = list(map(func, u))

    if (p=="xy"):
        x1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
        y1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
        x2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
        y2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
        x3 = np.array([ ir[2] * np.cos(iu) for ir, iu in zip(r,u) ])
        y3 = np.array([ ir[2] * np.sin(iu) for ir, iu in zip(r,u) ])
    else:
        y1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
        x1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
        y2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
        x2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
        y3 = np.array([ ir[2] * np.cos(iu) for ir, iu in zip(r,u) ])
        x3 = np.array([ ir[2] * np.sin(iu) for ir, iu in zip(r,u) ])

    writePolarPlotData (x1, y1, "1")
    writePolarPlotData (x2, y2, "2")
    writePolarPlotData (x3, y3, "3")
    print("b.create('curve', [dataX1,dataY1], {strokeColor:'red', strokeWidth: %d});" % (linewidth))
    print("b.create('curve', [dataX2,dataY2], {strokeColor:'green', strokeWidth: %d});" % (linewidth))
    print("b.create('curve', [dataX3,dataY3], {strokeColor:'blue', strokeWidth: %d});" % (linewidth))
    print('</script>')


################################################################################################

def dirVec(theta, phi):
    return [ math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta) ]

def dirVec1(theta, phi, chi):
    return [ math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta) ]

def dirVec2(theta, phi, chi):
    return [ math.cos(theta)*math.cos(phi)*math.cos(chi) - math.sin(phi)*math.sin(chi),
             math.cos(theta)*math.sin(phi)*math.cos(chi) + math.cos(phi)*math.sin(chi),
             - math.sin(theta)*math.cos(chi) ]


# Functions to minimize/maximize
def minimize(func, dim):
    if dim == 2:
        r = ((0, np.pi), (0, np.pi))
        n = 25
    elif dim == 3:
        r = ((0, np.pi), (0, np.pi), (0, np.pi))
        n = 10

    # TODO -- try basin hopping or annealing
    return optimize.brute(func, r, Ns = n, full_output = True, finish = optimize.fmin)[0:2]


def maximize(func, dim):
    res = minimize(lambda x: -func(x), dim)
    return (res[0], -res[1])


class Elastic:
    """An elastic tensor, along with methods to access it"""

    def __init__(self, s):
        """Initialize the elastic tensor from a string"""

        if not s:
            raise ValueError("no matrix was provided")

        # Argument can be a 6-line string, a list of list, or a string representation of the list of list
        try:
            if type(json.loads(s)) == list: s = json.loads(s)
        except:
            pass

        if type(s) == str:
            # Remove braces and pipes
            s = s.replace("|", " ").replace("(", " ").replace(")", " ")

            # Remove empty lines
            lines = [line for line in s.split('\n') if line.strip()]
            if len(lines) != 6:
                raise ValueError("should have six rows")

            # Convert to float
            try:
                mat = [list(map(float, line.split())) for line in lines]
            except:
                raise ValueError("not all entries are numbers")
        elif type(s) == list:
            # If we already have a list, simply use it
            mat = s
        else:
            raise ValueError("invalid argument as matrix")

        # Make it into a square matrix
        mat = np.array(mat)
        if mat.shape != (6,6):
            # Is it upper triangular?
            if list(map(len, mat)) == [6,5,4,3,2,1]:
                mat = [ [0]*i + mat[i] for i in range(6) ]
                mat = np.array(mat)

        # Is it lower triangular?
        if list(map(len, mat)) == [1,2,3,4,5,6]:
            mat = [ mat[i] + [0]*(5-i) for i in range(6) ]
            mat = np.array(mat)

        if mat.shape != (6,6):
            raise ValueError("should be a square matrix")

        # Check that is is symmetric, or make it symmetric
        if np.linalg.norm(np.tril(mat, -1)) == 0:
            mat = mat + np.triu(mat, 1).transpose()
        if np.linalg.norm(np.triu(mat, 1)) == 0:
            mat = mat + np.tril(mat, -1).transpose()
        if np.linalg.norm(mat - mat.transpose()) > 1e-3:
            raise ValueError("should be symmetric, or triangular")
        elif np.linalg.norm(mat - mat.transpose()) > 0:
            mat = 0.5 * (mat + mat.transpose())

        # Store it
        self.CVoigt = mat

        # Put it in a more useful representation
        try:
            self.SVoigt = np.linalg.inv(self.CVoigt)
        except:
            raise ValueError("matrix is singular")

        VoigtMat = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
        def SVoigtCoeff(p,q): return 1. / ((1+p//3)*(1+q//3))

        self.Smat = [[[[ SVoigtCoeff(VoigtMat[i][j], VoigtMat[k][l]) * self.SVoigt[VoigtMat[i][j]][VoigtMat[k][l]]
                         for i in range(3) ] for j in range(3) ] for k in range(3) ] for l in range(3) ]
        return

    def isOrthorhombic(self):
        def iszero(x): return (abs(x) < 1.e-3)
        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5]))

    def isCubic(self):
        def iszero(x): return (abs(x) < 1.e-3)
        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[3][3] - self.CVoigt[4][4]) and iszero(self.CVoigt[3][3] - self.CVoigt[5][5])
                and iszero(self.CVoigt[0][1] - self.CVoigt[0][2]) and iszero(self.CVoigt[0][1] - self.CVoigt[1][2]))

    def Young(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

    def Young_2(self,x,y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

    def LC(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

    def LC_2(self, x, y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

    def shear(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r = sum([ a[i]*b[j]*a[k]*b[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/(4*r)

    def Poisson(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r1 = sum([ a[i]*a[j]*b[k]*b[l] * self.Smat[i][j][k][l]
                   for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        r2 = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                   for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return -r1/r2

    def averages(self):
        A = (self.CVoigt[0][0] + self.CVoigt[1][1] + self.CVoigt[2][2]) / 3
        B = (self.CVoigt[1][2] + self.CVoigt[0][2] + self.CVoigt[0][1]) / 3
        C = (self.CVoigt[3][3] + self.CVoigt[4][4] + self.CVoigt[5][5]) / 3
        a = (self.SVoigt[0][0] + self.SVoigt[1][1] + self.SVoigt[2][2]) / 3
        b = (self.SVoigt[1][2] + self.SVoigt[0][2] + self.SVoigt[0][1]) / 3
        c = (self.SVoigt[3][3] + self.SVoigt[4][4] + self.SVoigt[5][5]) / 3

        KV = (A + 2*B) / 3
        GV = (A - B + 3*C) / 5

        KR = 1 / (3*a + 6*b)
        GR = 5 / (4*a - 4*b + 3*c)

        KH = (KV + KR) / 2
        GH = (GV + GR) / 2

        return [ [KV, 1/(1/(3*GV) + 1/(9*KV)), GV, (1 - 3*GV/(3*KV+GV))/2],
                 [KR, 1/(1/(3*GR) + 1/(9*KR)), GR, (1 - 3*GR/(3*KR+GR))/2],
                 [KH, 1/(1/(3*GH) + 1/(9*KH)), GH, (1 - 3*GH/(3*KH+GH))/2] ]

    def shear2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.shear([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.shear([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun))

    def shear3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.005
        def func1(z): return self.shear([x, y, z])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.shear([x, y, z])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun), float(r1.x), float(r2.x))

    def Poisson2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.Poisson([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.Poisson([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun))

    def poisson3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.005
        def func1(z): return self.Poisson([x, y, z])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.Poisson([x, y, z])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun), float(r1.x), float(r2.x))


class ElasticOrtho(Elastic):
    """An elastic tensor, for the specific case of an orthorhombic system"""

    def __init__(self, arg):
        """Initialize from a matrix, or from an Elastic object"""
        if type(arg) == str:
            Elastic.__init__(self, arg)
        elif isinstance(arg, Elastic):
            self.CVoigt = arg.CVoigt
            self.SVoigt = arg.SVoigt
            self.Smat = arg.Smat
        else:
            raise TypeError("ElasticOrtho constructor argument should be string or Elastic object")

    def Young(self, x):
        ct2 = math.cos(x[0])**2
        st2 = 1 - ct2
        cf2 = math.cos(x[1])**2
        sf2 = 1 - cf2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1/(ct2**2*s33 + 2*cf2*ct2*s13*st2 + cf2*ct2*s55*st2 + 2*ct2*s23*sf2*st2 + ct2*s44*sf2*st2 + cf2**2*s11*st2**2 + 2*cf2*s12*sf2*st2**2 + cf2*s66*sf2*st2**2 + s22*sf2**2*st2**2)

    def LC(self, x):
        ct2 = math.cos(x[0])**2
        cf2 = math.cos(x[1])**2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1000 * (ct2 * (s13 + s23 + s33) + (cf2 * (s11 + s12 + s13) + (s12 + s22 + s23) * (1 - cf2)) * (1 - ct2))

    def shear(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        sf2 = sf*sf
        cx = math.cos(x[2])
        cx2 = cx*cx
        sx = math.sin(x[2])
        sx2 = 1 - cx2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        r = (
            ct2*ct2*cx2*s44*sf2 + cx2*s44*sf2*st2*st2 + 4*cf**3*ct*cx*(-2*s11 + 2*s12 + s66)*sf*st2*sx
            + 2*cf*ct*cx*sf*(ct2*(s44 - s55) + (4*s13 - 4*s23 - s44 + s55 - 4*s12*sf2 + 4*s22*sf2 - 2*s66*sf2)*st2)*sx
            + s66*sf2*sf2*st2*sx2 + cf**4*st2*(4*ct2*cx2*s11 + s66*sx2)
            + ct2*(2*cx2*(2*s33 + sf2*(-4*s23 - s44 + 2*s22*sf2))*st2 + s55*sf2*sx2)
            + cf**2*(ct2*ct2*cx2*s55 + ct2*(-2*cx2*(4*s13 + s55 - 2*(2*s12 + s66)*sf2)*st2 + s44*sx2)
                     + st2*(cx2*s55*st2 + 2*(2*s11 - 4*s12 + 2*s22 - s66)*sf2*sx2))
            )
        return 1/r

    def Poisson(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        cx = math.cos(x[2])
        sx = math.sin(x[2])
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]

        return (
    (-(ct**2*cx**2*s33*st2) - cf**2*cx**2*s13*st2*st2 - cx**2*s23*sf**2*st2*st2 + ct*cx*s44*sf*st2*(ct*cx*sf + cf*sx) -
            ct**2*s23*(ct*cx*sf + cf*sx)**2 - cf**2*s12*st2*(ct*cx*sf + cf*sx)**2 - s22*sf**2*st2*(ct*cx*sf + cf*sx)**2 +
            cf*ct*cx*s55*st2*(cf*ct*cx - sf*sx) - cf*s66*sf*st2*(ct*cx*sf + cf*sx)*(cf*ct*cx - sf*sx) -
            ct**2*s13*(cf*ct*cx - sf*sx)**2 - cf**2*s11*st2*(cf*ct*cx - sf*sx)**2 - s12*sf**2*st2*(cf*ct*cx - sf*sx)**2)/
            (ct**4*s33 + 2*cf**2*ct**2*s13*st2 + cf**2*ct**2*s55*st2 + 2*ct**2*s23*sf**2*st2 + ct**2*s44*sf**2*st2 +
            cf**4*s11*st2*st2 + 2*cf**2*s12*sf**2*st2*st2 + cf**2*s66*sf**2*st2*st2 + s22*sf**4*st2*st2)
        )


################################################################################################

# Materials Project URL
urlBase = 'https://www.materialsproject.org/rest'


def queryMaterials(query, mapiKey):
    """Return a list of material IDs for a given query string"""

    # If the query is a material ID, return it
    if query[0:3] == "mp-": return [query]

    try:
        r = requests.get(f'{urlBase}/v2/materials/{query}/mids', headers={"X-API-KEY": mapiKey})
        resp = r.json()
    except Exception as e:
        print(str(e), file=sys.stderr)
        return []

    if (not resp['valid_response']): return []
    return resp['response']


def queryElasticityV2(mat, mapiKey):
    """Return elastic properties for a given material ID, using V2 MAPI"""

    data = { 'criteria': '{"task_id": "' + mat + '"}',
             'properties': '["formula", "pretty_formula", "material_id", "elasticity"]',
             'API_KEY': mapiKey }
    try:
        r = requests.post(f'{urlBase}/v2/query', data)
        resp = r.json()
    except Exception as e:
        print(str(e), file=sys.stderr)
        return None

    if not resp["valid_response"]: return None
    if len(resp["response"]) > 1: raise(Exception("Multiple results returned"))
    if len(resp["response"]) == 0: return None
    return resp["response"][0]


def ELATE_MaterialsProject(query, mapiKey):
  """Call ELATE with a query from the Materials Project"""

  # If we were directly given a material ID, or there is a simple match
  materials = queryMaterials(query, mapiKey)
  if len(materials) == 1:
    r = queryElasticityV2(query, mapiKey)

    if r and 'elasticity' in r:
      tensor = r["elasticity"]["elastic_tensor"]
      return ELATE(tensor, '%s (Materials Project id <a href="%s%s" target="_blank">%s</a>)' % (r["pretty_formula"], "https://www.materialsproject.org/materials/", r["material_id"], r["material_id"]))

  # Otherwise, run the MP query, list the matches and let the user choose
  sys.stdout = outbuffer = StringIO()
  printTitle(outbuffer, "ELATE: Elastic tensor analysis")
  print('<h2>Query from the Materials Project database</h2>')

  # Either there was no match, or a single match with no elastic data
  if len(materials) <= 1:
    print("""<p>
            Your query for <tt style="background-color: #e0e0e0;">%s</tt> from the <a href="https://materialsproject.org">Materials Project</a> database
            has returned a total of zero result. Or is it zero results? In any case, we are very sorry.</p>

            <p>If you wish, you can try another query here:
            <form name="elastic" action="/elate/mp" method="get">
              <input type="text" name="query" style="font-family:sans-serif; width: 20em;">
              <input type="submit" style="font-size: 100%%; color: #b02020;" value="Submit query">
            </form>
            or go back to our <a href="/elate">main page</a>.
            </p>""" % (query))
    return finishWebPage(outbuffer)


  print("""<p>Your query for <tt style="background-color: #e0e0e0;">%s</tt> from the <a href="https://materialsproject.org">Materials Project</a> database
           has returned %d results.""" % (query, len(materials)))

  if len(materials) > 10:
    materials = materials[0:10]
    print("Below is a table of the 10 first matches.")

  print("<table><tr><th>Identifier</th><th>Formula</th><th>Elastic data</th></tr>")
  for mat in materials:
    r = queryElasticityV2(mat, mapiKey)
    print('<tr><td><a href="https://www.materialsproject.org/materials/%s" target="_blank">%s</a></td><td>%s</td>' % (mat, mat, r["pretty_formula"]))
    if "elasticity" in r and r["elasticity"]:
      print('<td>Elastic data available, <a href="/elate/mp?%s" target="_blank">perform analysis</a></td></tr>' % (mat))
    else:
      print('<td>No elastic data available</td></tr>')
  print("</table>")

  return finishWebPage(outbuffer)


def ELATE(matrix, sysname):
  """ELATE performs the calculation and plots every property in 2D"""

  # Redirect output to out string buffer
  sys.stdout = outbuffer = StringIO()

  # Start timing
  print('<script type="text/javascript">var startTime = %g</script>' % time.perf_counter())

  printTitle(outbuffer, "Elastic analysis of " + removeHTMLTags(sysname))

  try:
    elas = Elastic(matrix)
  except ValueError as e:
    print('<div class="error">Invalid stiffness matrix: ')
    print(e.args[0])
    if matrix:
      print('<pre>' + str(matrix) + '</pre>')

    print('</div>')
    print('<input action="action" type="button" value="Go back" onclick="window.history.go(-1); return false;" />')
    return finishWebPage(outbuffer)

  if elas.isOrthorhombic():
    elas = ElasticOrtho(elas)
    print('<script type="text/javascript">var isOrtho = 1;</script>')

  print('<h2>Summary of the properties</h2>')

  print('<h3>Input: stiffness matrix (coefficients in GPa) of %s</h3>' % (sysname))
  print('<pre>')
  for i in range(6):
    print(("   " + 6*"%7.5g  ") % tuple(elas.CVoigt[i]))
  print('</pre>')

  avg = elas.averages()
  print('<h3>Average properties</h3>')

  print("<table><tr><th>Averaging scheme</th><th>Bulk modulus</th><th>Young's modulus</th><th>Shear modulus</th><th>Poisson's ratio</th></tr>")
  print(('<tr><td>Voigt</td><td><em>K</em><sub>V</sub> = %7.5g GPa</td><td><em>E</em><sub>V</sub> = %7.5g GPa</td>'
          + '<td><em>G</em><sub>V</sub> = %7.5g GPa</td><td><em>&nu;</em><sub>V</sub> = %.5g</td></tr>')
        % tuple(avg[0]))
  print(('<tr><td>Reuss</td><td><em>K</em><sub>R</sub> = %7.5g GPa</td><td><em>E</em><sub>R</sub> = %7.5g GPa</td>'
          + '<td><em>G</em><sub>R</sub> = %7.5g GPa</td><td><em>&nu;</em><sub>R</sub> = %.5g</td></tr>')
        % tuple(avg[1]))
  print(('<tr><td>Hill</td><td><em>K</em><sub>H</sub> = %7.5g GPa</td><td><em>E</em><sub>H</sub> = %7.5g GPa</td>'
          + '<td><em>G</em><sub>H</sub> = %7.5g GPa</td><td><em>&nu;</em><sub>H</sub> = %.5g</td></tr>')
        % tuple(avg[2]))
  print('</table>')


  print('''<h3>Eigenvalues of the stiffness matrix</h3>
  <table><tr>
  <th>&lambda;<sub>1</sub></th>
  <th>&lambda;<sub>2</sub></th>
  <th>&lambda;<sub>3</sub></th>
  <th>&lambda;<sub>4</sub></th>
  <th>&lambda;<sub>5</sub></th>
  <th>&lambda;<sub>6</sub></th>
  </tr><tr>''')
  eigenval = sorted(np.linalg.eig(elas.CVoigt)[0])
  print((6*'<td>%7.5g GPa</td>') % tuple(eigenval))
  print('</tr></table>')

  if eigenval[0] <= 0:
    print('<div class="error">Stiffness matrix is not definite positive, crystal is mechanically unstable<br/>')
    print('No further analysis will be performed.</div>')
    return finishWebPage(outbuffer)


  minE = minimize(elas.Young, 2)
  maxE = maximize(elas.Young, 2)
  minLC = minimize(elas.LC, 2)
  maxLC = maximize(elas.LC, 2)
  minG = minimize(elas.shear, 3)
  maxG = maximize(elas.shear, 3)
  minNu = minimize(elas.Poisson, 3)
  maxNu = maximize(elas.Poisson, 3)

  print("""<h3>Variations of the elastic moduli</h3>
            <table>
            <tr><td></td><th colspan="2">Young\'s modulus</th><th colspan="2">Linear compressibility</th>
            <th colspan="2">Shear modulus</th><th colspan="2">Poisson\'s ratio</th><th></th></tr>
            <td></td><th><em>E</em><sub>min</sub></th><th><em>E</em><sub>max</sub></th>
            <th>&beta;<sub>min</sub></th><th>&beta;<sub>max</sub></th><th><em>G</em><sub>min</sub></th><th><em>G</em><sub>max</sub></th>
            <th>&nu;<sub>min</sub></th><th>&nu;<sub>max</sub></th><th></th></tr>""")

  print(('<tr><td>Value</td><td>%8.5g GPa</td><td>%8.5g GPa</td>'
        + '<td>%8.5g TPa<sup>&ndash;1</sup></td><td>%8.5g TPa<sup>&ndash;1</sup></td>'
        + '<td>%8.5g GPa</td><td>%8.5g GPa</td>'
        + '<td>%.5g</td><td>%.5g</td><td>Value</td></tr>') % (minE[1], maxE[1], minLC[1], maxLC[1], minG[1], maxG[1], minNu[1], maxNu[1]))

  anisE = '%8.4g' % (maxE[1]/minE[1])
  anisLC = ('%8.4f' % (maxLC[1]/minLC[1])) if minLC[1] > 0 else "&infin;"
  anisG = '%8.4g' % (maxG[1]/minG[1])
  anisNu = ('%8.4f' % (maxNu[1]/minNu[1])) if minNu[1]*maxNu[1] > 0 else "&infin;"
  print(('<tr><td>Anisotropy</td>' + 4 * '<td colspan="2">%s</td>'
        + '<td>Anisotropy</td></tr>') % (anisE, anisLC, anisG, anisNu))

  print('<tr><td>Axis</td>')
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec(*minE[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec(*maxE[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec(*minLC[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec(*maxLC[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec1(*minG[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec1(*maxG[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec1(*minNu[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec1(*maxNu[0])))
  print('<td>Axis</td></tr>')

  print('<tr><td></td><td></td><td></td><td></td><td></td>')
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec2(*minG[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec2(*maxG[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec2(*minNu[0])))
  print('<td>%.4f<br />%.4f<br />%.4f</td>' % tuple(dirVec2(*maxNu[0])))
  print('<td>Second axis</td></tr></table>')

  print("<h2>Spatial dependence of Young's modulus</h2>")
  print("""<form id="elastic" action="/wait3D" method="post" target="_blank">
            <textarea name="matrix" style="display: none;">%s</textarea>
            <textarea name="sysname" style="display: none;">%s</textarea>
            <textarea name="job" style="display: none;">%s</textarea>
            <br /><input type="submit" style="font-size: 100%%; color: #b02020;" value="Visualize in 3D">
           </form>""" % (matrix, sysname, "young"))
  m = 1.2 * maxE[1]
  makePolarPlot(lambda x: elas.Young([np.pi / 2, x]), m, "Young's modulus in (xy) plane", "xy")
  makePolarPlot(lambda x: elas.Young([x, 0]), m, "Young's modulus in (xz) plane", "xz")
  makePolarPlot(lambda x: elas.Young([x, np.pi / 2]), m, "Young's modulus in (yz) plane", "yz")

  print("<h2>Spatial dependence of linear compressibility</h2>")
  print("""<form id="elastic" action="/wait3D" method="post" target="_blank">
            <textarea name="matrix" style="display: none;">%s</textarea>
            <textarea name="sysname" style="display: none;">%s</textarea>
            <textarea name="job" style="display: none;">%s</textarea>
            <br /><input type="submit" style="font-size: 100%%; color: #b02020;" value="Visualize in 3D">
           </form>""" % (matrix, sysname, "lc"))
  m = 1.2 * max(maxLC[1], abs(minLC[1]))
  makePolarPlotPosNeg(lambda x: elas.LC([np.pi / 2, x]), m, "linear compressibility in (xy) plane", "xy")
  makePolarPlotPosNeg(lambda x: elas.LC([x, 0]), m, "linear compressibility in (xz) plane", "xz")
  makePolarPlotPosNeg(lambda x: elas.LC([x, np.pi / 2]), m, "linear compressibility in (yz) plane", "yz")

  print("<h2>Spatial dependence of shear modulus</h2>")
  print("""<form id="elastic" action="/wait3D" method="post" target="_blank">
            <textarea name="matrix" style="display: none;">%s</textarea>
            <textarea name="sysname" style="display: none;">%s</textarea>
            <textarea name="job" style="display: none;">%s</textarea>
            <br /><input type="submit" style="font-size: 100%%; color: #b02020;" value="Visualize in 3D">
           </form>""" % (matrix, sysname, "shear"))
  m = 1.2 * maxG[1]
  makePolarPlot2(lambda x: elas.shear2D([np.pi / 2, x]), m, "Shear modulus in (xy) plane", "xy")
  makePolarPlot2(lambda x: elas.shear2D([x, 0]), m, "Shear modulus in (xz) plane", "xz")
  makePolarPlot2(lambda x: elas.shear2D([x, np.pi / 2]), m, "Shear modulus in (yz) plane", "yz")

  print("<h2>Spatial dependence of Poisson's ratio</h2>")
  print("""<form id="elastic" action="/wait3D" method="post" target="_blank">
            <textarea name="matrix" style="display: none;">%s</textarea>
            <textarea name="sysname" style="display: none;">%s</textarea>
            <textarea name="job" style="display: none;">%s</textarea>
            <br /><input type="submit" style="font-size: 100%%; color: #b02020;" value="Visualize in 3D">
           </form>""" % (matrix, sysname, "poisson"))
  m = 1.2 * max(abs(maxNu[1]), abs(minNu[1]))
  makePolarPlot3(lambda x: elas.Poisson2D([np.pi / 2, x]), m, "Poisson's ratio in (xy) plane", "xy")
  makePolarPlot3(lambda x: elas.Poisson2D([x, 0]), m, "Poisson's ratio in (xz) plane", "xz")
  makePolarPlot3(lambda x: elas.Poisson2D([x, np.pi / 2]), m, "Poisson's ratio in (yz) plane", "yz")

  print("</div>")
  return finishWebPage(outbuffer)


def wait3D(matrix, sysname, job):
    """Display a waiting page while we calculate a 3D plot"""

    sys.stdout = outbuffer = StringIO()
    writeHeader(outbuffer, "Young 3D for " + removeHTMLTags(sysname))

    print("""
    <div class="content">
    <img src="/loading.gif" alt="[loading]" />
    <p>Please wait while your 3D graph is loading… (it can take from 15 seconds up to a minute)</p>
    """)

    # Pass arguments
    print("""
    <form id="elastic" action="/plot3D" method="post" style="display: none;">
        <textarea name="matrix">%s</textarea>
        <textarea name="sysname">%s</textarea>
        <textarea name="job">%s</textarea>
        <input type="submit" value="">
    </form>""" % (matrix, sysname, job))

    # Reload immediately
    print("""
    <script type="text/javascript">
        window.onload = function(){
        setTimeout(function () {
            document.getElementById("elastic").submit();
                }, 100);
            };
    </script>""")

    return finishWebPage(outbuffer)


def plot3D(matrix, sysname, job):
    """Display a 3D plot"""

    # Dispatch to the specific function depending on type
    functions = {'young': YOUNG3D, 'lc': LC3D, 'shear': SHEAR3D, 'poisson': POISSON3D}
    return functions[job](matrix, sysname)


# ELATE : basic usage of the tool, only 2D plots
# YOUNG3D : visualize Young's modulus in 3D
# LC3D : visualize Linear compressiblity in 3D
# SHEAR3D : visualize Shear modulus in 3D
# POISSON3D : visualize Poisson ratio in 3D
################################################################################################


def YOUNG3D(matrix, sysname):

    sys.stdout = outbuffer = StringIO()
    writeHeader(outbuffer, "Young 3D for " + removeHTMLTags(sysname))

    # Start timing
    print('<script type="text/javascript">var startTime = %g</script>' % time.perf_counter())
    print('<div class="content">')

    print("<h1> 3D Visualization of Young's modulus </h1>")
    elas = Elastic(matrix)
    if elas.isOrthorhombic():
        elas = ElasticOrtho(elas)
        print('<script type="text/javascript">var isOrtho = 1;</script>')

    make3DPlot(lambda x, y: elas.Young_2(x, y), "Young's modulus")

    print('<h3>Input: stiffness matrix (coefficients in GPa) of %s</h3>' % (sysname))
    print('<pre>')
    for i in range(6):
        print(("   " + 6 * "%7.5g  ") % tuple(elas.CVoigt[i]))
    print('</pre></div>')
    return finishWebPage(outbuffer)


def LC3D(matrix, sysname):

    sys.stdout = outbuffer = StringIO()
    writeHeader(outbuffer, "LC 3D for " + removeHTMLTags(sysname))

    # Start timing
    print('<script type="text/javascript">var startTime = %g</script>' % time.perf_counter())
    print('<div class="content">')

    print("<h1> 3D Visualization of Linear compressiblity </h1>")
    elas = Elastic(matrix)
    if elas.isOrthorhombic():
        elas = ElasticOrtho(elas)
        print('<script type="text/javascript">var isOrtho = 1;</script>')

    make3DPlotPosNeg(lambda x, y: elas.LC_2(x, y), "Linear compressiblity")

    print('<h3>Input: stiffness matrix (coefficients in GPa) of %s</h3>' % (sysname))
    print('<pre>')
    for i in range(6):
        print(("   " + 6 * "%7.5g  ") % tuple(elas.CVoigt[i]))
    print('</pre></div>')
    return finishWebPage(outbuffer)


def SHEAR3D(matrix, sysname):

    sys.stdout = outbuffer = StringIO()
    writeHeader(outbuffer, "Shear 3D for " + removeHTMLTags(sysname))

    # Start timing
    print('<script type="text/javascript">var startTime = %g</script>' % time.perf_counter())
    print('<div class="content">')

    print("<h1> 3D Visualization of Shear modulus </h1>")
    elas = Elastic(matrix)
    if elas.isOrthorhombic():
        elas = ElasticOrtho(elas)
        print('<script type="text/javascript">var isOrtho = 1;</script>')

    make3DPlot2(lambda x, y, g1, g2: elas.shear3D(x, y, g1, g2), "Shear modulus")

    print('<h3>Input: stiffness matrix (coefficients in GPa) of %s</h3>' % (sysname))
    print('<pre>')
    for i in range(6):
        print(("   " + 6 * "%7.5g  ") % tuple(elas.CVoigt[i]))
    print('</pre></div>')
    return finishWebPage(outbuffer)


def POISSON3D(matrix, sysname):

    sys.stdout = outbuffer = StringIO()
    writeHeader(outbuffer, "Poisson 3D for " + removeHTMLTags(sysname))

    # Start timing
    print('<script type="text/javascript">var startTime = %g</script>' % time.perf_counter())
    print('<div class="content">')

    print("<h1> 3D Visualization of Poisson's ratio </h1>")
    elas = Elastic(matrix)
    if elas.isOrthorhombic():
        elas = ElasticOrtho(elas)
        print('<script type="text/javascript">var isOrtho = 1;</script>')

    make3DPlot3(lambda x, y, g1, g2: elas.poisson3D(x, y, g1, g2), "Poisson's ratio")

    print('<h3>Input: stiffness matrix (coefficients in GPa) of %s</h3>' % (sysname))
    print('<pre>')
    for i in range(6):
        print(("   " + 6 * "%7.5g  ") % tuple(elas.CVoigt[i]))
    print('</pre></div>')
    return finishWebPage(outbuffer)


# Added by MG
INTERNAL_CSS = """
<style>
a {
  color: #3030d0;
}

h1 a {
  text-decoration: none;
  color: #b02020;
}

table, th, td
{
  padding: 10px;
  text-align: center;
  empty-cells: hide;
}

th { background-color: #FFDD88; }
tr:nth-child(even) { background: #F0F0F0; }
tr:nth-child(odd)  { background: #F0F0FF; }

div.error {
  margin-top: 30px;
  padding: 10px;
  color: red;
  background-color: #ffd0d0;
  font-weight: bold;
}

div.plot {
  display: inline-block;
  margin: 10px;
  font-size: 80%;
  color: grey;
  text-align: center;
}

div.plot3D {
  display: block;
  font-size: 80%;
  color: grey;
  text-align: center;
  padding-top: 50px;
  padding-bottom: 50px;
}

div#footer {
  margin-top: 50px;
  font-size: 70%;
  padding-top: 10px;
  border-top: 1px solid grey;
}
</style>
"""

