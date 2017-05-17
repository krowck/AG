# Get this figure: fig = py.get_figure("https://plot.ly/~seanmac06/6/")
# Get this figure's data: data = py.get_figure("https://plot.ly/~seanmac06/6/").get_data()
# Add data to this figure: py.plot(Data([Scatter(x=[1, 2], y=[2, 3])]), filename ="Dendrogram", fileopt="extend")
# Get y data of first trace: y1 = py.get_figure("https://plot.ly/~seanmac06/6/").get_data()[0]["y"]

# Get figure documentation: https://plot.ly/python/get-requests/
# Add data documentation: https://plot.ly/python/file-options/

# If you're using unicode in your file, you may need to specify the encoding.
# You can reproduce this figure in Python with the following code!

# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *
trace1 = {
  "x": [1, 1, 2, 2], 
  "y": [0, 0.8, 0.8, 0], 
  "marker": {"color": "rgb(61,153,112)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "9fe5ba", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace2 = {
  "x": [1, 1, 2, 2], 
  "y": [0, 0.8, 0.8, 0], 
  "marker": {"color": "rgb(61,153,112)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "db5284", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace3 = {
  "x": [4, 4, 5, 5], 
  "y": [0, 2.7, 2.7, 0], 
  "marker": {"color": "rgb(255,65,54)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "bff71b", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace4 = {
  "x": [1.5, 1.5, 4.5, 4.5], 
  "y": [0.8, 3.8, 3.8, 3.8], 
  "marker": {"color": "rgb(255,65,54)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "9ed573", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace5 = {
  "x": [3, 3, 3, 3], 
  "y": [0, 3.8, 3.8, 0], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "fe92fc", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace6 = {
  "x": [3, 3, 2, 2], 
  "y": [0, 3.8, 3.8, 0], 
  "marker": {"color": "rgb(35,205,205)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "2d03e7", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace7 = {
  "x": [4, 4, 3, 3], 
  "y": [0, 5.7, 5.7, 0], 
  "marker": {"color": "rgb(133,20,75)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "6d9e58", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace8 = {
  "x": [4, 4, 3, 3], 
  "y": [0, 5.7, 5.7, 0], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "8b861d", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace9 = {
  "x": [105, 105, 120, 120], 
  "y": [0, 31.5234332069, 31.5234332069, 26.4549357718], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "7adf4f", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace10 = {
  "x": [135, 135, 145, 145], 
  "y": [0, 25.4921745996, 25.4921745996, 0], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "6ed3ca", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace11 = {
  "x": [155, 155, 165, 165], 
  "y": [0, 31.748077587, 31.748077587, 0], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "3b6ff5", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace12 = {
  "x": [140, 140, 160, 160], 
  "y": [25.4921745996, 35.5366510598, 35.5366510598, 31.748077587], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "16f5e5", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace13 = {
  "x": [185, 185, 195, 195], 
  "y": [0, 26.3123066366, 26.3123066366, 0], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "9f7f80", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace14 = {
  "x": [175, 175, 190, 190], 
  "y": [0, 43.0139782528, 43.0139782528, 26.3123066366], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "104f9c", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace15 = {
  "x": [150, 150, 182.5, 182.5], 
  "y": [35.5366510598, 57.5617726333, 57.5617726333, 43.0139782528], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "117b52", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace16 = {
  "x": [112.5, 112.5, 166.25, 166.25], 
  "y": [31.5234332069, 59.1893157854, 59.1893157854, 57.5617726333], 
  "marker": {"color": "rgb(255,220,0)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "5fe6e5", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace17 = {
  "x": [215, 215, 225, 225], 
  "y": [0, 27.1521415541, 27.1521415541, 0], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "4e204f", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace18 = {
  "x": [255, 255, 265, 265], 
  "y": [0, 20.0515918243, 20.0515918243, 0], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "efe1c4", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace19 = {
  "x": [245, 245, 260, 260], 
  "y": [0, 27.5758635995, 27.5758635995, 20.0515918243], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "440a6d", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace20 = {
  "x": [235, 235, 252.5, 252.5], 
  "y": [0, 41.8133183872, 41.8133183872, 27.5758635995], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "9e3e0e", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace21 = {
  "x": [220, 220, 243.75, 243.75], 
  "y": [27.1521415541, 47.4525930617, 47.4525930617, 41.8133183872], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "d2cf06", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace22 = {
  "x": [205, 205, 231.875, 231.875], 
  "y": [0, 50.4718839242, 50.4718839242, 47.4525930617], 
  "marker": {"color": "rgb(40,35,35)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "7e9273", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace23 = {
  "x": [285, 285, 295, 295], 
  "y": [0, 19.2926252199, 19.2926252199, 0], 
  "marker": {"color": "rgb(61,153,112)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "b1f837", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace24 = {
  "x": [275, 275, 290, 290], 
  "y": [0, 26.7385370449, 26.7385370449, 19.2926252199], 
  "marker": {"color": "rgb(61,153,112)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "42e830", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace25 = {
  "x": [315, 315, 325, 325], 
  "y": [0, 20.494871671, 20.494871671, 0], 
  "marker": {"color": "rgb(255,65,54)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "95da30", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace26 = {
  "x": [305, 305, 320, 320], 
  "y": [0, 45.8593636162, 45.8593636162, 20.494871671], 
  "marker": {"color": "rgb(255,65,54)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "ad1b21", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace27 = {
  "x": [282.5, 282.5, 312.5, 312.5], 
  "y": [26.7385370449, 63.1884559196, 63.1884559196, 45.8593636162], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "86b634", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace28 = {
  "x": [218.4375, 218.4375, 297.5, 297.5], 
  "y": [50.4718839242, 66.5790164593, 66.5790164593, 63.1884559196], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "73556d", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace29 = {
  "x": [139.375, 139.375, 257.96875, 257.96875], 
  "y": [59.1893157854, 69.1184493533, 69.1184493533, 66.5790164593], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "0b7e14", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace30 = {
  "x": [90, 90, 198.671875, 198.671875], 
  "y": [32.2392430018, 74.0163384507, 74.0163384507, 69.1184493533], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "28e368", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace31 = {
  "x": [70, 70, 144.3359375, 144.3359375], 
  "y": [35.3429593137, 77.1086406846, 77.1086406846, 74.0163384507], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "dc3433", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace32 = {
  "x": [27.5, 27.5, 107.16796875, 107.16796875], 
  "y": [63.7396717612, 88.2854179415, 88.2854179415, 77.1086406846], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "uid": "312974", 
  "xaxis": "x", 
  "yaxis": "y"
}
data = Data([trace1, trace3, trace4])
#data = Data([trace1, trace2, trace3, trace4, trace5, trace6, trace7, trace8, trace9, trace10, trace11, trace12, trace13, trace14, trace15, trace16, trace17, trace18, trace19, trace20, trace21, trace22, trace23, trace24, trace25, trace26, trace27, trace28, trace29, trace30, trace31, trace32])
layout = {
  "autosize": False, 
  "height": 960, 
  "hovermode": "closest", 
  "margin": {
    "r": 100, 
    "t": 100, 
    "b": 100, 
    "l": 100, 
    "pad": 20
  }, 
  "showlegend": False, 
  "title": "Dendrogram", 
  "width": 1280, 
  "xaxis": {
    "autorange": True, 
    "mirror": "allticks", 
    "range": [0, 325], 
    "rangemode": "tozero", 
    "showgrid": False, 
    "showline": False, 
    "showticklabels": True, 
    "tickmode": "array", 
    "ticks": "", 
    "ticktext": ["<span style=\"fill:#00FF00\">du145r4</span>", "<span style=\"fill:#00FF00\">du145r2</span>", "<span style=\"fill:#00FF00\">du145r3</span>", "<span style=\"fill:#0000FF\">bph1r1</span>", "<span style=\"fill:#0000FF\">bph1r3</span>", "<span style=\"fill:#0000FF\">bph1r4</span>", "<span style=\"fill:#E7298A\">lncapc1r2</span>", "<span style=\"fill:#E7298A\">lncapc1r3</span>", "<span style=\"fill:#008080\">vcapr2</span>", "<span style=\"fill:#008080\">vcapr3</span>", "<span style=\"fill:#008000\">pc3r4</span>", "<span style=\"fill:#008000\">pc3r1</span>", "<span style=\"fill:#008000\">pc3r2</span>", "<span style=\"fill:#00FFFF\">lncapcontrolr3</span>", "<span style=\"fill:#00FFFF\">lncapcontrolr5</span>", "<span style=\"fill:#000080\">lncapc9r3</span>", "<span style=\"fill:#00FFFF\">lncapcontrolr1</span>", "<span style=\"fill:#800080\">rpwe1r3</span>", "<span style=\"fill:#800080\">rpwe1r4</span>", "<span style=\"fill:#800080\">rpwe1r5</span>", "<span style=\"fill:#008080\">vcapr1</span>", "<span style=\"fill:#FFA500\">lncapr1</span>", "<span style=\"fill:#FFA500\">lncapr2</span>", "<span style=\"fill:#E7298A\">lncapc1r1</span>", "<span style=\"fill:#000080\">lncapc9r2</span>", "<span style=\"fill:#FFA500\">lncapr4</span>", "<span style=\"fill:#000080\">lncapc9r1</span>", "<span style=\"fill:#E41A1C\">a22rv1r3</span>", "<span style=\"fill:#E41A1C\">a22rv1r1</span>", "<span style=\"fill:#E41A1C\">a22rv1r2</span>", "<span style=\"fill:#999999\">pwre1r5</span>", "<span style=\"fill:#999999\">pwre1r2</span>", "<span style=\"fill:#999999\">pwre1r4</span>"], 
    "tickvals": [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145], 
    "type": "linear", 
    "zeroline": False
  }, 
  "yaxis": {
    "autorange": True, 
    "mirror": False, 
    "range": [0, 92.9320188858], 
    "rangemode": "tozero", 
    "showgrid": False, 
    "showline": True, 
    "showticklabels": True, 
    "ticks": "outside", 
    "type": "linear", 
    "zeroline": False
  }
}
fig = Figure(data=data, layout=layout)
plot_url = py.plot(fig)