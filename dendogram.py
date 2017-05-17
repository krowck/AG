import plotly.plotly as py
from plotly.graph_objs import *

py.sign_in('krowck', '642031qwe')
trace1 = {
  "x": [5.0, 5.0, 15.0, 15.0], 
  "y": [0.0, 0.870574522944, 0.870574522944, 0.0], 
  "marker": {"color": "rgb(61,153,112)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace2 = {
  "x": [25.0, 25.0, 35.0, 35.0], 
  "y": [0.0, 1.25980157168, 1.25980157168, 0.0], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace3 = {
  "x": [45.0, 45.0, 55.0, 55.0], 
  "y": [0.0, 1.08917399896, 1.08917399896, 0.0], 
  "marker": {"color": "rgb(255,65,54)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace4 = {
  "x": [65.0, 65.0, 75.0, 75.0], 
  "y": [0.0, 1.21799014774, 1.21799014774, 0.0], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace5 = {
  "x": [50.0, 50.0, 70.0, 70.0], 
  "y": [1.08917399896, 1.34376337203, 1.34376337203, 1.21799014774], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace6 = {
  "x": [30.0, 30.0, 60.0, 60.0], 
  "y": [1.25980157168, 1.4436758639, 1.4436758639, 1.34376337203], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
trace7 = {
  "x": [10.0, 10.0, 45.0, 45.0], 
  "y": [0.870574522944, 1.58094908204, 1.58094908204, 1.4436758639], 
  "marker": {"color": "rgb(0,116,217)"}, 
  "mode": "lines", 
  "type": "scatter", 
  "xaxis": "x", 
  "yaxis": "y"
}
data = Data([trace1, trace2, trace3, trace4, trace5, trace6, trace7])
layout = {
  "autosize": False, 
  "height": 800, 
  "hovermode": "closest", 
  "showlegend": False, 
  "title": "Basket Dendrogram 2013", 
  "width": 800, 
  "xaxis": {
    "mirror": "allticks", 
    "rangemode": "tozero", 
    "showgrid": False, 
    "showline": True, 
    "showticklabels": True, 
    "tickmode": "array", 
    "ticks": "outside", 
    "ticktext": ["Analogfotos", "Sonstige", "Fotos", "Wandbilder", "Fotobuch", "Fotokalender", "Fotogeschenke", "Grußkarten"], 
    "tickvals": [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0], 
    "type": "linear", 
    "zeroline": False
  }, 
  "yaxis": {
    "mirror": "allticks", 
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