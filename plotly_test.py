import plotly as pl
from plotly.graph_objs import Scatter, Layout

def plotly_test():
	pl.offline.plot({
    "data": [Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])],
    "layout": Layout(title="hello world")
})
