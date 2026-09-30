
```{note}
AbiPy provides two different APIs to produce figures, with either matplotlib or plotly.
In this tutorial, we use the plotly API as much as possible, although it should be noted
that not all the plotting methods have been ported to plotly yet.

AbiPy uses a relatively simple rule to differentiate between the two plotting libraries:
if an object provides an `obj.plot` method producing a matplotlib plot, the corresponding
native plotly version (if available) is named `obj.plotly`.
Note that plotly requires a web browser, hence the matplotlib version is still valuable if you need to
visualize results on machines on which only an X server is available.

In AbiPy versions greater than 0.9, you can try to convert the matplotlib figures produced by `plot` methods
into Plotly figures using the optional argument `plotly=True`, although the result is not always guaranteed.
```
