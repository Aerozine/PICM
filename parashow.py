from paraview.simple import *

reader = OpenDataFile('results/simulation.pvd')
Show(reader)
Render()
# rescale time 
#GetDisplayProperties().RescaleTransferFunctionToDataRangeOverTime()
GetAnimationScene().Play()
