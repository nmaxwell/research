
from setup import *


def render( data, fname ):
    write_png( data, fname, center=mean(data), major_scale=std(data)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )




grid = grid2d( a1 = -500, b1 = 3500, a2=-300, b2=1200, n1=1024, n2=384)
u = grid.evaluate( lambda x,y: exp(1.0E-6*(-(x-1500)**2-(y-800)**2)))
count = 1   
fname = "/tmp/research/" + str(count)  
scipy.io.savemat( fname, {'m':u } )
data = scipy.io.loadmat( fname  )
v = data['m']
w = v.copy()
render(u, "/workspace/output/scratch/" + "u" + ".png" )
render(v, "/workspace/output/scratch/" + "v" + ".png" )
render(w, "/workspace/output/scratch/" + "w" + ".png" )





