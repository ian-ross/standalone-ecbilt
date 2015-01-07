function main(args1)
**
**
** this script makes plots of timeseries of basin mean fluxes
** between ocean and atmosphere
** be sure to open the ctl file of your mean basin data prior to
** running this script
**
**

rec   = args1

nt      = subwrd(rec,1)

'enable print totfx.gx'

rm = 0.0
tm = 0.5
dx = 3.2
dy = 2.5


var.1 = 'rlat'
var.2 = 'sen'
var.3 = 'dlrad'
var.4 = 'ulrad'
var.5 = 'solar'
var.6 = 'hefx'
var.7 = 'pp'
var.8 = 'runo'
var.9 = 'brine'
var.10 = 'safx'
var.11 = 'icev'
var.12 = 'icea'



flux.1 = 1
flux.2 = 2
flux.3 = 3

flux.4 = 4
flux.5 = 5
flux.6 = 6

flux.7 = 7
flux.8 = 8
flux.9 = 9

flux.10 = 10
flux.11 = 11
flux.12 = 12


fxname.1 = 'latent heat flux' 
fxname.2 = 'sensible heatflux' 
fxname.3 = 'downward longwave' 
fxname.4 = 'upward longwave' 
fxname.5 = 'solar radiation' 
fxname.6 = 'total heatflux' 
fxname.7 = 'precipitation' 
fxname.8 = 'runoff' 
fxname.9 = 'brine rejection' 
fxname.10 = 'total salt flux' 
fxname.11 = 'ice volume' 
fxname.12 = 'ice area' 


'set t 1 'nt''
'set lat 1'
'set lon 1'

p = 1
it = 1

while ( p <= 2)

  'clear'
    
    fx = 1
    while (fx <= 3)
     fy = 1
     while (fy <= 2)
      
**      Set subpage's coordinate.

        x1 = rm + (fx-1)*dx
        x2 = x1 + dx
        y1 = 8.5 - tm - fy*dy
        y2 = y1 + dy

        'set vpage 'x1' 'x2' 'y1' 'y2''
        'set grads off'
        'set cmark 0'
        'd 'var.it''

**      Set the size of character string.
        'set strsiz 0.1'

**      Draw title for the plot


        'draw title  'fxname.it' '

         it = it +1 

         fy = fy + 1
        endwhile
      fx = fx + 1
     endwhile
    p = p + 1
   'print'
  
endwhile

** Make hardcopy.
  'disable print'
say 'Plotting is finished'


