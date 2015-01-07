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

'enable print tottrnsp.gx'

rm = 0.0
tm = 0.5
dx = 3.2
dy = 2.5


var.1 = 'vta'
var.2 = 'vsa'
var.3 = 'vtp'
var.4 = 'vsp'



trname.1 = 'meridional heat transport Atlantic' 
trname.2 = 'meridional salt transport Atlantic' 
trname.3 = 'meridional heat transport Pacific' 
trname.4 = 'meridional salt transport Pacific' 


'set t 1 'nt''
'set lat 1'
'set lon 1'

p = 1
it = 1

while ( p <= 1)

  'clear'
    
    fx = 1
    while (fx <= 2)
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


        'draw title  'trname.it' '

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


