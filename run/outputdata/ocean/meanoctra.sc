function main(args1)
**
**
** this script makes plots of timeseries of basin mean ocean tracers
** be sure to open the ctl file of your mean basin data prior to
** running this script
**
**

rec   = args1

var     = subwrd(rec,1)
nt      = subwrd(rec,2)

'enable print 'var'.gx'

rm = 0.5
tm = 0.5
dx = 3.2
dy = 2.5

bas.1.1 = 1
bas.1.2 = 3
bas.1.3 = 5

bas.2.1 = 2
bas.2.2 = 4
bas.2.3 = 6

bas.3.1 = 7
bas.3.2 = 8
bas.3.3 = 9

bname.1 = 'North Atlantic' 
bname.2 = 'North Pacific' 
bname.3 = 'Midlat Atlantic' 
bname.4 = 'Midlat Pacific' 
bname.5 = 'Tropical Atlantic' 
bname.6 = 'Tropical Pacific' 
bname.7 = 'Indian ocean' 
bname.8 = 'Midlat Southern' 
bname.9 = 'Antartic circ polar' 

lev.1 = 1
lev.2 = 6
lev.3 = 12

'set t 1 'nt''

p = 1

while ( p <= 3)
  'clear'
  b = 1
  while (b <= 3)
    l = 1
    while (l <= 3)
      
**    Set ocean basin

      'set lat 'bas.p.b''

**    Set vertical level

      'set lon 'lev.l''

**    Set subpage's coordinate.

      x1 = rm + (b-1)*dx
      x2 = x1 + dx
      y1 = 8.5 - tm - l*dy
      y2 = y1 + dy

      'set vpage 'x1' 'x2' 'y1' 'y2''
      'set grads off'

      'set cmark 0'
      'd 'var''

**    Set the size of character string.
      'set strsiz 0.1'

**    Draw title for the plot

      i = bas.p.b

      'draw title  'var' of 'bname.i' at level 'lev.l''

      l = l + 1
    endwhile
    b = b + 1
  endwhile
  p = p + 1
  'print'
endwhile

** Make hardcopy.
  'disable print'
say 'Plotting is finished'


