function main(args1)
**
**
**  This script plots lat-height plot of the var (variable) field.
**  it plots 3 plots on one page, each for one time step.
**  It can make multiople pages depending on the number that you choose. 
**  How to choose the number see below. 
**
**  To run this script in grads, type:
**  'run pxlev.sc2 ctl-file variablename level1 level2 n'.
**  For example: 'run vplot.sc WIND.ctl u 1 3 4'. The script 
**  will use WIND.ctl as control file for grads and plot lat-height plot
**  of u field from level 1 to level 3, it will produce 4 pages plots, on each 
**  page 3 time step's field are plotted. So in this case in total 12 time steps
**  data field are plotted.
 
**  This script also can be run in batch mode.
**  The hardcopy it produces named var.plotvn. 

**  written by Xueli Wang,  March 1995. Last revision 30 May, 1995.
**  !!!!!!!!
**  Warning !!!: to run this script, no control file need to be opened before 
**  hands. If you have control file already opened in this grads session then 
**  be sure that the control file you type in will be the default file.


rec   = args1

var     = subwrd(rec,1)
nt      = subwrd(rec,2)

'enable print 'var'.gx'

rm = 0.5
tm = 0.5
dx = 3.2
dy = 2.5

bas.1.1 = 1
bas.1.2 = 2
bas.1.3 = 3

bas.2.1 = 4
bas.2.2 = 5
bas.2.3 = 6

bas.3.1 = 7
bas.3.2 = 8

bname.1 = 'Arctic ocean' 
bname.2 = 'Hudson Bay' 
bname.3 = 'Bering sea' 
bname.4 = 'Ochotsk sea' 
bname.5 = 'Mediteranean' 
bname.6 = 'South Chinese sea' 
bname.7 = 'Timor sea' 
bname.8 = 'Gulf of Mexico' 

lev.1 = 1
lev.2 = 3
lev.3 = 6

pag.1 = 3
pag.2 = 3
pag.3 = 2

'set t 1 'nt''

p = 1

while ( p <= 3)
  'clear'
  b = 1
  while (b <= pag.p)
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


