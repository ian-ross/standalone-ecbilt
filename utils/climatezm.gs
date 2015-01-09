function main(args)

fname=subwrd(args,1)
trun=subwrd(args,2)

'open 'fname'grdstl.ctl'

'clear'
'set lev 1000'

'enable print 'fname'zm.gx'

'set grads off'


************************************ T2m ***********************************
var=t
des=T2m
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ T650 ***********************************
var=t
des=T650
il1=650
il2=650

plot(var,des,il1,il2,trun)
************************************ T350 ***********************************
var=t
des=T350
il1=350
il2=350

plot(var,des,il1,il2,trun)
************************************ teta ***********************************
var='((t(lev=350)+273.15)/0.741-(t(lev=650)+273.15)/0.884)*0.5'
des=teta
il1=350
il2=350

plot(var,des,il1,il2,trun)

************************************ T(p) ***********************************
var=t
des='temperature'
il1=1000
il2=350

plot(var,des,il1,il2,trun)

************************************ W(p) ***********************************
var=omega
des='omega'
il1=1000
il2=350

plot(var,des,il1,il2,trun)

************************************ prec ***********************************

var=pp
des=prec.
il1=1000
il2=1000

plot(var,des,il1,il2,trun)

************************************ evap ***********************************

var=evap
des=evap.
il1=1000
il2=1000

plot(var,des,il1,il2,trun)

************************************ eminp ***********************************

var='evap-pp'
des=eminp
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ pwat ***********************************

var=q%'*1000'
des=prec%'.water(kg/m2)'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ olr ***********************************

var=ttr
des='OLR (W/m2)'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ slwr ***********************************

var=str
des='SLWR (W/m2)'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ tsr ***********************************

var=tsr
des='TSR (W/m2)'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ ssr ***********************************

var=ssr
des='SSR (W/m2)'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
************************************ alb ***********************************

var=alb
des='Surface albedo'
il1=1000
il2=1000

plot(var,des,il1,il2,trun)
'close 1'

'open 'fname'grdsul.ctl'
'set x 1'

************************************ u800 ***********************************
var=u
des='zonal wind (m/s)'
il1=800
il2=200

plot(var,des,il1,il2,trun)

'disable print'
'!gxps.sc 'fname'zm.gx 'fname'zm.ps'
'!psfix 'fname'zm.ps 8'
'!mv 'fname'zmfix.ps 'fname'zm.ps'
'!rm in.ps '
'quit'
return

***************************************************************************

function plot(var,des,il1,il2,trun)

'clear'
iyear=2
'set x 1 64'
'define win=0'
'set lev 'il1' 'il2

while (iyear < 9)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win=win+'var''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win=win/'iyear''

dy=5.2
ds=0.5
yl=0
x1=0
x2=8.5
y2=yl+2*dy
y1=y2-dy-0.5*ds

'set vpage 'x1' 'x2' 'y1' 'y2
'set grads off'
'set x 1'
'define winz=ave(win,x=1,x=64)'
'd winz'
'draw title 'des' djf 'trun''

iyear=2
'set x 1 64'
'define sum=0'

while (iyear < 9)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum=sum+'var''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum=sum/'iyear''

y2=yl+dy+0.5*ds
y1=yl
'set vpage 'x1' 'x2' 'y1' 'y2
'set grads off'
'set x 1'
'define sumz=ave(sum,x=1,x=64)'
'd sumz'
'draw title 'des' jja 'trun''
'print'
'set x 1 64'
*pull return
return
  
