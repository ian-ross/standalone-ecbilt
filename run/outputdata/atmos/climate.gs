function main(args)


fname=subwrd(args,1)
trun=subwrd(args,2)

'open 'fname'grdstl.ctl'

'clear'
'set lev 1000'

'enable print 'fname'.gx'

'set lon -180 180'
'set grads off'


************************************ T2m ***********************************
var=t
des=T2m
ilev=1000
cint=5

plot(var,des,ilev,trun,cint)

************************************ T650 ***********************************
var=t
des=T650
ilev=650
cint=5

plot(var,des,ilev,trun,cint)

************************************ T350 ***********************************
var=t
des=T350
ilev=350
cint=5

plot(var,des,ilev,trun,cint)


************************************ teta ***********************************
var='((t(lev=350)+273.15)/0.741-(t(lev=650)+273.15)/0.884)*0.5'
des=teta
ilev=350
cint=2

plot(var,des,ilev,trun,cint)


************************************ prec ***********************************

var=pp
des='prec. [cm/yr]'
ilev=1000
cint=50

plot(var,des,ilev,trun,cint)

************************************ eminp **********************************

var='evap-pp'
des='evap. [cm/yr]'
ilev=1000
cint=50

plot(var,des,ilev,trun,cint)

************************************ pwat ***********************************

var=q%'*1000'
des=prec%'.water(kg/m2)'
ilev=1000
cint=5

plot(var,des,ilev,trun,cint)

************************************ olr ***********************************

var=ttr
des='OLR (W/m2)'
ilev=1000
cint=25

plot(var,des,ilev,trun,cint)

************************************ slwr ***********************************

var=str
des='SLWR (W/m2)'
ilev=1000
cint=10

plot(var,des,ilev,trun,cint)

************************************ tsr ***********************************

var=tsr
des='TSR (W/m2)'
ilev=1000
cint=25

plot(var,des,ilev,trun,cint)

************************************ ssr ***********************************

var=ssr
des='SSR (W/m2)'
ilev=1000
cint=25

plot(var,des,ilev,trun,cint)

************************************ tau ***********************************
var1=ustress
var2=vstress
des='windstress'
ilev=1000
cint=5

*pltvec(var1,var2,des,ilev,trun,cint)

'close 1'

'open 'fname'grdsul.ctl'
'set lon -180 180'

************************************ u800 ***********************************
var=u
des='800hPa u(m/s)'
ilev=800
cint=3

plot(var,des,ilev,trun,cint)

************************************ u500 ***********************************
var=u
des='500hPa u(m/s)'
ilev=500
cint=5

plot(var,des,ilev,trun,cint)

************************************ u200 ***********************************
var=u
des='200hPa u(m/s)'
ilev=200
cint=5

plot(var,des,ilev,trun,cint)

************************************ psi 500*********************************
var='psi*6187+550'
des='500hPa Geopotential (dm)'
ilev=500
cint=10

'set mproj nps'
'set lat 20 90'

plot(var,des,ilev,trun,cint)
'disable print'
'!gxps -i 'fname'.gx -o 'fname'.ps'
'!psfix 'fname'.ps 4'
'!mv 'fname'fix.ps 'fname'.ps'
'!rm in.ps '
'close 1'
'reinit'
'quit'
return

***************************************************************************

function plot(var,des,ilev,trun,cint)

'clear'
iyear=2
'define win=0'
'set lev 'ilev

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
'set cint 'cint
'd win'
'draw title 'des' djf 'trun''

iyear=2
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
'set cint 'cint
'd sum'
'draw title 'des' jja 'trun''
'print'
*pull return
return
  
***************************************************************************

function pltvec(var1,var2,des,ilev,trun,cint)

'clear'
iyear=2
'define win1=0'
'define win2=0'
'set lev 'ilev

while (iyear < 9)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win1=win1+'var1''
  'define win2=win2+'var2''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win1=win1/'iyear''
'define win2=win2/'iyear''

dy=5.2
ds=0.5
yl=0
x1=0
x2=8.5
y2=yl+2*dy
y1=y2-dy-0.5*ds

'set vpage 'x1' 'x2' 'y1' 'y2
'set grads off'
'set cint 'cint
'd win1;win2'
'draw title 'des' djf 'trun''

iyear=2
'define sum1=0'
'define sum2=0'

while (iyear < 9)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum1=sum1+'var1''
  'define sum2=sum2+'var2''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum1=sum1/'iyear''
'define sum2=sum2/'iyear''

y2=yl+dy+0.5*ds
y1=yl
'set vpage 'x1' 'x2' 'y1' 'y2
'set grads off'
'set cint 'cint
'd sum1;sum2'
'draw title 'des' jja 'trun''
'print'
*pull return
return
  
