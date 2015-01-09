function main(args)

fname=subwrd(args,1)
nm=subwrd(args,2)

'open 'fname'grdstl.ctl'

say "Global Mean Budgets "fname
say "-------------------------------"
say " "

fname1=fname
fname=fname'GB'

code=write(''fname'.dat','Global Mean Budgets 'fname1'')
code=write(''fname'.dat',"-------------------------------")
code=write(''fname'.dat'," ")

************************************ olr ***********************************

var=ttr
des='OLR (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ slwr ***********************************

var=str
des='SLWR (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ sulwr ***********************************
var=dumt1
des='SULWR (W/m2)'
ilev=350

*plot(var,des,ilev,fname,nm)
************************************ sdlwr ***********************************

var=dumt1
des='SDLWR (W/m2)'
ilev=650

*plot(var,des,ilev,fname,nm)
************************************ tsr ***********************************

var=tsr
des='TSR (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ ssr ***********************************

var=ssr
des='SSR (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)


************************************ shf ***********************************

var=shf
des='SHF (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ lhf ***********************************

var=lhf
des='LHF (W/m2)'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ palb **********************************

var=palb
des='planetary albedo'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ alb **********************************

var=alb
des='surface albedo'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ prec **********************************

var=pp
des='precipitation'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ evap **********************************

var=evap
des='evaporation'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ t10 **********************************

var=t
des='t10'
ilev=1000

plot(var,des,ilev,fname,nm)

************************************ pwat **********************************

var='q*1000'
des='prec.water [kg/m^2]'
ilev=1000

plot(var,des,ilev,fname,nm)

******************************net surface budget****************************

des='Net SFC (W/m2)'

plotnet1(des,fname,nm)


******************************net toa budget****************************

des='Net TOA (W/m2)'

plotnet2(des,fname,nm)

'close 1'
'quit'
return


***************************************************************************


function plot(var,des,ilev,fname,nm)

iyear=2
'define win=0'
'set lev 'ilev

while (iyear < nm)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win=win+'var''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win=win/'iyear''
'define wina=aave(win,x=1,x=64,y=1,y=32)'

'd wina'
res=des' DJF    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

iyear=2
'define sum=0'

while (iyear < nm)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum=sum+'var''
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum=sum/'iyear''
'define suma=aave(sum,x=1,x=64,y=1,y=32)'

'd suma'
res=des' JJA    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')
tm=(nm-1)*4
'define an=0'
'define ana=0'
'define an=ave('var',t=1,t='tm')'
'define ana=aave(an,x=1,x=64,y=1,y=32)'

'd ana'
res=des' annual 'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

say "************"
code=write(''fname'.dat',"************")

return

***************************************************************************


function plotulw(des,fname,nm)


iyear=2
'define win=0'
'set lev '1000

while (iyear < nm)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win=win+dumt1(lev=350)+str'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win=win/'iyear''
'define wina=aave(win,x=1,x=64,y=1,y=32)'

'd wina'
res=des' DJF    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

iyear=2
'define sum=0'

while (iyear < nm)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum=sum+dumt1(lev=350)+str'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum=sum/'iyear''
'define suma=aave(sum,x=1,x=64,y=1,y=32)'

'd suma'
res=des' JJA    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

'define an=0'
'define ana=0'
'define an=ave(dumt1(lev=350)+str,t=5,t='(nm-1)*12')'
'define ana=aave(an,x=1,x=64,y=1,y=32)'

'd ana'
res=des' annual 'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

say "************"
code=write(''fname'.dat',"************")

return
***************************************************************************


function plotnet1(des,fname,nm)


iyear=2
'define win=0'
'set lev '1000

while (iyear < nm)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win=win+str+shf+lhf-ssr'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win=win/'iyear''
'define wina=aave(win,x=1,x=64,y=1,y=32)'

'd wina'
res=des' DJF    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

iyear=2
'define sum=0'

while (iyear < nm)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum=sum+str+shf+lhf-ssr'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum=sum/'iyear''
'define suma=aave(sum,x=1,x=64,y=1,y=32)'

'd suma'
res=des' JJA    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')
tm=(nm-1)*4
'define an=0'
'define ana=0'
'define an=ave(str+shf+lhf-ssr,t=5,t='tm')'
'define ana=aave(an,x=1,x=64,y=1,y=32)'

'd ana'
res=des' annual 'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

say "************"
code=write(''fname'.dat',"************")

return
***************************************************************************


function plotnet2(des,fname,nm)


iyear=2
'define win=0'
'set lev '1000

while (iyear < nm)
  t1=(iyear-1)*4+1
  'set t 't1
  'define win=win+ttr-tsr'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define win=win/'iyear''
'define wina=aave(win,x=1,x=64,y=1,y=32)'

'd wina'
res=des' DJF    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

iyear=2
'define sum=0'

while (iyear < nm)
  t1=(iyear-1)*4+3
  'set t 't1
  'define sum=sum+ttr-tsr'
  iyear=iyear+1
endwhile

iyear=iyear-2  
  
'define sum=sum/'iyear''
'define suma=aave(sum,x=1,x=64,y=1,y=32)'

'd suma'
res=des' JJA    'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')
tm=(nm-1)*4
'define an=0'
'define ana=0'
'define an=ave(ttr-tsr,t=5,t='tm')'
'define ana=aave(an,x=1,x=64,y=1,y=32)'

'd ana'
res=des' annual 'subwrd(result,4)
say res
code=write(''fname'.dat',''res'')

say "************"
code=write(''fname'.dat',"************")

return
