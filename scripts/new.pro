set_plot,'ps'
!p.font=0
for i = 1,12 do begin
n=i
num=strtrim(n,2)
file = 'ryujones'+num+'.txt'

print, file
rESTORE, 'myPlotTemplate.sav'  

plot_ascii= read_ascii(file , template=myplottemplate)

density=plot_Ascii.field01
velx=plot_Ascii.field02
vely=plot_Ascii.field03
velz=plot_Ascii.field04
energy=plot_Ascii.field05
bx=plot_Ascii.field06
by=plot_Ascii.field07
bz=plot_Ascii.field08
pressure=plot_Ascii.field09

;device, decomposed=0, retain=2
;Window, 0, Title='Density Plot',xsize=800,ysize=800


;rescol





;white = GetColor('White', 1)
;black = GetColor('Black', 2)
;LoadCT, 33, NColors=64

   !Y.OMARGIN=[2, 10]
!p.multi=[0,3,3]
plot, density, ytitle='density' ;,psym=4
plot, pressure, ytitle='Pressure' ;,psym=4
plot, energy, ytitle='Energy' ;,psym=4
plot, velx, ytitle='Velx' ;,psym=4
plot, vely, ytitle='Vely' ;,psym=4
plot, velz, ytitle='Velz' ;,psym=4
plot, Bx, ytitle='Bx' ;,psym=4
plot, By, ytitle='By' ;,psym=4
plot, Bz, ytitle='Bz' ;,psym=4
   XYOUTS, 0.5, 0.93, ALIGN=0.5, CHARSIZE=1.0, /NORMAL, $
	      'Ryu-Jones test suite #'+ num+' '+systime()


!p.font=-1

endfor

device,/close
END
