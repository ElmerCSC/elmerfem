sqrt=math.sqrt

function normalload(x,t)
 radius = sqrt(tx[0]*tx[0])
 loadperunit =  -100.0*9.81*917.0
 switchtime = 20.0*365.25*24.0*3600.0
 if (t < switchtime) then
   if (radius <= 3.00e04) then
     return loadperunit
   else
     return 0.0
   end
 else
   return 0.0
 end
end 
