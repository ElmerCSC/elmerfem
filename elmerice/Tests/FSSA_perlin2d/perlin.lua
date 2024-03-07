-- constants
-- #########
betamin=0.01
betamax=1
dbeta=10.0
delta=2.0
sigma=delta/4.59511985013458992685
alpha = 1.0
beta = 3.0
LX=8000.0
LY=2500.0

-- the vertical offset to the base topography
-- ##########################################
function parbl(Y)
  aux= LY - math.sqrt(LY*LY - 0.25*Y*Y)
  return aux
end

-- this is our accumulation rate
-- #############################
function accum(X)
  aux=alpha-beta*X/LX
  if (aux > 0.0) then
    return aux
  else
    return 0.0
  end  
end

-- this is our accumulation rate in 3D
-- ###################################
function accum3D(X,Y)
  fade=1.0 - ((Y-2500.0)/LY)^2.0
  return accum(X)*fade
end

-- linear sliding coefficient
-- ##########################
function lslidingcoeff(X,H,D)
 height = H-D
 fact = 1.0/(1.0 + math.exp((height - hmin)/(0.5*sigma)))
 return (betamax - betamin)/(1.0 + math.exp((X - 0.5*LX)/sigma)) + betamin + fact*betamax
end


-- sliding with fixed target velocity ub agnostic of driving stress
-- ################################################################
function slidingcoeffv(X,H,D)
  return lslidingcoeff(X,H,D)/(ub^(mw - 1.0))
end

-- sliding that tries to match the linear sliding speed for given driving stress (2nd argument)
-- ############################################################################################
function slidingcoeff(X,Tau,H,D)
  ub=Tau/lslidingcoeff(X,H,D)
  return lslidingcoeff(X,H,D)/(ub^(mw - 1.0))
end


 