
function deriv = stocheck(t, x, param)
sGEprime = param.sGEprime;
Sto = x(1);
dSto_dt = ppval(sGEprime, t);
deriv = dSto_dt;
end