function dXdt = rate_func02(t,X)
dXdt = [0,-1;1,0]*X;
end
function X = solution02(t)
X = [cos(t);sin(t)];
end