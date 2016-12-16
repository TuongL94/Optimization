[t1 y1] = data1;
[t2 y2] = data2;
tol = 10^-3;
xstart = [1;2;3;4];

xmin = gaussnewton(@phi2,t1,y1,xstart,tol,1,0,0);
