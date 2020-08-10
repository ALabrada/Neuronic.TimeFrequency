function R = evalwav(wname)
Y = wavefun(wname, 10);
psi = Y(2);
x = Y(end);
csvwrite(['wavefun_', wname, '.csv'], [psi; x]);
R = [x(1), x(end), size(psi, 2)]
end