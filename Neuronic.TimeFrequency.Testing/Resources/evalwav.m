function R = evalwav(wname)
[phi, psi, x] = wavefun(wname, 10);
csvwrite(['wavefun_', wname, '.csv'], [psi; x]);
R = [x(1), x(end), size(psi, 2)]
end