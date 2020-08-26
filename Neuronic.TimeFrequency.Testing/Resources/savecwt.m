function savecwt(S,wname)
scales = 1:10;
c = cwt(S,scales,wname);
R = [S; c];
save(['cwt_', wname, '.mat'], 'R');
end