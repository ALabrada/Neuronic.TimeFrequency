function mycwt(y,T,wname)

[wtype,fname,family,bounds] =  ...
    wavemngr('fields',wname,'type','file','fn','bounds');      

a = 1:10;
n = length(y);
t = ((0:(n-1)) - n/2)*T;
oms = 2*pi/T;
rows=length(a);
yhat=fft(y);
%Matrix-Initialization
matrix=zeros(rows,n);
%Loop for increasing scale factors
for i=1:rows,
    %feval(fname,lb,ub,np,wname)
    psi_scale=conj(feval(fname,-t(end)/a(i),-t(1)/a(i),n,wname));
    psi_scale_hat=fft(fliplr(psi_scale));
    %Time translation such that minimal time=0
    trans=exp((-1i*t(1)*(0:(n-1))*oms/n));
    %Fourier transform of wavelet transform;
    conv_hat=((yhat.*psi_scale_hat).*trans)/sqrt(a(i));
    matrix(i,:)=ifft(conv_hat);
end

csvwrite(['mycwt_', wname, '.csv'], [t; y; matrix]);
end