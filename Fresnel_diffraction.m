%
% Optical Imaging and Spectroscopy
%
% David J. Brady
% Duke University
% www.opticalimaging.org
%
% Figure 4.4
%
%
% Fresnel diffraction
%
%
figure(1);set(gcf,'color','white');
X=1000; %microns
lambda=1; %micron
N=1024; % samples per axis
xrange=linspace(-X,X,N);
urange=(-N/(4*X)):(1/(2*X)):(N/2-1)/(2*X);
[x,y]=meshgrid(xrange);
[u,v]=meshgrid(urange);
% input pattern is modulated by 1 mm gaussian chirp pattern varying from
% 100 micron period to 10 micron period
%
f=exp(-pi*(x.^2+y.^2)/250^2).*(1+cos(pi*.05*x));
subplot(2,2,1);imagesc(xrange,xrange,f);colormap 'gray';axis 'square';title('f(x,y)');zoom(2);
% d=1000; % microns
% fresnelT=exp(i*pi*d*lambda*(u.^2+v.^2));
% ft=ifft2(fftshift(fresnelT).*fft2(f));
% subplot(2,2,2);imagesc(xrange,xrange,abs(ft));colormap 'gray';axis 'square';title('d=1 mm');zoom(2);
d=2000; % microns
fresnelT=exp(i*pi*d*lambda*(u.^2+v.^2));
ft=ifft2(fftshift(fresnelT).*fft2(f));
subplot(2,2,3);imagesc(xrange,xrange,abs(ft));colormap 'gray';axis 'square';title('d=2 mm');zoom(2);
d=10000; % microns
fresnelT=exp(i*pi*d*lambda*(u.^2+v.^2));
ft=ifft2(fftshift(fresnelT).*fft2(f));
subplot(2,2,4);imagesc(xrange,xrange,abs(ft));colormap 'gray';axis 'square';title('d=10 mm');zoom(2);

