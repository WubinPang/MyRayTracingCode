% Optical Imaging and Spectroscopy

% David J. Brady
% Duke University
% www.opticalimaging.org
%
% HW 4.4 The Fraunhofer Diffraction
% Adapted by Wubin Pang

%
clear all;
figure(1);set(gcf,'color','white');
X=1000; %microns
lambda=0.633; %micron
N=1024; % samples per axis
FL=1000; %micron focal length of lens
xrange=linspace(-X,X,N);
fxrange=(-N*FL*lambda/(4*X)):(FL*lambda/(2*X)):(N/2-1)*FL*lambda/(2*X);
[x,y]=meshgrid(xrange);
[fx,fy]=meshgrid(fxrange);

pinh=circle(x,y,0,0,50); %% diameter of pinhole is 50 micron
subplot(2,2,1);imagesc(xrange,xrange,pinh);colormap 'gray';axis 'square';title('Pinhole'); xlabel('x');ylabel('y');
fraunhoferT=fftshift(fft2(pinh));
subplot(2,2,2);imagesc(fxrange,fxrange,abs(fraunhoferT));colormap 'gray';axis 'square';title('Fraunhofer Diffraction of Pinhole');xlabel('fx');ylabel('fy'); zoom(4);


X=50; %microns
lambda=0.633; %micron
N=1024; % samples per axis
FL=1000; %micron focal length of lens
xrange=linspace(-X,X,N);
fxrange=(-N*FL*lambda/(4*X)):(FL*lambda/(2*X)):(N/2-1)*FL*lambda/(2*X);
[x,y]=meshgrid(xrange);
[fx,fy]=meshgrid(fxrange);
urange=(-N/(4*X)):(1/(2*X)):(N/2-1)/(2*X);
[u,v]=meshgrid(urange);

hair=1-rect(x,y,0,0,0,50); %% diameter of hair is 100 micron
subplot(2,2,3);imagesc(xrange,xrange,hair);colormap 'gray';axis 'square';title('hair strand'); xlabel('x');ylabel('y');
fraunhoferT=fftshift(fft(hair));
fresnelT=exp(1i*pi*FL*lambda*(u.^2));
ft=ifft(fftshift(fresnelT).*fft(fraunhoferT'));
fraunhoferT=ft';

subplot(2,2,4);imagesc(xrange,fxrange,abs(fraunhoferT));colormap 'gray';axis 'square';title('Fraunhofer Diffraction of hair strand');xlabel('fx');ylabel('fy'); zoom(40);
