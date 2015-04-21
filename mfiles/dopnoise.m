function [y,iflaw]=dopnoise(N,Fs,f0,d,v,t0,c);
%DOPNOISE Complex noisy doppler signal.
%	[Y,IFLAW]=DOPNOISE(N,FS,F0,D,V,T0,C) generates a complex 
%	noisy doppler signal, normalized so as to be of unit energy. 
%
%	 N  : number of points.  
%	 FS : sampling frequency (in Hertz).  
%	 F0 : target frequency   (in Hertz).  
%	 D  : distance from the line to the observer (in meters).  
%	 V  : target velocity    (in m/s).  
%	 T0 : time center                  (default : N/2).  
%	 C  : wave velocity      (in m/s)  (default : 340).
%	 Y  : Output signal.
%	 IFLAW : Model used as instantaneous frequency law.
%
%	Example :
%	 [z,iflaw]=dopnoise(500,200,60,10,70,128);
%	 subplot(211); plot(real(z)); 
%	 subplot(212); plot(iflaw); hold; ifl=instfreq(z,11:478,10);
%	 plot(ifl,'g'); hold; sum(abs(z).^2)
%	 axis([0 468 -0.01 0.51])

%	F. Auger, July 94, August 95.
%	Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

if (nargin <= 4),
 error ( 'At least 5 parameters are required' ); 
elseif (nargin == 5),
 t0=N/2; c=340.0;
elseif (nargin == 6),
 c=340.0;
end;

if (N <= 0),
 error ('The signal length N must be strictly positive' );
elseif (Fs<0),
 error ('The sampling frequency FS must be positive');
elseif (t0<1)|(t0>N),
 error ('T0 must be between 1 and N');
elseif (f0<0)|(f0>Fs/2),
 error ('F0 must be between 0 and FS/2');
elseif (v<0),
 error ('V must be positive');
else
 r=0.9; rr=r*r; r2=r*2; vv=v*v;
 x=randn(2*N,1);
 tmt0=((1:2*N)'-t0-N)/Fs ; dist=sqrt(d^2+(v*tmt0).^2);
 iflaw=(1-vv*tmt0./dist/c)*f0/Fs; 
 y=zeros(2*N-2,1);
 for t=3:2*N,
  y(t)=x(t)-rr*(x(t-2)+y(t-2))+r2*cos(2.0*pi*iflaw(t))*y(t-1);
 end; 
 y=hilbert(y(N+1:2*N) ./ sqrt(dist(N+1:2*N))); 
 y=y/sqrt(sum(abs(y).^2)); iflaw=iflaw(N+1:2*N);
end;
