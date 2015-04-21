function x=anasing(N,t0,h) ;
%ANASING Lipschitz singularity.
%	X=ANASING(N,T0,H) generates the N-points Lipschitz singularity 
%	centered around T=T0 : X(T) = |T-T0|^H. 
%
%	N  : number of points in time
%	T0 : time localization of the singularity  (default : N/2)
%	H  : strenght of the Lipschitz singularity (positive or
%	     negative)                             (default : 0.0)
%	X  : the time row vector containing the signal samples
%
%	Example :  
%	 x=anasing(128); plot(real(x));
%
%	See also anastep, anapulse, anabpsk, doppler.

%	P. Goncalves - September 1995
%	Copyright (c) 1995 Rice University
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


if (nargin == 0),
 error ( 'The number of parameters must be at least 1.' );
elseif (nargin == 1),
 t0=N/2; h=0.0;
elseif (nargin ==2),
 h=0.0;
end;

if (N<=0),
 error('N must be greater or equal to 1.');
elseif (t0<=0) | (t0>N),
 error('t0 must be between 1 and N.');
else
 if h <= 0 
   f = (1/N:1/N:0.5-1/N) ;
   y = zeros(1,N/2);
   y(2:N/2) = (f.^(-1-h)).*exp(-i*2*pi*f.*(t0-1));
   x = real(ifft(y,N)) ;
   x = x./max(x); 
   x = x.' - sign(min(x))*abs(min(x)) ;
 else
   t = 1:N;
   x = abs(t-t0).^h;
   x = max(x)-x.' ;
 end
end

x=hilbert(x);