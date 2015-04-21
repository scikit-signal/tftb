function [x,gpd,nu]=gdpower(N,k,c)
%GDPOWER Signal with power-law group delay.
%	[X,GPD,F]=GDPOWER(N,K,C) generates a signal with a
%	power-law group delay of the form tx(f) = T0 + C*f^(K-1).
%	The output signal is of unit energy.
%
%	N   : number of points in time		(must be even)
%	K   : degree of the power-law		(default : 0)
%	C   : rate-coefficient of the power-law group delay.  
%	      C must be non-zero. 		(default : 1)	
%	X   : time row vector containing the modulated signal samples
%	GPD : group delay	      (length : round(N/2)-1)
%	F   : frequency bins
%
%	Examples :   
%	 sig=gdpower(64,0); tfrbert(sig,1:64,0.01,0.32,256,1);
%	 [sig,gpd,f]=gdpower(256,1/2); plot(gpd,f); 
%
%	See also FMPOWER.

%	O. Lemoine - April, July 1996
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

if (nargin < 1),
 error ( 'At least one parameter.' );
elseif nargin==1,
 k=0; c=1;
elseif nargin==2,
 c=1;
end

if (N <= 0),
 error ('The signal length N must be strictly positive' );
end  

t0=0;

Lnu=round(N/2);
nu=linspace(0,0.5,Lnu+1); nu=nu(2:Lnu+1);
AM=nu.^((k-2)/6);

if c==0,
 error('C must be non-zero');
end

t=1:N;
TFx=zeros(1,N);

if (k<1 & k~=0),
 d=N^k*c; t0=N/10;
 TFx(1:Lnu)=exp(-j*2*pi*(t0*(nu)+d*(nu).^k/k)).*AM;
 x=ifft(TFx).';
elseif k==0,
 d=c; t0=N/10;
 TFx(1:Lnu)=exp(-j*2*pi*(t0*(nu)+d*log(nu))).*AM;
 x=ifft(TFx).';
elseif k==1,
 t0=N;
 x=anapulse(N,t0);
elseif (k>1),
 d=N*2^(k-1)*c; 
 TFx(1:Lnu)=exp(-j*2*pi*(t0*(nu)+d*(nu).^k/k)).*AM;
 x=ifft(TFx).';
end

if k~=1,
 gpd=t0+abs(sign(c)-1)/2*(N+1)+d*nu.^(k-1);
else
 gpd=t0*ones(1,N/2);
end

x=x-mean(x);
x=x/norm(x);
