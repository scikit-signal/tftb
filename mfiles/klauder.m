function x=klauder(N,lambda,f0)
%KLAUDER Klauder wavelet in time domain.
%	X=KLAUDER(N,LAMBDA,F0) generates the KLAUDER wavelet 
%	in the time domain
%	K(f) = e^{-2.pi.LAMBA.f} f^{2.pi.LAMBDA.F0-1/2} 
%
%	N      : number of points in time   
%	LAMBDA : attenuation factor or the envelope (default : 10)
%	F0     : central frequency of the wavelet (default : 0.2)
%	X      : time row vector containing the klauder samples.
%
%	Example :   
%	 x=klauder(128); plot(x);
%
%	See also  ALTES, ANASING, DOPPLER, ANAFSK, ANASTEP. 

%	P. Goncalves 9-95
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
 f0=0.2; lambda=10;
elseif (nargin ==2),
 f0=0.2;
end;

if (N<=0),
 error('N must be greater or equal to 1.');
elseif (f0 > 0.5) | (f0 < 0),
 error('f0 must be between 0 and 0.5') ;
else
 f = linspace(0,0.5,N/2+1) ;
 mod = exp(-2*pi*lambda*f).*f.^(2*pi*lambda*f0-0.5) ;
 wave = mod; wave(1)=0 ; 
 wave = [wave(1:N/2) fliplr(wave(2:N/2+1))] ;
 wavet=ifft(wave) ;
 wavet = fftshift(wavet) ; 
 x = real(wavet).'/norm(wavet) ; 
end

