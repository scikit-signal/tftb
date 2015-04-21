function [y,iflaw]=fmlin(N,fnormi,fnormf,t0);
%FMLIN	Signal with linear frequency modulation.
%	[Y,IFLAW]=FMLIN(N,FNORMI,FNORMF,T0) generates a linear frequency  
%	modulation.
%	The phase of this modulation is such that Y(T0)=1.
%
%	N       : number of points
%	FNORMI  : initial normalized frequency (default: 0.0)
%	FNORMF  : final   normalized frequency (default: 0.5)
%	T0      : time reference for the phase (default: N/2).
%	Y       : signal
%	IFLAW   : instantaneous frequency law  (optional).
%
%	Example : 
%	 z=amgauss(128,50,40).*fmlin(128,0.05,0.3,50); 
%	 plot(real(z));
%
%	see also FMCONST, FMSIN, FMODANY, FMHYP, FMPAR, FMPOWER.

%	F. Auger, July 1995.
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

if (nargin == 0),
 error ( 'The number of parameters must be at least 1.' );
elseif (nargin == 1),
 fnormi=0.0; fnormf=0.5; t0= round(N/2);
elseif (nargin == 2),
 fnormf=0.5; t0 = round(N/2);
elseif (nargin == 3),
 t0 = round(N/2);
end;

if (N <= 0),
 error ('The signal length N must be strictly positive' );
elseif (abs(fnormi) > 0.5) | (abs(fnormf) > 0.5),
 error ( 'fnormi and fnormf must be between -0.5 and 0.5' ) ;
else
 y = (1:N)';
 y = fnormi*(y-t0) + ((fnormf-fnormi)/(2.0*(N-1))) * ((y-1).^2 - (t0-1).^2);
 y = exp(j*2.0*pi*y) ;
 y=y/y(t0);
 if (nargout==2), iflaw=linspace(fnormi,fnormf,N).'; end;
end ;

