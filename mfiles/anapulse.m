function y=anapulse(N,ti);
%ANAPULSE Analytic projection of unit amplitude impulse signal.
%	y=ANAPULSE(N,TI) returns an analytic N-dimensional signal 
%	whose real part is a Dirac impulse at t=TI.
%
% 	N  : number of points.
%	TI : time position of the impulse	(default : round(N/2)). 
%	
%	Example :
%	 signal=2.5*anapulse(512,301);plot(real(signal));
%
%	See also ANASTEP, ANASING, ANABPSK, ANAFSK.

%	O. Lemoine - June 1995, F. Auger, August 1995.
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

if (nargin==0),
 error('at least one parameter required');
elseif (nargin==1),
 ti=round(N/2);
end;

if (N<=0),
 error('N must be greater than zero');
else
 t=(1:N)';
 y=hilbert(t==ti);
end;
