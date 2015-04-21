function som = integ(y,x) 
%INTEG	Approximate integral.
%	SOM=INTEG(Y,X) approximates the integral of vector Y
%	according to X.
%
%	Y   : N-row-vector (or MxN-matrix) to be integrated 
%	      (along each row).  
%	X   : N-row-vector containing the integration path of Y
%				(default : 1:N)
%	SOM : value (or Mx1 vector) of the integral
%
%	Example :    
%	 Y = altes(256,0.1,0.45,10000)'; X = (0:255);
%	 SOM = integ(Y,X)
%
%	See also INTEG2D.

%	P. Goncalves, October 95
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

[M,N]=size(y);

if nargin<1,
 error('At least one parameter required');
elseif nargin==1,
 x=1:N;
end

[Mx,Nx]=size(x);
if (Mx~=1),
 error('X must be a row-vector');
elseif (N~=Nx),
 error('Y must have as many columns as X');
elseif (N==1 & M>1),
 error('Y must be a row-vector or a matrix');
end
 
dy = y(:,1:N-1) + y(:,2:N) ;
dx = (x(2:N)-x(1:N-1))/2 ;
som = dy*dx';
