function som=integ2d(mat,x,y);
%INTEG2D Approximate 2-D integral.
%	SOM=INTEG2D(MAT,X,Y) approximates the 2-D integral of  
%	matrix MAT according to abscissa X and ordinate Y.
%
%	MAT : MxN matrix to be integrated
%	X   : N-row-vector indicating the abscissa integration path 	
%				(default : 1:N)
%	Y   : M-column-vector indicating the ordinate integration path 
%				(default : 1:M)	
%	SOM : result of integration
%
%	Example :      
%	 S = altes(256,0.1,0.45,10000) ;
%	 [TFR,T,F] = tfrscalo(S,21:190,8,'auto') ;
%	 E = integ2d(TFR,T,F)
%
%	See also INTEG.

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

[M,N]=size(mat);

if nargin<1,
 error('At least one parameter required');
elseif nargin==1,
 x=1:N; y=(1:M)';
elseif nargin==2,
 y=(1:M)';
end

[xr,xc]=size(x);
[yr,yc]=size(y);

if (xr>xc & xr~=1),
 error('X must be a row-vector');
elseif (yc>yr & yc~=1),
 error('Y must be a column-vector');
elseif (N~=xc),
 error('MAT must have as many columns as X');
elseif (M~=yr),
 error('MAT must have as many rows as Y');
end

mat = (sum(mat.').'-mat(:,1)/2-mat(:,N)/2).*(x(2)-x(1)) ;
dmat = mat(1:M-1)+mat(2:M) ;
dy = (y(2:M)-y(1:M-1))/2 ;
som = sum(dmat.*dy) ;

