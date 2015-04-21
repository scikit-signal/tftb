function y=lambdak(u,k);
%LAMBDAK Evaluate lambda function for Affine Wigner distribution.
%	Y=LAMBDAK(U,K) evaluates the parametrization lambda function
%	involved in the affine smoothed pseudo Bertrand distribution.
%	
%	 LAMBDAK(U,0) = -U/(exp(-U)-1) for K = 0
%	 LAMBDAK(U,1) = exp(1+U exp(-U)/(exp(-U)-1)) for K = 1
%	 LAMBDAK(U,K) = (K (exp(-U)-1)/(exp(-KU)-1))^(1/(K-1)) otherwise
%
%	U : real vector
%	Y : value of LAMBDAD at point(s) U
%
%	See also LAMBDAB, LAMBDAU, LAMBDAD.

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

if (k~=0) & (k~=1)
  ind = find(u ~= 0) ;
  y = ones(size(u)) ;
  y(ind) = (k*(exp(-u(ind))-1)./(exp(-k*u(ind))-1)).^(1/(k-1)) ;
elseif k==0
  ind = find(u ~= 0) ;
  y = ones(size(u)) ;
  y(ind) = -u(ind)./(exp(-u(ind))-1) ;
elseif k==1
  ind = find(u ~= 0) ;
  y = ones(size(u)) ;
  y(ind) = exp(1+u(ind).*exp(-u(ind))./(exp(-u(ind))-1)) ;
end;
