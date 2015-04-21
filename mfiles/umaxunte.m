function y=umaxunte(u);
%UMAXUNTE Determination of the maximum value of u for Unterberger distribution.
%	Y=UMAXUNTE(u) is the function Y(u)=(H(u)+u/2)/(H(u)-u/2)-fmax/fmin. 
%	Doing UMAX = fzero('umaxunte',0); gives the maximum value for U in the
%	computation of the Unterberger distribution. For this distribution, 
%	 	 	H(u) = sqrt(1+(u/2)^2).
%
%	U     : real vector
%	Y     : value of the function (H(u)+u/2)/(H(u)-u/2)-fmax/fmin.

%	P. Goncalves, October 95 - O. Lemoine, July 1996. 
%	Copyright (c) 1995 Rice University - CNRS (France).
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

global ratio_f

y = (sqrt(1+(u/2)^2)+u/2)/(sqrt(1+(u/2)^2)-u/2)-ratio_f;