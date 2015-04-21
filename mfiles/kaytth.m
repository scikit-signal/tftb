function H=kaytth(length);
%KAYTTH	 Kay-Tretter filter computation. 
%	H=KAYTTH(length); Kay-Tretter filter computation.
% 
%	See also INSTFREQ

%	F. Auger, March 1994.
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

pp1=length*(length+1);
den=2.0*length*(length+1)*(2.0*length+1.0)/3.0;
i=1:length; H=pp1-i.*(i-1);

H=H ./ den;

