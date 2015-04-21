function Dh=dwindow(h);
%DWINDOW Derive a window.
%	DH=DWINDOW(H) derives a window H.
%
%	Example : 
%	 plot(dwindow(tftb_window(210,'hanning')))
%
%	See also WINDOW.

%	F. Auger, August 1994, July 1995.
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
 error('one parameter required'); 
end;
[hrow,hcol]=size(h); 
if (hcol~=1),
 error('h must have only one column');
end;

Lh=(hrow-1)/2;
step_height=(h(1)+h(hrow))/2;
ramp=(h(hrow)-h(1))/(hrow-1);
h2=[0;h-step_height-ramp*(-Lh:Lh).';0]; 
Dh=(h2(3:hrow+2)-h2(1:hrow))/2 + ramp; 
Dh(1)   =Dh(1)   +step_height; 
Dh(hrow)=Dh(hrow)-step_height;

