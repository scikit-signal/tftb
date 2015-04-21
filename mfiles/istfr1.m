function flag=istfr1(method);
% flag=istfr1(method) returns true is method is a
% time frequency representation of type 1 (positive and negative frequencies).
%	See also istfr2, istfraff.

%	F. Auger, may 98
%	Copyright (c) CNRS - France 1998. 
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

method=upper(method);
if strcmp(method,'TFRPMH'  )| strcmp(method,'TFRRPMH' )| ...
   strcmp(method,'TFRSP'   )| strcmp(method,'TFRRSP'  )| ...
   strcmp(method,'TFRPPAGE')| strcmp(method,'TFRRPPAG')| ...
   strcmp(method,'TFRMHS'  )| strcmp(method,'TFRRGAB' )| ...
   strcmp(method,'TFRMH'   )| strcmp(method,'TFRMMCE' )| ...
   strcmp(method,'TFRRMSC' )| strcmp(method,'TFRPAGE' )| ...
   strcmp(method,'TFRGABOR')| strcmp(method,'TFRRI'   )| ...
   strcmp(method,'TFRMSC'  )| strcmp(method,'TYPE1'   )| ...
   strcmp(method,'TFRSTFT' ),
 flag=1;
else
 flag=0;
end;

