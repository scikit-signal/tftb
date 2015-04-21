function flag=istfr2(method);
% flag=istfr2(method) returns true is method is a
% time frequency representation of type 2 (only positive frequencies).
%	See also istfr1, istfraff.

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
if strcmp(method,'TFRWV'   ) | strcmp(method,'TFRPWV'  ) | ...
   strcmp(method,'TFRSPWV' ) | strcmp(method,'TFRCW'   ) | ...
   strcmp(method,'TFRZAM'  ) | strcmp(method,'TFRBJ'   ) | ...
   strcmp(method,'TFRBUD'  ) | strcmp(method,'TFRGRD'  ) | ...
   strcmp(method,'TFRRSPWV') | strcmp(method,'TFRRPWV' ) | ...
   strcmp(method,'TFRRIDB' ) | strcmp(method,'TFRRIDH' ) | ...
   strcmp(method,'TFRRIDT' ) | strcmp(method,'TFRASPW' ) | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPAW' ) | ...
   strcmp(method,'TFRRIDBN') | strcmp(method,'TFRUNTER') | ...
   strcmp(method,'TFRBERT' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TFRSPBK' ) | strcmp(method,'TYPE2'   ),
 flag=1;
else
 flag=0;
end;

