function disprog(i,N,steps);
%DISPROG Display progression of a loop.
%	DISPROG(i,N,steps) displays the progression of a loop.
%
%	I     : loop variable
%	N     : final value of i
%	STEPS : number of displayed steps.
%
%	Example:
%	 N=16; for i=1:N, disprog(i,N,5); end;

%	F. Auger, August, December 1995.
%       from an idea of R. Settineri.
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

global begin_time_disprog ;

if (i==1),
 begin_time_disprog=cputime;
end;

if (i==N),
 fprintf('100 %% complete in %g seconds.\n', cputime-begin_time_disprog);
 clear begin_time_disprog;
elseif (floor(i*steps/N)~=floor((i-1)*steps/N)),
 fprintf('%g ', floor(i*steps/N)*ceil(100/steps));
end;


