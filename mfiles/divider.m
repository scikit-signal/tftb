function [N,M]=divider(N1);
%DIVIDER Find dividers of an integer.  
%	[N,M]=DIVIDER(N1) find two integers N and M such that M*N=N1 and
%	M and N as close as possible from sqrt(N1).
%
%	Example :
%	 N1=258; [N,M]=divider(N1)

%	F. Auger - November 1995.
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

N=floor(sqrt(N1));
TFTBcontinue=1.0;
while TFTBcontinue,
 Nold=N;
 M=ceil(N1/N);
 N=floor(N1/M);
 TFTBcontinue=(N~=Nold);
end;
