function [ti,fi]=midpoint(t1,f1,t2,f2,k)
%MIDPOINT Mid-point construction used in the interference diagram. 
%	[TI,FI]=MIDPOINT(T1,F1,T2,F2,K) gives the coordinates in the
%	time-frequency plane of the interference-term corresponding to
%	the points (T1,F1) and (T2,F2), for a distribution in the
%	affine class perfectly localized on power-law group-delays of 
%	the form : tx(nu)=t0+c nu^(K-1).
%
%	T1 : time-coordinate of the first point
%	F1 : frequency-coordinate of the first point (>0)
%	T2 : time-coordinate of the second point
%	F2 : frequency-coordinate of the second point (>0)
%	K  : power of the group-delay law
%	  K = 2    : Wigner-Ville 
%	  K = 1/2  : D-Flandrin
%	  K = 0    : Bertrand (unitary) 
%	  K = -1   : Unterberger (active)
%	  K = inf  : Margenau-Hill-Rihaczek
%	TI : time-coordinate of the interference term
%	FI : frequency-coordinate of the interference term
%
%	See also PLOTSID.

%	P. Flandrin, September 1995 - F. Auger, April 1996.
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

if f1<=0 | f2<=0,
 error('F1 and F2 must be >0');
end
[rt1,ct1]=size(t1);
[rt2,ct2]=size(t2);
[rf1,cf1]=size(f1);
[rf2,cf2]=size(f2);
if (rt1~=rt2|rt1~=rf1|rt1~=rf2) | (ct1~=ct2|ct1~=cf1|ct1~=cf2), 
 error('T1, T2, F1 and F2 must have the same size');
end
if rt1>ct1,
 error('T1 must be a row-vector');
elseif rt2>ct2,
 error('T2 must be a row-vector');
elseif rf2>cf2,
 error('F2 must be a row-vector');
elseif rf1>cf1,
 error('F1 must be a row-vector');
end
 
if (k==2),
 fi=(f1+f2)/2;
 ti=(t1+t2)/2;
elseif (k==inf),
 ti=[t1;t2];
 fi=[f2;f1];
else
 I=find(abs(f1-f2)>sqrt(eps));
 if length(I)~=0, 
  if (k==1),
   fi(I)=exp( (f1(I).*(log(f1(I))-1)-f2(I).*(log(f2(I))-1)) ./ ...
 	(f1(I)-f2(I))); 
   ti(I)=(t1(I).*f1(I)-t2(I).*f2(I)) ./ (f1(I)-f2(I)) - ...
 	(t1(I)-t2(I)) ./ (log(f1(I))-log(f2(I)));
  elseif (k==0),
   fi(I)=(f1(I)-f2(I))./(log(f1(I))-log(f2(I)));
   ti(I)=(t1(I).*f1(I)-t2(I).*f2(I)) ./ (f1(I)-f2(I)) + ...
         f1(I) .* f2(I) .* (t2(I)-t1(I)) .* ...
    	(log(f1(I))-log(f2(I))) ./ (f2(I)-f1(I)).^2; 
  else
   t0(I)=(t1(I).*f2(I).^(k-1)-t2(I).*f1(I).^(k-1)) ./ ...
 	(f2(I).^(k-1)-f1(I).^(k-1));
   fi(I)=((f1(I).^k-f2(I).^k) ./ (f1(I)-f2(I))/k).^(1/(k-1));
   ti(I)=t0(I)+(t2(I)-t1(I)) ./ (f2(I).^(k-1)-f1(I).^(k-1)) .*fi(I).^(k-1);
  end
 end
end
