function tfrsave(name,tfr,method,sig,t,f,p1,p2,p3,p4,p5);
%TFRSAVE Save the parameters of a time-frequency representation.
%	TFRSAVE(NAME,TFR,METHOD,SIG,T,F,P1,P2,P3,P4,P5) saves the 
%	parameters of a time-frequency representation in the file
%	NAME.mat. Two additional parameters are saved : TfrQView and
%	TfrView. If you load the file 'name.mat' and do eval(TfrQView), you
%	will restart the display session under tfrqview ; if you do
%	eval(TfrView), you will display the representation by means of
%	tfrview. 
% 
%	NAME   : name of the mat-file (less than 8 characters).   
%	TFR    : time-frequency representation (MxN).
%	METHOD : chosen representation.	
%	SIG    : signal from which the TFR was obtained 
%	T      : time instant(s)	   (default : (1:N)).
%	F      : frequency bins		   (default : linspace(0,0.5,M)).
%	P1..P5 : optional parameters : run the file TFRPARAM(METHOD) 
%		 to know the meaning of P1..P5 for your method.  
%
%	Example : 
%	 sig=fmlin(64); tfr=tfrwv(sig);
%	 tfrsave('wigner',tfr,'TFRWV',sig,1:64);  
%	 clear; load wigner; eval(TfrQView);
%
%	See also TFRQVIEW, TFRVIEW, TFRPARAM.

%	F. Auger, August 1994 - O. Lemoine, June 1996.
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

[M,N]=size(tfr);

if (nargin < 3),
 error('At least 3 parameters are required'); 
elseif nargin==3,
 sig=[]; t=1:N; f=(0.5*(0:M-1)/M); 
elseif nargin==4,
 t=1:N; f=(0.5*(0:M-1)/M); 
elseif nargin==5,
 f=(0.5*(0:M-1)/M);
end;

method=upper(method); 
namedflt=name;

while (length(name)>8),
 disp('The name must have less than 8 characters');
 namedflt=name(1:8);
 nameStr=[' Name of the MAT file [',namedflt,'] : '];
 name=input(nameStr,'s'); 
end

if name==[],
 name=namedflt;
end

if (nargin==3),
 t=1:N;
 TfrQView=['tfrqview(tfr,[],t,method)'];
 TfrView =['clf;tfrview(tfr,sig,t,method,param,map)'];
 eval(['save ',name,' tfr t f method TfrQView TfrView']);
elseif (nargin==4),
 TfrQView=['tfrqview(tfr,[],t,method)'];
 TfrView =['clf;tfrview(tfr,sig,t,method,param,map)'];
 eval(['save ',name,' tfr t f method TfrQView TfrView']);
elseif (nargin==5),
 TfrQView=['tfrqview(tfr,sig,t,method)'];
 TfrView =['clf;tfrview(tfr,sig,t,method,param,map)'];
 eval(['save ',name,' tfr sig t f method TfrQView TfrView']);
elseif (nargin==6),
 TfrQView=['tfrqview(tfr,sig,t,method,p1)'];
 TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1)'];
 eval(['save ',name,' tfr sig t f method p1 TfrQView TfrView']);
elseif (nargin==7),
 TfrQView=['tfrqview(tfr,sig,t,method,p1,p2)'];
 TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2)'];
 eval(['save ',name,' tfr sig t f method p1 p2 TfrQView TfrView']);
elseif (nargin==8),
 TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3)'];
 TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3)'];
 eval(['save ',name,' tfr sig t f method p1 p2 p3 TfrQView TfrView']);
elseif (nargin==9),
 TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3,p4)'];
 TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4)'];
 eval(['save ',name,' tfr sig t f method p1 p2 p3 p4 TfrQView TfrView']);
elseif (nargin==10),
 TfrQView=['tfrqview(tfr,sig,t,method,p1,p2,p3,p4,p5)'];
 TfrView =['clf; tfrview(tfr,sig,t,method,param,map,p1,p2,p3,p4,p5)'];
 eval(['save ',name,' tfr sig t f method p1 p2 p3 p4 p5 TfrQView TfrView']);
end;
