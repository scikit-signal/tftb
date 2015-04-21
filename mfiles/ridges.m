function [pointst,pointsf]=ridges(tfr,hat,t,method,trace);
%RIDGES	Extraction of ridges.
%	[POINTST,POINTSF]=RIDGES(TFR,HAT,T,METHOD,TRACE) extracts the
%	ridges of a time-frequency distribution. These ridges are some
%	particular sets of curves deduced from the stationary points of
%	their  reassignment operators.
%
%	TFR    : time-frequency representation
%	HAT    : complex matrix of the reassignment vectors.
%	T      : the time instant(s).
%	METHOD : the chosen representation (default: 'tfrrsp'). 
%	TRACE  : if nonzero, the progression of the algorithm is shown
%					   (default : 0).
%
%	POINTST,POINTSF are two vectors for the time and frequency 
%	coordinates of the stationary points of the reassignment. 
%	Therefore, PLOT(POINTST,POINTSF,'.') shows the squeleton of the 
%	representation.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); g=tftb_window(21,'kaiser'); 
%	 h=tftb_window(47,'Kaiser'); t=1:2:127; 
%	 figure(1), [tfr,rtfr,hat]=tfrrspwv(sig,t,128,g,h); 
%	 ridges(tfr,hat,t,'tfrrspwv',1);
%	 figure(2), [tfr,rtfr,hat]=  tfrrsp(sig,t,128,h);   
%	 ridges(tfr,hat,t,'tfrrsp',1);
%
%	See also : FRIEDMAN.

%	F. Auger, August 1994, December 1995.
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

if (nargin<2),
 error('at least 2 parameters required'); 
end;

[tfrrow,tfrcol]=size(tfr);
[hatrow,hatcol]=size(hat);

if (nargin==2),
 t=1:tfrcol; method='tfrrsp'; trace=0;
elseif (nargin==3),
 method='tfrrsp'; trace=0; 
elseif (nargin==4),
 trace=0; 
end;

[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row'); 
elseif (tfrrow~=hatrow)|(tfrcol~=hatcol),
 error('TFR and HAT must have the same size');
end;

Nf=tfrrow; frequencies=(1:Nf)'; 
threshold=sum(sum(tfr))*0.5/(tfrrow*tfrcol);

pointst=[]; pointsf=[];

if trace, fprintf ('\nRidge extraction: \n'); end;

method=upper(method);
if strcmp(method,'TFRRPWV') | strcmp(method,'TFRRPMH'),
 for icol=1:tfrcol, ti=t(icol); 
  if trace, disprog(icol,tfrcol,10); end;
  indices=find((tfr(:,icol)>threshold)&(hat(:,icol)-frequencies==0)); 
  nbindices=length(indices);
  if (nbindices>0), 
   pointst=[pointst;ones(nbindices,1)*ti];
   pointsf=[pointsf;indices/(2.0*Nf)]; 
  end;
 end;
elseif strcmp(method,'TFRRSPWV'),
 for icol=1:tfrcol, ti=t(icol); 
  if trace, disprog(icol,tfrcol,10); end;
  indices=find((real(hat(:,icol))-frequencies==0)&...
               (imag(hat(:,icol))-icol==0)&...
               (tfr(:,icol)>threshold)); 
  nbindices=length(indices);
  if (nbindices>0), 
   pointst=[pointst;ones(nbindices,1)*ti];
   pointsf=[pointsf;indices/(2.0*Nf)]; 
  end;
 end;
elseif strcmp(method,'TFRRSP')|strcmp(method,'TYPE1')
 for icol=1:tfrcol, ti=t(icol); 
  if trace, disprog(icol,tfrcol,10); end;
  indices=find((real(hat(:,icol))-frequencies==0)&...
               (imag(hat(:,icol))-icol==0)&...
               (tfr(:,icol)>threshold)); 
  nbindices=length(indices);
  if (nbindices>0), 
   pointst=[pointst;ones(nbindices,1)*ti];
   pointsf=[pointsf;indices/Nf]; 
  end;
 end;
else 
 error('unknown representation');
end;

if (nargout==0),
 clf;
 plot(pointst,pointsf,'.')  
 axis([min(t) max(t) 0 0.5]); 
end;
