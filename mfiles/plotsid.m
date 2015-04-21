function plotsid(t,iflaws,k); 
%PLOTSID Schematic interference diagram of FM signals.  
%	PLOTSID(T,IFLAWS,K) plots the schematic interference diagram of 
%	(analytic) FM signals.  
% 
%	T : time instants, 
%	IFLAWS : matrix of instantaneous frequencies, 
%	         with as may columns as signal components.  
%	K : distribution		(default : 2): 
%	  K = 2     : Wigner-Ville 
%	  K = 1/2   : D-Flandrin
%	  K = 0     : Bertrand (unitary) 
%	  K = -1    : Unterberger (active)
%	  K = inf   : Margenhau-Hill-Rihaczek
% 
%	Example : 
%	 Nt=90; [y,iflaw]=fmlin(Nt,0.05,0.25); 
%	 [y2,iflaw2]=fmconst(50,0.4); 
%	 iflaw(:,2)=[NaN*ones(10,1);iflaw2;NaN*ones(Nt-60,1)]; 
%	 plotsid(1:Nt,iflaw,0); 
% 
%	See also PLOTIFL, MIDPOINT.

%	P. Flandrin, September 1995.
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

if (nargin==1),
 error('at least 2 parameters : t and iflaws');
elseif (nargin==2),
 k=2; % Wigner-Ville
end;

indices=find(1-isnan(iflaws));
if (min(iflaws(indices))<0)|(max(iflaws(indices))>0.5),
 error ('each element of IFLAWS must be between 0 and 0.5');
end;

[iflawrow,iflawcol]=size(iflaws);
tcol=length(t);
clf; figure(gcf); plotifl(t,iflaws); 
hold on;
col=['y','m','c','r'];

% auto-terms
for j=1:iflawcol,
 indices=find(1-isnan(iflaws(:,j)));
 Nbpoints=length(indices);
 for i=1:Nbpoints-1,
  ta=       t(indices(i))  *ones(1,Nbpoints-i); 
  fa=  iflaws(indices(i),j)*ones(1,Nbpoints-i);
  tb= t(indices(i+1:Nbpoints));
  fb=iflaws(indices(i+1:Nbpoints),j)';
  [ti,fi]=midscomp(ta,fa,tb,fb,k);
  plot(ti,fi,['.',num2str(col(rem(j-1,4)+1))]);
 end;
end;

% cross-terms
for j1=1:iflawcol,
 indices1=find(1-isnan(iflaws(:,j1)));
 Nbpoints1=length(indices1);
 for j2=j1+1:iflawcol,
  indices2=find(1-isnan(iflaws(:,j2)));
  Nbpoints2=length(indices2);
  for i=1:Nbpoints1,
   ta=       t(indices1(i))   *ones(1,Nbpoints2); 
   fa=  iflaws(indices1(i),j1)*ones(1,Nbpoints2);
   tb= t(indices2);
   fb=iflaws(indices2,j2)'; 
   [ti,fi]=midscomp(ta,fa,tb,fb,k);
   plot(ti,fi,'.g')   
  end;
 end;
end;

hold off
axis([t(1) t(tcol) 0 0.5]);
grid
if k==2,
 dist=' of the Wigner-Ville distribution';
elseif k==1/2,
 dist=' of the D-Flandrin distribution';
elseif k==0,
 dist=' of the (unitary) Bertrand distribution';
elseif k==-1,
 dist=' of the (active) Unterberger distribution';
elseif k>1/sqrt(eps),
 dist=' of the Margenhau-Hill-Rihaczek distribution';
else
 dist='';
end

title(['Interference diagram',dist,' (k = ',num2str(k),')']);
