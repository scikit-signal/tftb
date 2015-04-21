function tifd=friedman(tfr,hat,t,method,trace);
%FRIEDMAN Instantaneous frequency density.
%	TIFD = FRIEDMAN(TFR,HAT,T,METHOD,TRACE) computes the
%	time-instantaneous frequency density (defined by Friedman [1])
%	of a reassigned time-frequency representation.
% 
%	TFR   : time-frequency representation, (N,M) matrix.
%	HAT   : complex matrix of the reassignment vectors.
%	T     : the time instant(s)	(default : (1:M)).
%	METHOD: chosen representation	(default : 'tfrrsp').  
%	TRACE : if nonzero, the progression of the algorithm is shown
%					(default : 0).
%	TIFD  : time instantaneous-frequency density. When called without 
%	        output arguments, FRIEDMAN runs TFRQVIEW.
%
%	WARNING : TIFD is not an energy distribution, but an estimated 
%	-------        probability distribution !
%
%	Example : 
%	 sig=fmlin(128,0.1,0.4); h=tftb_window(47,'Kaiser');
%	 t=1:2:127; [tfr,rtfr,hat]=tfrrpwv(sig,t,128,h);
%	 friedman(tfr,hat,t,'tfrrpwv',1); 
%
%	See also : RIDGES.

%	F. Auger, August 1994, Decembre 1995.
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
%
%	[1] : D. H. Friedman, "Instantaneous Frequency vs Time : An
%	      Interpretation of the Phase Structure of Speech", Proc. IEEE
%	      ICASSP, pp. 29.10.1-4, Tampa, 1985.	

if (nargin<3),
 error('At least 2 parameters required'); 
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
 error('tfr and hat must have the same size');
end;

tifd=zeros(tfrrow,tfrcol);
bins=0.5+(0:tfrrow-1);
threshold=sum(sum(tfr))*0.5/(tfrrow*tfrcol);

if trace, fprintf ('\nFriedman distribution: \n'); end;

for j=1:tfrcol,
 if trace, disprog(j,tfrcol,10); end;
 indices=find(tfr(:,j)>threshold);
 if (length(indices)>=1),
  [occurences,trash]=hist(real(hat(indices,j)),bins);
  tifd(:,j)=occurences';
 end;
end; 
tifd=tifd/sum(sum(tifd));

method=upper(method);
if nargout==0,
 tfrqview(tifd,[],t,method);
end
