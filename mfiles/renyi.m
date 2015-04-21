function r=renyi(tfr,t,f,alpha);
%RENYI	Measure Renyi information.
%	R=RENYI(TFR,T,F,ALPHA) measures the Renyi information relative 
%	to a 2-D density function TFR (which can be eventually a TF
%	representation).
%
%	TFR : (M,N) 2-D density function (or mass function). Eventually
%	     TFR can be a time-frequency representation, in which case
%	     its first row must correspond to the lower frequencies
%	T : abscissa vector parametrizing the TFR matrix. T can be a
%	    non-uniform sampled vector (eventually a time vector)
%						(default : (1:N)).	
%	F : ordinate vector parametrizing the TFR matrix. F can be a
%	    non-uniform sampled vector (eventually a frequency vector)
%						(default : (1:M)).	
%	ALPHA : rank of the Renyi measure	(default : 3).
%	R : the alpha-rank Renyi measure (in bits if TFR is a time- 
%	    frequency matrix) :
%	     R=log2[Sum[TFR(Fi,Ti)^ALPHA dFi.dTi]/(1-ALPHA)]
%                  Fi,Ti
%
%	Example :
%	 s=atoms(64,[32,.3,16,1]); [TFR,T,F]=tfrsp(s); R=renyi(TFR,T,F,3) 
%	 s=atoms(64,[16,.2,10,1;40,.4,12,1]); [TFR,T,F]=tfrsp(s); 
%	 R=renyi(TFR,T,F,3) 

%	P. Goncalves, October 95
%	Copyright (c) 1995 Rice University.
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

if (nargin == 0),
 error('At least one parameter required');
end;
[M,N] = size(tfr);

if (nargin == 1),
 t=1:N; f=(1:M)'; alpha=3;
elseif (nargin == 2),
 f=(1:M)'; alpha=3;
elseif (nargin == 3),
 alpha=3;
end;

f=sort(f);
tfr = tfr./integ2d(tfr,t,f);
if alpha == 1
 if (min(min(tfr))<0),
  error('distribution with negative values => alpha=1 not allowed');
 else
  r=-integ2d(tfr.*log2(tfr+eps),t,f);
 end
else
 r=log2(integ2d(tfr.^alpha,t,f)+eps)/(1-alpha) ;
end

