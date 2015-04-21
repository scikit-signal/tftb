function h=holder(tfr,f,n1,n2,t,pl)
%HOLDER	Estimate the Holder exponent through an affine TFR.
%	H=HOLDER(TFR,F,N1,N2,T) estimates the Holder exponent of a 
%	function through an affine time-frequency representation of it, 
%	and plots the frequency marginal and the regression line.  
%
%	TFR : affine time-frequency representation.  
%	F   : frequency values of the spectral analysis. 
%	N1  : indice of the  minimum frequency for the linear regression.  
%                                         (default : 1).
%	N2  : indice of the  maximum frequency for the linear regression.  
%                                         (default : length(F)).
%	T : time vector. If T is omitted, the function returns the
%	    global estimate of the Holder exponent. Otherwise, it
%	    returns the local estimates H(T) at the instants specified
%	    in T.  
%	H : output value (if T omitted) or vector (otherwise) containing
%	    the Holder estimate(s).
%
%	Example :
%	 S=altes(128); [TFR,T,F]=tfrscalo(S,1:128,8);
%	 H=holder(TFR,F,1,length(F));

%	P. Goncalves, October 1995
%	Copyright (c) 1995 Rice University
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

clf;
tfr = abs(tfr) ;

if nargin<2,
 error('There must be at least 2 input parameters');
elseif nargin==2,
 n1=1; n2=length(f);
elseif nargin==3,
 n2=length(f);
end
if nargin<=5,
 pl=1;
else
 pl=0;
end

if nargin <= 4
  fmarg = mean(tfr.').' ;
  if pl, plot(log(f),log(fmarg)) , hold on; end;
  p = polyfit(log(f(n1:n2)),log(fmarg(n1:n2)),1) ;
  drt = p(1)*log(f)+p(2) ;
  if pl,
   plot(log(f),drt,'--g') ; 
   legend('Frequency marginal','Regression line');
   plot(log(f([n1 n2])),log(fmarg([n1 n2])),'+r') ;
   xticklabels = round(logspace(log10(min(f)),log10(max(f)),4)*1000)./1000 ;
   set(gca,'XTick',log(xticklabels)) ; 
   xlabel('frequency (logarithmically spaced)') ;
   hold off ;
  end
  h = (-p(1)-1)/2 ;
else
  if length(t) == 1
    if pl, plot(log(f),log(tfr(:,t))) ; hold on; end
    p = polyfit(log(f(n1:n2)),log(tfr(n1:n2,t)),1) ;
    drt = p(1)*log(f)+p(2) ;
    if pl, 
     plot(log(f),drt,'--g') ;
     legend('Frequency marginal','Regression line');
     plot(log(f([n1 n2])),log(tfr([n1 n2],t)),'+r') ;
     xticklabels = round(logspace(log10(min(f)),log10(max(f)),4)*1000)./1000 ;
     set(gca,'XTick',log(xticklabels)) ; 
     xlabel('frequency (logarithmically spaced)') ;
     hold off ;
    end
    h = (-p(1)-1)/2 ;
  elseif length(t) > 1 
    [yt,xt] = size(t) ; if yt>xt, t = t.' ; end ;
    j = 1 ;
    for k = t,
      p = polyfit(log(f(n1:n2)),log(tfr(n1:n2,k)),1) ;
      h(j) = (-p(1)-1)/2 ; j = j + 1 ;
    end
    if pl, 
     plot(t,h) ;grid
     title('Holder estimates at time instants T');
    end
  end,
end
