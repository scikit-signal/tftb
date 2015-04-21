function [x,iflaw]=fmpar(N,p1,p2,p3)
%FMPAR	Parabolic frequency modulated signal.
%	[X,IFLAW]=FMPAR(N,P1,P2,P3) generates a signal with
%	parabolic frequency modulation law.
%	X(T) = exp(j*2*pi(A0.T + A1/2.T^2 +A2/3.T^3)) 
%
%	N  : the number of points in time
%	P1 : if NARGIN=2, P1 is a vector containing the three 
%	    coefficients [A0 A1 A2] of the polynomial instantaneous phase.
%	    If NARGIN=4, P1 (as P2 and P3) is a time-frequency point of 
%	    the form [Ti Fi].
%	    The coefficients (A0,A1,A2) are then deduced such that  
%	    the frequency modulation law fits these three points.
%	P2,P3 : same as P1 if NARGIN=4.       (optional)
%	X     : time row vector containing the modulated signal samples 
%	IFLAW : instantaneous frequency law
%
%	Examples :   
%	 [X,IFLAW]=fmpar(128,[1 0.4],[64 0.05],[128 0.4]);
%	 subplot(211);plot(real(X));subplot(212);plot(IFLAW);
%	 [X,IFLAW]=fmpar(128,[0.4 -0.0112 8.6806e-05]);
%	 subplot(211);plot(real(X));subplot(212);plot(IFLAW);
%
%	See also FMCONST, FMHYP, FMLIN, FMSIN, FMODANY, FMPOWER.

%	P. Goncalves - October 1995, O. Lemoine - November 1995
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

if (nargin <= 1),
 error ( 'The number of parameters must be at least 2.' );
elseif (N <= 0),
 error ('The signal length N must be strictly positive' );
elseif nargin == 2 ;
  if length(p1) ~= 3
    error('Bad number of coefficients for P1');
  end
  a0 = p1(1) ; a1 = p1(2) ; a2 = p1(3) ;
elseif nargin == 4 ;
  if (length(p1) ~= 2) |(length(p2) ~= 2) |(length(p3) ~= 2),
    error('Bad number of coefficients for P1, P2, P3');
  end
  if p1(1)>N | p1(1)<1,
   error ('P1(1) must be between 1 and N');
  elseif p2(1)>N | p2(1)<1,
   error ('P2(1) must be between 1 and N');
  elseif p3(1)>N | p3(1)<1,
   error ('P3(1) must be between 1 and N');
  elseif p1(2)<0,
   error ('P1(2) must be > 0');
  elseif p2(2)<0,
   error ('P2(2) must be > 0');
  elseif p3(2)<0,
   error ('P3(2) must be > 0');
  end
  Y = [p1(2) p2(2) p3(2)] ;
  X = [1 1 1;p1(1) p2(1) p3(1);p1(1)^2 p2(1)^2 p3(1)^2] ;
  coef = Y*inv(X) ; 
  a0 = coef(1) ;
  a1 = coef(2) ;
  a2 = coef(3) ;
end  

t=1:N;

phi = 2*pi*(a0*t + a1/2*t.^2 + a2/3*t.^3) ;
iflaw = (a0 + a1*t + a2*t.^2).' ;

aliasing = find(iflaw<0 | iflaw>0.5) ;
if isempty(aliasing) == 0
  disp(['!!! WARNING: signal is undersampled or has negative frequencies']) ;
end

x = exp(i*phi).';

