function y = tffilter(tfr,x,t,trace);
%TFFILTER Time frequency filtering of a signal.
%       Y=TFFILTER(TFR,X,T,TRACE) filters the signal X
%	with a non stationary filter. 
%	
% 
%	X     : input signal (must be analytic).
%	T     : time instant(s)          (default : 1:length(X)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : Wigner-Ville distribution of the filter
%	        frequency axis is graduated from 0.0 to 0.5.
%
%	Example 1:
%        Nt=128; t=1:Nt; sig=fmlin(Nt,0.05,0.3)+fmlin(Nt,0.2,0.45); 
%        sig(Nt/2)=sig(Nt/2)+8; figure(1);tfrwv(sig,t);
%        Nf=128;freqs=0.5*(0:Nf-1).'/Nf; 
%        for tloop=1:Nt, 
%        rate=0.2*(tloop-1)/(Nt-1);
%        H(:,tloop)=(0+rate<freqs).*(freqs<0.1+rate); 
%        end;
%        y=tffilter(H,sig,t,1);figure(2); tfrwv(y,t);
%
%       Example 2:
%
%        Nt=128; t=1:Nt; sig=atoms(128,[64 0.25 sqrt(128) 1]);
%        figure(1);tfrwv(sig,t);
%        Nf=64;
%        H=zeros(Nf,Nt);H(Nf/4+(-15:15),Nt/2+(-15:15))=ones(31);
%        y=tffilter(H,sig,t,1);figure(2); tfrwv(y,t);
%
% 
%	See also: tfristft

%	F. Auger, Jan 1997.
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

% this function is based on the following reference :
% W. Kozek, F. Hlawatsch, A comparative study of linear and non-linear
% time-frequency filters, proc IEEE int symp on time-frequency and 
% time-scale analysis, pp 163-166, victoria, canada, 1992.

if (nargin<3),
 error('At least 3 parameters required');
elseif (nargin==3),
 trace=0;
end;

[N,NbPoints]=size(tfr);
[trow,tcol] =size(t);
[xrow,xcol] =size(x);

if (trow~=1),
 error('T must only have one row'); 
elseif (xcol~=1),
 error('X must only have one column'); 
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
elseif (xrow~=NbPoints)
 error('tfr should have as many columns as X has rows.');
elseif (min(t)<1)|(max(t)>NbPoints),
 error('the values of T must be between 1 and xrow.');
end; 

if trace, disp('non stationary signal filtering'); end;

tfr=ifft(tfr);
y=zeros(tcol,1);

for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 ti=t(icol); 
 valuestj=max([1 2-ti ti-N/2+1]):min([NbPoints 2*NbPoints-ti ti+N/2]);
 %fprintf('%f %f %f\n',ti, min(valuestj),max(valuestj));

 for tj=valuestj,
  tmid = fix(0.5*(ti+tj)); 
  tdiff= fix(ti-tj); indice= rem(N+tdiff,N)+1;
  %fprintf('%g %g %g\n',tj,tmid,indice);
  y(icol,1)=y(icol,1)+tfr(indice,tmid)*x(tj);
 end;
end;

if trace, fprintf('\n'); end;
