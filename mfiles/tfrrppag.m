function [tfr,rtfr,hat] = tfrrppag(x,t,N,h,trace);
%TFRRPPAG Reassigned pseudo Page time-frequency distribution.
%	[TFR,RTFR,HAT] = TFRRPPAG(X,T,N,H,TRACE) computes the pseudo 
%	Page distribution and its reassigned version.
%	When invoked without left hand arguments, TFRRPPAG runs TFRQVIEW.
% 
%	X     : analysed signal,
%	T     : time instant(s)          (default : 1:length(X)),
%	N     : number of frequency bins (default : length(X)).
%	H     : frequency smoothing window, H(0) being forced to 1
%	                                 (default : Hamming(N/4)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR,  : time-frequency representation and its reassigned
%	RTFR    version. When called without output arguments, 
%	        TFRRPMH runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example:
%	 sig=fmconst(128,0.2); t=1:128; 
%	 h=tftb_window(65,'gauss'); tfrrppag(sig,t,128,h,1);
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-July 1994, July 1995.
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

if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);

if (nargin <= 2),
 N=xrow;
end;

hlength=floor(N/4);
if (rem(hlength,2)==0),
 hlength=hlength+1;
end;

if (nargin == 1),
 t=1:xrow; h = tftb_window(hlength); trace=0;
elseif (nargin == 2)|(nargin == 3),
 h = tftb_window(hlength); trace=0;
elseif (nargin == 4),
 trace = 0;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must only have one column');
elseif (trow~=1),
 error('T must have only one row'); 
elseif (2^nextpow2(N)~=N & nargin==5),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
end;

if (tcol==1),
 Dt=1; 
else
 Deltat=t(2:tcol)-t(1:tcol-1); 
 Mini=min(Deltat); Maxi=max(Deltat);
 if (Mini~=Maxi),
  error('The time instants must be regularly sampled.');
 else
  Dt=Mini;
 end;
 clear Deltat Mini Maxi;
end;

tfr= zeros(N,tcol); tf2= zeros(N,tcol);
if trace, disp('Pseudo Page distribution'); end;
Dh=dwindow(h); % .*[zeros(Lh,1);ones(Lh+1,1)]);
for icol=1:tcol,
 ti= t(icol); tau=0:min([N-1,Lh,ti-1]);
 indices= rem(N+tau,N)+1;
 if trace, disprog(icol,tcol,10); end;
 tfr(indices,icol)= h(Lh+1+tau).*x(ti).*conj(x(ti-tau));
 tf2(indices,icol)=Dh(Lh+1+tau).*x(ti).*conj(x(ti-tau));
 tf2(1,icol)=tf2(1,icol)+x(ti).*conj(x(ti));
end ;

tfr= fft(tfr); tf2=fft(tf2);
tfr=tfr(:);tf2=tf2(:);
avoid_warn=find(tfr~=0); 
tf2(avoid_warn)= round(N*tf2(avoid_warn)./tfr(avoid_warn)/(2.0*pi)); 
tfr=real(tfr); tf2=imag(tf2);
if trace, fprintf ('\nreassignment: \n'); end;
tfr=reshape(tfr,N,tcol);
tf2=reshape(tf2,N,tcol);

rtfr= zeros(N,tcol); 
Ex=mean(abs(x(min(t):max(t))).^2); Threshold=1.0e-6*Ex;
for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 for jcol=1:N,
  if abs(tfr(jcol,icol))>Threshold,
   jcolhat= jcol - tf2(jcol,icol);
   jcolhat=rem(rem(jcolhat-1,N)+N,N)+1;
   rtfr(jcolhat,icol)=rtfr(jcolhat,icol) + tfr(jcol,icol) ;
   tf2(jcol,icol)=jcolhat;
  else 
   tf2(jcol,icol)=inf;
   rtfr(jcol,icol)=rtfr(jcol,icol) + tfr(jcol,icol) ;
  end;
 end;
end;

if trace, fprintf('\n'); end;
if (nargout==0),
 TFTBcontinue=1;
 while (TFTBcontinue==1),
  choice=menu ('Choose the representation:',...
               'stop',...
               'pseudo Page distribution',...
               'reassigned pseudo Page distribution');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   tfrqview(tfr,x,t,'tfrppage',h);
  elseif (choice==3),
   tfrqview(rtfr,x,t,'tfrrppag',h);
  end;
 end;
elseif (nargout>2),
 hat=tf2;
end;
