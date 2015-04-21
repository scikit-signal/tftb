function [tfr,t,f] = tfrparam(x,t,N,p,L,trace,Rxtype,method,q);
% [tfr,t,f]=tfrparam(x,t,N,p,L,trace,Rxtype,method,q) parametric time-frequency 
% representation.
%
% X      : signal
% T      : time instant(s)          (default : 1:length(X)).
% N      : number of frequency bins (default : max(256,length(X)) ).
% p      : last autocorrelation lag (default : 2 ).
% L      : length of the window around the analyzed time sample. Must be odd.
%                                   (default : max(51,round(xrow/10)) ).
% Rxtype : choice of the correlation matrix algorithm (default : 'fbhermitian')
%          possible values are 'hermitian', 'fbhermitian', 'burg' or 'fbburg'
% method : can be either 'ar', 'periodogram', 'capon', 'capnorm', 'lagunas',
%          or 'genlag'.
% q      : parameter for the generalized Lagunas method.
%
% example:
%
%          x=fmconst(256,0.1)+fmlin(256,0.15,0.30)+fmlin(256,0.2,0.45);
%          figure(1); tfrparam(x,1:256,256,3,21,1,'fbhermitian','ar');
%          figure(2); tfrparam(x,1:256,256,15,61,1,'fbhermitian','periodogram');
%          figure(3); tfrparam(x,1:256,256,3,21,1,'fbhermitian','capon');
%          figure(4); tfrparam(x,1:256,256,3,21,1,'fbhermitian','lagunas');

% F. Auger, july, november 1998, april 99.
if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);
if (nargin==1),
 t=1:xrow; N=max(256,xrow); p=2; L=max(51,round(xrow/10)); trace=0; Rxtype='fbhermitian'; method='ar';
elseif (nargin==2),
 N=max(256,xrow); p=2; L=max(51,round(xrow/10)); trace=0; Rxtype='fbhermitian'; method='ar';
elseif (nargin==3),
 p=2; L=max(51,round(xrow/10)); trace=0; Rxtype='fbhermitian'; method='ar';
elseif (nargin==4),
 L=max(51,round(xrow/10)); trace=0; Rxtype='fbhermitian'; method='ar';
elseif (nargin==5),
 trace=0; Rxtype='fbhermitian'; method='ar';
elseif (nargin==6),
 Rxtype='fbhermitian'; method='ar';
elseif (nargin==7),
 method='ar';
end;
 
[trow,tcol] = size(t);

if xcol>1,
 error('x must be a column vector');
elseif p<1,
 error('p must be greater than 0');
end;

if (rem(L,2)==0), % L must be odd
 L=L+1; fprintf('L corrected to L+1\n');
end;
Lhalf=(L-1)/2;

tfr= zeros (N,tcol) ;  
if trace, disp('Tfrparam'); end;
for icol=1:tcol,
 ti= t(icol);
 if ti-Lhalf >= 1,    timin=ti-Lhalf; else timin=1; end;
 if ti+Lhalf <= xrow, timax=ti+Lhalf; else timax=xrow; end;
 if trace, disprog(icol,tcol,10); end;
 Realp=min (p,timax-timin); 
 Rx=correlmx(x(timin:timax),Realp,Rxtype);
 if strcmp(upper(method),'GENLAG'),
  [tfr(:,icol),f]=parafrep(Rx,N,method,q);
 else
  [tfr(:,icol),f]=parafrep(Rx,N,method);
 end
end;
if trace, fprintf('\n'); end;

if (nargout==0),
 imagesc(t,f,10.0*log10(tfr)); axis('xy'); 
 title(['tfrparam, ' method, ...
        ', p=', num2str(p), ...
        ', L=', num2str(L), ', ', ...
        Rxtype, ' autocorrelation matrix']); drawnow;
end;


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
