function Rx=correlmx(x,p,Rxtype);
% Rx=correlmx(x,p,Rxtype) correlation matrix of a signal
% 
%        Rx : correlation matrix (p+1) x (p+1)
%         x : analyzed signal
%         p : last autocorrelation lag
%    Rxtype : computation algorithm (default : 'fbhermitian')
%             possible values : 'hermitian', 'fbhermitian', 'burg' or 'fbburg'
%
% example :
%
% N=100; sig=real(fmconst(N,0.1))+0.4*randn(N,1); 
% Rx=correlmx(sig,2,'burg'); [v,d] = eig(Rx), acos(-0.5*v(2,1)/v(1,1))/(2*pi)
% Rx=correlmx(sig,2,'hermitian'); [v,d] = eig(Rx), acos(-0.5*v(2,1)/v(1,1))/(2*pi)

% F. Auger, july 1998.

if (nargin<2),
 error('At least two parameters required');
elseif (nargin==2),
 Rxtype='fbhermitian';
end;

[L,xcol]=size(x);
if xcol>1,
 error('x must be a column vector');
elseif p>L,
 error('L must be greater than p');
elseif p<1,
 error('p must be greater than 0');
end;

Rxtype=upper(Rxtype);
if strcmp(Rxtype,'HERMITIAN')|strcmp(Rxtype,'FBHERMITIAN'),
 vector=x(p+1-(0:p)); Rx=conj(vector) * vector.';
 for t=p+2:L,
  vector=x(t-(0:p)); Rx=Rx+conj(vector) * vector.';
 end;
 Rx=Rx/(L-p);

elseif strcmp(Rxtype,'BURG')|strcmp(Rxtype,'FBBURG'),
 R0=sum(abs(x).^2)/L; % variance
 Rpos=zeros(1,p); Rneg=zeros(1,p);
 for n=1:p, 
  Rpos(n)=sum(x(n+1:L).*conj(x(1:L-n)))/(L-n); 
  Rneg(n)=sum(x(1:L-n).*conj(x(n+1:L)))/(L-n);
 end;
 Rx=toeplitz([R0 Rpos],[R0 Rneg]);
else error(['unknown algorithm name' Rxtype]); 
end;

if strcmp(Rxtype(1:2),'FB'),
 Rx=0.5*(Rx+Rx');
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
