function plotifl(t,iflaws,signal);
%PLOTIFL Plot normalized instantaneous frequency laws.
% 	PLOTIFL(T,IFLAWS,SIGNAL) plot the normalized instantaneous frequency 
% 	laws of each signal component.
%
% 	T      : time instants,
% 	IFLAWS : (M,P)-matrix where each column corresponds to
%                the instantaneous frequency law of an (M,1)-signal,
%	         These P signals do not need to be present at the same 
%                time instants.
%                The values of IFLAWS must be between -0.5 and 0.5.
%       SIGNAL : if the signal is precised, display it also 
%
% 	Example : 
%        N=140; t=0:N-1; [x1,if1]=fmlin(N,0.05,0.3); 
%        [x2,if2]=fmsin(70,0.35,0.45,60);x1(35+(1:70))=x1(35+(1:70))+2*x2;
%        if2=[zeros(35,1)*NaN;if2;zeros(N-70-35,1)*NaN];
%        figure(1); clf; plotifl(t,[if1 if2]);
%        figure(2); clf; plotifl(t,[if1 if2],x1);
%
%	See also TFRIDEAL, PLOTSID.

% 	F. Auger, August 94, August 95, May 1998.
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
 error('2 parameters are required'); 
end;

[trow,tcol] = size(t);
[ifrow,ifcol]=size(iflaws); 

indices=find(1-isnan(iflaws));
maxif=max(max(iflaws(indices))); 
minif=min(min(iflaws(indices)));

if (trow~=1),
 error('t must only have one row'); 
elseif (tcol~=ifrow),
 error('T must have as many lines as iflaws has columns');
elseif (maxif > 0.5) | (minif < -0.5),
 disp('Each element of IFLAWS should be between -0.5 and 0.5'); 
end;

if (nargin==3),
 [Nsig,Ncol]=size(signal);
 if Ncol~=1, 
  error('SIGNAL must have one column'); 
 elseif (Nsig~=tcol),
  error('signal must have as many lines as t has rows');
 end
 clf; set(gcf,'Resize','On','NextPlot','Add');
 axsig = axes('Units','normal','Visible','on','Box','On','position', [0.10 0.69 0.80 0.25]);
 axes(axsig);  plot(t,real(signal)); 
 set(gca, 'XLim',[min(t),max(t)], 'XGrid','on', 'YGrid','on'); title('signal');

 axtfr = axes('Units','normal','Visible','on','Box','On',...
              'position', [0.10 0.12 0.80 0.45]);
 axes(axtfr); 
else
 clf; axes(gca);
end;

plot(t,iflaws); 
if (minif>=0),
 set(gca,'XLim',[min(t),max(t)], 'XGrid','on','YLim',[0 0.5], 'YGrid','on');
else
 set(gca,'XLim',[min(t),max(t)], 'XGrid','on','YLim',[-0.5 0.5], 'YGrid','on');
end;

xlabel('Time');
ylabel('Normalized frequency');
title('Instantaneous frequency law(s)');
