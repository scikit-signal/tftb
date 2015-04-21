function [sig,locatoms]= atoms(N,coord,display);
% ATOMS Linear combination of elementary Gaussian atoms.
%	[SIG,LOCATOMS] = ATOMS(N,COORD,DISPLAY) 
%	generates a signal consisting in a linear combination of elementary
%	gaussian wave packets. The locations of the time-frequency centers 
% 	of the different atoms are either fixed by the input parameter COORD 
%	or successively defined by clicking with the mouse (if NARGIN==1).  
% 
%	N        : number of points of the signal
%	COORD    : matrix of time-frequency centers, of the form
%		   [t1,f1,T1,A1;...;tM,fM,TM,AM]. (ti,fi) are the 
%		   time-frequency coordinates of atom i, Ti is its time 
%		   duration and Ai its amplitude. Frequencies f1..fM should 
%		   be normalized (between 0 and 0.5). 
%		   If nargin==1, the location of the atoms will be defined
%		   by clicking with the mouse, with the help of a menu. The 
%		   default value for Ti is N/4.
%       DISPLAY    display switch. if DISPLAY=1 a figure is displayed,
%                  otherwise nothing is displayed. default value is 1
%	SIG      : output signal.
%	LOCATOMS : matrix of time-frequency coordinates and durations of the
%		   atoms. 
%
%	Example : 
%	 sig=atoms(128);
%	 sig=atoms(128,[32,0.3,32,1;56,0.15,48,1.22;102,0.41,20,0.7]); 
%	 sig=atoms(128,[32,0.3,32,1;56,0.15,48,1.22;102,0.41,20,0.7],0); 

%	P. Flandrin, May 1995 - O. Lemoine, February 1996.
%	F. Auger - O. Lemoine, June 1996.
%       E. Chassande-Mottin, F. Auger, May 1998.
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

if ( nargin < 1 ),
 error ( 'At least one parameter required' ) ;
end

comp=computer; % so as to know the running computer
MatlabVersion=version; MatlabVersion=str2num(MatlabVersion(1));
%if MatlabVersion<5,
% error('Unfortunately, this version does not run on matlab 4.');
%end;

if (nargin<3), display=1; end;

sig=(1+j)*zeros(N,1);
t=linspace(0,2*pi,100); 
locatoms=[];
Natoms=0; choice=1;
T=N/4; A=1;

if (display==1),
 clf; set(gcf,'Resize','On','NextPlot','Add');

 axsig  = axes('Units','normal','Visible','on','Box','On',...
               'position', [0.10 0.65 0.80 0.25],...
               'XLim', [1 N], 'XGrid','on', ...
               'YLim', [-1 1],'YGrid','on');

 axtfr  = axes('Units','normal','Visible','on','Box','On',...
               'position', [0.10 0.12 0.80 0.45],...
               'XLim',[1 N],  'XGrid','on',...
               'YLim',[0 0.5],'YGrid','on');

 axes(axtfr);
 xlabel('Time'); ylabel('Normalized frequency'); 
 hold on
end;

if (nargin==1),
 fprintf(' Default value for the time-duration : %f\n',T);
 fprintf(' Default value for the amplitude     : %f\n',A);
 while choice~=5,
  MaxSig=max(abs(real(sig)));
  axes(axsig); plot(1:N, real(sig),'g');
  if MaxSig==0.0,
   set(axsig,'XLim',[1 N], 'XGrid','on', 'YLim', [-1 1], 'YGrid','on');
  else
   set(axsig,'XLim',[1 N], 'XGrid','on', 'YLim', [-MaxSig MaxSig],'YGrid','on');
  end;
  title([int2str(Natoms),' Gaussian atom(s)'])

  choice=menu('ATOMS MENU',...
              'Add a gaussian atom',...
              'Delete the last atom',...
              'Change the time-duration',...
              'Change the amplitude',...
              'Stop');
  if choice==1, 
   axes(axtfr); [t0,f0]=ginput(1);             % add a gaussian atom
   t0=round(max(min(t0,N),1)); f0=max(min(f0,0.5),0.0);
   locatoms=[locatoms; t0 f0 T A];
   axes(axtfr); plot(t0,f0,'x'); plot((t0+j*f0)+(0.5*T*cos(t)+j*(2/(T*pi))*sin(t)))
   sig=sig + A*amgauss(N,t0,T) .* fmconst(N,f0,t0);
   Natoms=Natoms+1;
  elseif (choice==2 & Natoms>=1),              % delete last atom
   t0=locatoms(Natoms,1);
   f0=locatoms(Natoms,2);
   Told =locatoms(Natoms,3);
   Aold =locatoms(Natoms,4);
   axes(axtfr); AxtfrChildren=get(gca,'Children'); delete(AxtfrChildren(1:2)); 
   if (Natoms==1)
    Natoms=0; locatoms=[]; sig=(1+j)*zeros(N,1);
   else
    sig=sig - Aold*amgauss(N,t0,Told) .* fmconst(N,f0,t0);
    Natoms=Natoms-1;
    locatoms=locatoms(1:Natoms,:);
   end;
  elseif choice==3,
   fprintf(' Old time duration : %f\n', T);
   Told=T; T=input(' New time duration : ');
   if isempty(T), T=Told; end;
  elseif choice==4,
   fprintf(' Old amplitude : %f\n', A);
   Aold=A; A=input(' New amplitude : ');
   if isempty(A), A=Aold; end;
  end
 end;
elseif (nargin>=2),
 [Natoms,ccoord]=size(coord);
 if (ccoord~=4),
  error('Bad dimension for COORD');
 end;
 for k=1:Natoms,
  t0=round(max(min(coord(k,1),N),1));
  f0=max(min(coord(k,2),0.5),0.0);
  T=coord(k,3); A=coord(k,4);
  if t0~=coord(k,1),
   disp('Warning : ti should be between 1 and N');
  elseif f0~=coord(k,2),
   disp('Warning : fi should be between 0 and 0.5');
  elseif T<0,
   error('T must be positive');
  elseif A<0,
   error('A must be positive');
  else
   sig=sig+A*amgauss(N,t0,T) .* fmconst(N,f0,t0); 
   if (display==1),
    axes(axtfr); plot(t0,f0,'x');
    plot((t0+j*f0)+(0.5*T*cos(t)+j*(2/(T*pi))*sin(t)))
   end;
  end
 end
 locatoms=coord;
end

if (display==1),
 hold off
 MinSig=min(real(sig));
 MaxSig=max(real(sig));
 axes(axsig); plot(1:N, real(sig),'g');
 if MaxSig==0.0,
  set(axsig,'XLim',[1 N], 'XGrid','on', 'YLim', [-1 1], 'YGrid','on');
 else
  set(axsig,'XLim',[1 N], 'XGrid','on', 'YLim', [-MaxSig MaxSig],'YGrid','on');
 end;
 title([int2str(Natoms),' Gaussian atom(s)'])
end;

