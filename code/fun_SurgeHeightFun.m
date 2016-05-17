function [surgeheight]=fun_SurgeHeightFun(h0,u10,m,check_plot)
%  
% Code to estimate storm surge by using one dimensional sureg model       %
% The storm surge model explained in the book titled:                     %
% "water wave mechanics for engineers and scientists"                     %
% by Dean Dalrymple method 1992                                           %
% Ver 1                                                                   %
%                                                     by: Arash Karimpour %
%                                                   www.arashkarimpour.com%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% 
%INPUT---------------------------------------------------------------------
% Modeling the storm surge using Dean Dalrymple method 1992
% 
% h0=40;               % deep water depth in (m), It was 42 m for Katrina
% u10=40;              % wind velocity 10 meter above the ground in (m/s), It was 40 m/s for Katrina 
% m=0.00084;           % bed slope = h0/l; 
% 
%OUTPUT--------------------------------------------------------------------
% 
% surgeheight         : Surge Height in (m)
% 
%--------------------------------------------------------------------------
% 
% Borja G. Reguero 20160406 adapted for climada coastal module 
%-- 
%%
%--------------------------------------------------------------------------
% initial values

n=1.3; %shear stress coefficient
Row=1028; %Seawater density in kg/m3
Roa=1.204; %air density in kg/m3

if m==0
    l=50000; % length of the continental shelf, It was 50000 m for Katrina
else
    l=h0/m;
end

%--------------------------------------------------------------------------

x(:,1)=linspace(0,l,1000); %x axis points

%using Wu (1982), method used by ADCIRC
k=(.80+.065*40)/1000;
k(k>0.0034)=0.0034;
Tauwind=Roa*k*u10.^2;
A=n*Tauwind*l/(Row*9.81*h0^2);

% surge calculation

if m==0 % for bed without slope
    
    %claculating eta for bed with zero slope
    eta=h0.*(sqrt(1+2*A.*x/l)-1);
    
else %for bed with slpoe
    
    h=-m.*x+h0; % h=h0.*(1-x./l);
    
    %claculating eta for bed with slope
    eta(1,1)=0;
%     tic
    for i=2:length(x(:,1))
        etaini=eta(i-1); %initial value for eta from previous step
        f=@(eta1)(x(i,1)/l-((1-(h(i,1)+eta1)/h0)-A*log(((h(i,1)+eta1)/h0-A)/(1-A)))); %x/l=(1-(h+eta)/h0)-A*log(((h+eta)/h0-A)/(1-A))
        eta(i,1)=fzero(f,etaini);
    end
%     toc
end

surgeheight=eta(end,1);

%--------------------------------------------------------------------------
if check_plot 
%plotting
close
figure 

if m==0
    x1=[0;l+l/10;l+l/10;l;l-l/10;0];
    y1=[-h0/10;-h0/10;h0-h0/10;h0-h0/10;0;0];
    fill(x1,y1,1/255*[102,51,0])
    
    x2=[0;l-l/10;l;l+l/10;l+l/10;0];
    y2=[0;0;h0-h0/10;h0-h0/10;h0;h0];
    hold on
    fill(x2,y2,1/255*[153,204,255])
    
    x33=flipud(x);
    x3=[x;x33];
    y331(1:length(x(:,1)),1)=h0;
    y332=flipud(eta)+h0;
    y3=[y331;y332];
    fill(x3,y3,1/255*[255,0,0])
    
    xlim([0 l+l/10])
    
    hleg1=legend('Sea Bed','Water','Storm Surge');
    set(hleg1,'Location','SouthEast')    
    
else
    
    x1=[0;l+l/10;l+l/10;l;0];
    y1=[-h0/10;-h0/10;h0-h0/10;h0-h0/10;0];
    fill(x1,y1,1/255*[102,51,0])
    
    x2=[0;l;l+l/10;l+l/10;0];
    y2=[0;h0-h0/10;h0-h0/10;h0;h0];
    hold on
    fill(x2,y2,1/255*[153,204,255])
    
    x33=flipud(x);
    x3=[x;x33];
    y331(1:length(x(:,1)),1)=h0;
    y332=flipud(eta)+h0;
    y3=[y331;y332];
    fill(x3,y3,1/255*[255,0,0])
    
    xlim([0 l+l/10])
    %ylim([-h0/10 ;:])
    
    hleg1=legend('Sea Bed','Water','Storm Surge');
    set(hleg1,'Location','SouthEast')    

end

end
%--------------------------------------------------------------------------
% no need to clear workspace because is a function now 

% % % clear h0 l u10 n Row Roa m
% % % clear k Tauwind A h
% % % clear etaini i f x eta
% % % clear x1 y1 x2 y2 x3 x33 y3 y33 y331 y332