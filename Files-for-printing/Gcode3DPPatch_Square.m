%Description: This script generates a txt file that can be used as gcode
%to generate the route of the nozzle to print a sinusoidal fiber patch.
%Written by Diego Quevedo-Moreno (07/08/2020)
%Modified by Claudia Varela (2/23/2021)

%% G-code generation
clc, clear all, close all
yv=linspace(0,10*pi,500);
Ns=1;%Number of Samples
A=0;%Amplitude(mm)
P=1;%Density 
WL=5;%Wave Length(2*pi/WL=wavelength)(mm) 
W=.150;%Fiber Width(mm)
SL=7;%Sample Length(mm)
SW=7;%Sample Width(mm)
nL=5;%#Layers
Lt=.020;%Layer thickness(mm)
iL=0;%Initial Layer Thickness(mm)
dp=.5;%Distance Between points
F= 275;%Printing Speed(mm/min)
nf=(SW/W)*P;%#Fibers
dL=0.0002; %Layer height decrease between layers
%df=(SW-(2*A)-W)/(nf-1);%Distance between fibers
df=(SW-W)/(nf-1);%Distance between fibers version 2 
Offset=[-20 0];%%Offset if needed [X Y] (mm)-20 normally in Y
x=0:dp:SL;
Vy=0:0.01:(SW-(2*A));

%Sinusoidal Generator
for j=1:nf
    for i=1:length(x)
            y(i,j)=A*sin(x(i)*(pi)/WL)+(j-1)*df;
    end
end
%Sinusoidal Pattern Generator
for j=0:nf-1
    if rem(j,2)==0 %Odd direction 
        for i=1:length(x)
            xt(i+(length(x)*j))=x(i);
            yt(i+(length(x)*j))=y(i,(j+1));
        end
    else %Even direction
        for i=1:length(x)
            xt(i+(length(x)*j))= x((SL/dp+1)-i+1);
            yt(i+(length(x)*j))=y(((SL/dp+1)-i+1),(j+1));
        end
    end
end

%%Skirt 
svx=[-4 -4 -3 -3 -2 -2 -1 -1];
svy=[0 SW SW 0 0 SW SW 0];
svz=[iL iL iL iL iL iL iL iL];

xtf=flip(xt);
ytf=flip(yt);
for k=1:nL
    if rem(k,2)==1
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xt(i);
            gy(i+length(xt)*(k-1))=yt(i);
            gz(i+length(xt)*(k-1))=((k-1)*Lt)+iL-(dL*(k-1));
        end
    else
        for i=1:length(xt)
            gx(i+length(xt)*(k-1))=xtf(i);
            gy(i+length(xt)*(k-1))=ytf(i);
            gz(i+length(xt)*(k-1))=((k-1)*Lt)+iL-(dL*(k-1));
        end
    end
end

%Grip Zone
gzvx=0:W/2:10;
for k=1:2:length(gzvx)
gzx(k)=gzvx(k);
gzy(k)=SW+A;
gzz(k)=iL+Lt*nL;
if k>1
gzx(k-1)=gzvx(k);
gzy(k-1)=0-A;
gzz(k-1)=iL+Lt*nL;
end
end

%Vector with Grip Zone 
% gx=[svx gx gzx gzx+SL-10];
% gy=[svy gy gzy gzy];
% gz=[svz gz gzz gzz];

% %Vector without Grip Zone 
gx=[svx gx];
gy=[svy gy];
gz=[svz gz];

%Vector without Grip Zone and Skirt
% gx=[gx];
% gy=[gy];
% gz=[gz];


 
%%Offset (if needed)
gx=gx + Offset(1);
gy=gy + Offset(2);



for j=1:Ns
    for i=1:length(gx)
        gtx(i+(length(gx)*(j-1)))=gx(i);
        gty(i+(length(gx)*(j-1)))=(gy(i)+(SW+3)*(j-1));
        gtz(i+(length(gx)*(j-1)))=gz(i);
    end 
end
figure 
axis equal 
grid off
axis on
hold on
plot3(gtx,gty,gtz,'-','LineWidth',2);
% set(gca,'Color',[217 217 217]/255)
xlabel('mm');
ylabel('mm');
zlabel('mm');
box on


%% File saving
% Nomenclature
% p(density)_a(amplitude)_size_F(printing speed)_dp(distance between points)_L(# of layers)_#samples_fw(fiber width)_Lt(Layer Thickness)
currentFolder = pwd;
path=strcat(currentFolder,'\SampleTest_p100_a0_7x7mm_F275_dp.5mm_L5_1_fw150_Lt.02mm.gcode');
fid = fopen(path,'w');
fprintf(fid,'F%8.6f\r\n',F);

for j = 1:length(gtx)

    fprintf(fid,'LINEAR X%8.6f Y%8.6f Z%8.6f\r\n',gtx(j),gty(j),gtz(j));
    
end
    
fclose('all');
