%% Plotting a family of curves for bending of Cosserat circular bars
% This family is for different values of characteristic length of bending.
% The same procedure can be carried for whatever variable you wish to
% govern the family of curves. Just change the names of variables
% accordingly.  

clear,clc, close all

%% These are the stubby vert structures
diam = [7.78 18.95 37.76 56.95 75.86]';     %Enter the diameter of the specimens from smallest to largest in mm
E_measured = [92.4 20.7 8.51 8.38 4.72]';   %Enter the BVS bending modulus in MPa of the specimens from smallest to largest specimen.

%% Known Parameters
Ec = 3.14;                          %Asymptotic modulus (MPa)
Output = E_measured./Ec;            %Relative Stiffness
v = 0.05;               %For plotting curves of different lb, Poisson's ratio was required to stay constant 
                                    %which is why it is set here.
lbend = [4.4 8.76 13];              %This is the vector of lb values of which I want curves. More or less can be calculated
                                    %but the number of data sets in the data sets section and the number of curves being 
                                    %plotted in the plotting section must
                                    %match.
N = 0.99;               %For plotting curves of different lb, N was required to stay constant which is why it is set here.
p = 0.5;                            %beta/gamma set to be the same for all curves
error = [2.97 2.27 1.57 0.87 0.17]';    %This is for calculating error bars.  Error on each side of the data point. Error was symmetrical in this case.
errorOut = error./Ec;               %Changing error into necessary scale for plotting

%% Cosserat Parameters
del = 1./lbend;     %Calculating parameters necessary for use of the fitting equation below.


%% Plotting Equation
fitFunc = @(del,x) (1+(8.*N^2./(v+1)).*((1-p.^2)./(del.*x./2).^2+((p+v).^2./((8.*N.^2.*(1-v)) ...
            +((del.*x./2).^2.*(del.*(x./2).*besseli(0,(del.*x./2))-besseli(1,(del.*x./2)))./ ...
            (del.*x./2.*besseli(0,(del.*x./2))-2.*besseli(1,(del.*x./2))))))));
        
%% Making Data Sets
%Add or subtract data sets in the same notation as shown here to make the
%number of desired curves.  Each set corresponds to a specific element in a
%vector of values for which you want to plot. Follow the variables names.

xx = (1:0.1:100)';
yy0 = ones(size(xx,1));
yy1 = fitFunc(del(1),xx);
yy2 = fitFunc(del(2),xx);
yy3 = fitFunc(del(3),xx);

%% Plotting
%Add or subtract plotting lines to achieve the number of desired curves.
figure 
hold on
a1 = plot(diam,Output,'k.','MarkerSize',30);
a2 = plot(xx,yy0,'g--','LineWidth',2);
a3 = plot(xx,yy1,'b:','LineWidth',2);
a4 = plot(xx,yy2,'k-','LineWidth',2);
a5 = plot(xx,yy3,'c-.','LineWidth',2);
errorbar(diam,Output,errorOut,'r.','LineWidth',2);
% legend('Experimental Ponits','Classical Elasticity','l_t = 5.0 mm','l_t = 9.4 mm','l_t = 14 mm');
xlabel('Diameter (mm)','FontSize',16)
ylabel('Relative Stiffness \Omega','FontSize',16)
xlim([0 100])
ylim([0 40])
set(gca,'FontSize',14);
hold off