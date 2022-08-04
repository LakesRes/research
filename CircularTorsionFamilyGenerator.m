%% Plotting a family of curves for Torsion of Cosserat Circular bars
% This family is for different values of characteristic length of torsion.
% The same procedure can be carried for whatever variable you wish to
% govern the family of curves. Just change the names of variables
% accordingly.  

clear,clc, close all

%% These are the stubby vert structures
diam = [7.78 18.95 37.76 56.95 75.89]';     %Enter the diameter of the specimens from smallest to largest
Measured_G = [39.6 6.49 2.97 2.18 1.8]';    %Enter the BVS shear modulus of the specimens from smallest to largest specimen.

%% Known Parameters - Change as necessary
var = 1.1;                          %Asymptotic shear modulus (MPa)
Output = Measured_G./var;           %Relative Stiffness
psi = 1.0;                          %Psi values
ltors = [5.0 9.4 14];               %This is the vector of lt values of which I want curves. More or less can be calculated
                                    %but the number of data sets in the data sets section and the number of curves being 
                                    %plotted in the plotting section must
                                    %match.
N = 0.99983;          %For plotting curves of different lt, N was required to stay constant which is why it is set here.

Spread = [1.22 0.93 0.63 0.335 0.04]';      %This is for calculating error bars.  Error on each side of the data point. Error was symmetrical in this case.
OutSpread = Spread./var;                    %Changing error into necessary scale for plotting

%% Cosserat Parameters
bplusg = (ltors./1000).^2.*2.*var.*10^6;
alpha = (bplusg./1.5)-bplusg;
kappa = (-2.*var.*(10^6)*N.^2)./(N.^2-1);
p = -((2.*kappa)./(alpha+bplusg)).^(1/2);

%% Plotting Equation
fitFunc = @(lt,p,x) 1+6*((lt/1000)./(x./2000)).^2.*((1-4/3*psi.* ...
            (besseli(1,(p.*x./2000))./((p.*x./2000).*besseli(0,(p.* ...
            x./2000)))))./(1-psi.*(besseli(1,(p.*x./2000))./((p.*x./ ...
            2000).*besseli(0,p.*x./2000)))));
        
%% Making Data Sets
%Add or subtract data sets in the same notation as shown here to make the
%number of desired curves.  Each set corresponds to a specific element in a
%vector of values for which you want to plot. Follow the variables names.
xx = (1:0.1:100)';
yy0 = ones(size(xx,1));
yy1 = fitFunc(ltors(1),p(1),xx);
yy2 = fitFunc(ltors(2),p(2),xx);
yy3 = fitFunc(ltors(3),p(3),xx);

%% Plotting
%Add or subtract plotting lines to achieve the number of desired curves.
figure 
hold on
a1 = plot(diam,Output,'k.','MarkerSize',30);
a2 = plot(xx,yy0,'g--','LineWidth',2);
a3 = plot(xx,yy1,'b:','LineWidth',2);
a4 = plot(xx,yy2,'k-','LineWidth',2);
a5 = plot(xx,yy3,'c-.','LineWidth',2);
errorbar(diam,Output,OutSpread,'r.','LineWidth',2);
% legend('Experimental Ponits','Classical Elasticity','l_t = 5.0 mm','l_t = 9.4 mm','l_t = 14 mm');
xlabel('Diameter (mm)','FontSize',16)
ylabel('Relative Stiffness \Omega','FontSize',16)
xlim([0 100])
ylim([0 40])
set(gca,'FontSize',14);
hold off