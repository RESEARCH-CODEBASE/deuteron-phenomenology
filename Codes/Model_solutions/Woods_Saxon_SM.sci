//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.

clear
clc
close(winsid())

mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn)//reduced mass
hbar = 1.054D-34
e=1.602D-19

//Parameters in the Woods-Saxon model
a=0.429D-15 //diffusvity
V0=(37.199D+6)*e //Depth of the potential well
R=2.002D-15 //Width of the potential well


pmin=10^-3
pmax=10D-15/a
pmatch=(3/4)*pmax
pf=[pmin:0.001:pmatch]
n1=length(pf)
pb=[pmax:-0.001:pmatch]
n2=length(pb)

tolmin=1D-15
tolfn=1
E=(-5D+6)*e

//Coupled Differential equations
function f1=f_ode(p,u)
    f1=[u(2);((1 -(1/(epsilon*(1+exp(p-p0))))).*u(1))./(k)]
endfunction

function tol=f(E)
    p0=R/a
    k=hbar^2/(2*m*a^2*abs(E))
    epsilon=abs(E)/V0
    uf=ode([1D-3;1],pmin,pf,f_ode)
    ub=ode([1D-5;-9.7D-2],pmax,pb,f_ode)
    ufm=uf(1,n1)
    dufm=uf(2,n1)
    ubm=ub(1,n2)
    dubm=ub(2,n2)
    tol=(ufm*dubm)-(ubm*dufm)
endfunction

countmax=500
count=0

while(abs(tolfn)>tolmin&count<countmax)
    count=count+1
    tolfn=f(E)
    E_new=E-(f(E)/numderivative(f,E))
    E=E_new
    E_iterated(count)=E/((10^6)*e)
    tol_iterated(count)=tolfn   
end

p0=R/a
k=hbar^2/(2*m*a^2*abs(E))
epsilon=abs(E)/V0
uf=ode([1D-3;1],pmin,pf,f_ode)
ub=ode([1D-5;-9.7D-2],pmax,pb,f_ode)

disp("ufmatch: "+string(uf(1,n1)))
disp("ubmatch: "+string(ub(1,n2)))
disp("dufmatch: "+string(uf(2,n1)))
disp("dubmatch: "+string(ub(2,n2)))
[val,index]=max(uf(1,:).*uf(1,:))

//results
disp("Binding energy in MeV")
disp(E/((10^6)*e))
disp("Size of the nucleus in m")
disp(pf(index)*a)

//plots

//Forward solution
scf(1) ;clf(1)
plot(pf,uf(1,:),'r','LineWidth',2)
xtitle('Forward solution','pf','uf')

//Backward solution
scf(2);clf(2)
plot(pb,ub(1,:),'b','LineWidth',2)
xtitle('Backward solution','pb','ub')

//Total wavefunction
scf(3);clf(3)
plot(pf,uf(1,:),'r','LineWidth',2)
plot(pb,ub(1,:),'b','LineWidth',2)
xtitle('Shooting-Matching u for Woods-Saxon','Dimensionless variable p','Wavefunction u')
legend('Forward solution','Backward solution')
aa=gca
aa.x_label.font_size=4
aa.y_label.font_size=4
aa.title.font_size=4
aa.data_bounds=[0,-0.05;pmax,max(uf(1,:))]

//Radial wavefunction
scf(4);clf(4)
plot(pf,uf(1,:)./pf,'r','LineWidth',2)
plot(pb,ub(1,:)./pb,'b','LineWidth',2)
xtitle('Shooting matching R for Woods-Saxon','Dimensionless variable p','Radial wavefunction R')
legend('Forward solution','Backward solution')
aa=gca
aa.data_bounds=[0,-0.1;max(pb),1.1]
aa.x_label.font_size=4
aa.y_label.font_size=4
aa.title.font_size=4

//Derivative of the wavefunction
scf(5);clf(5)
plot(pf,uf(2,:),'r','LineWidth',2)
plot(pb,ub(2,:),'b','LineWidth',2)
xtitle('Shooting matching du/dp','p','du/dp')
legend('Forward solution','Backward solution')
aa=gca
aa.x_label.font_size=4
aa.y_label.font_size=4
aa.title.font_size=4

//Probability density
scf(6);clf(6)
plot(pf,uf(1,:).*uf(1,:),'r','LineWidth',2)
plot(pb,ub(1,:).*ub(1,:),'b','LineWidth',2)
xtitle('Probability density plot','p','Probability density u^2')
legend('Forward solution','Backward solution')
aa=gca
aa.data_bounds=[0,-1;max(pb),10]
aa.x_label.font_size=4
aa.y_label.font_size=4
aa.title.font_size=4
