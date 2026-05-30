//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.

clear;clc;
close(winsid())

//Parameters in Morse potential
am=3.94E-16
rm=12.43E-16
V0=9.4525E+07

mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn)//reduced mass
e=1.602*10^-19
hbar=1.0545*(10^-34)

N1=1001
N=N1-2
T=zeros(N,N)
V=zeros(N,N)

a=10*10^-15
h=a/(N+1)
x(1)=h

//For Hamiltonian matrix
for i=1:N
    T(i,i)=-2
    if(i<N)
        x(i+1)=x(i)+h
        T(i,i+1)=1
        T(i+1,i)=1
    end
    V(i,i)=V0*(exp(-2*(x(i)-rm)/am)-(2*exp(-(x(i)-rm)/am)))
    V1(i)=V(i,i) 
end

H=((-hbar^2/(2*m*e*h^2))*T)+V
[vec,E]=spec(H)

//Bound states
count=0
for i=1:N
    if E(i,i)<0
        count=count+1
        E1(count)=E(i,i)
        phi(:,count)=vec(:,i)
    end
end

//For normalisation of the wavefunctions
for(i=1:count)
    N1=1/sqrt(inttrap(x,phi(:,i).*phi(:,i)))
    psi(:,i)=N1*phi(:,i)
end

//results
disp("no. of bound states = "+string(count))
disp("Bound state energies in MeV")
disp(E1/(10^6))
[val,index]=max((psi(:,1).*psi(:,1)))
disp("The most probable x in m")
disp(x(index))


//For the plots
//Potential plot
scf(0)
plot(x/10^-15,V1/10^6,'LineWidth',4)
xtitle('Morse potential','r (fm)','Energy in MeV')
aa=gca
aa.data_bounds=[0,-150;5,200]
aa.x_label.font_size = 4
aa.y_label.font_size = 4
aa.title.font_size = 4
aa.font_size=4

//Wavefunction plot
scf(1)
clf(1)
plot(x/1D-15,psi(:,1)/(3.162D+7),'LineWidth',2)
xlabel("$\text{radial coordinate r (fm)}$", "fontsize", 4)
ylabel("$\text{u(r)}\ (fm )^{(-1/2)}$", "fontsize", 4)
title("u for Morse potential")
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
aaa.font_size=3

//Probability density plot
scf(2)
clf(2)
plot(x/1D-15,psi(:,1).*psi(:,1)/(1D+15),'LineWidth',2)
xlabel("$\text{radial coordinate r (fm)}$", "fontsize", 4)
ylabel("$\text{u^2(r)}\ (fm )^{(-1)}$", "fontsize", 4)
title("Probability Density for Morse potential")
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
aaa.font_size=3
