//Posted to GitHub, with Creative Commons, CC BY 4.0 license.
//Code release date: May 30, 2026.
//We welcome all constructive suggestions/feedback.

clear;
clc;
close(winsid())

//Parameters in the Woods-Saxon model
a0= 0.429D-15//diffusivity in m
L=2.002D-15 //Width of the potential well in m
V0=-37.199D+6//Depth of the potential well in eV


mp=1.6726*10^-27
mn=1.6749*10^-27
m=(mp*mn)/(mp+mn)//reduced mass of deuteron
hbar=1.054D-34
e=1.602*10^-19

N1=1001
N=N1-2
T=zeros(N,N)
V=zeros(N,N)
a=10D-15
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
    V(i,i)=V0/(1+exp((x(i)-L)/a0))
    V1(i)=V(i,i)
end

H=(((-hbar^2/(2*m*e*h^2))*T)+V)/(10^6) //Hamiltonian
[vec,E]=spec(H)

//For bound states
count=0
for i=1:N
    if E(i,i)<0
        count=count+1
        E1(count)=E(i,i)
        phi(:,count)=vec(:,i)
    end
end

//normalisation of the wavefunctions
for(i=1:count)
    N1=1/sqrt(inttrap(x,phi(:,i).*phi(:,i)))
    psi(:,i)=N1*phi(:,i)
end

//results
disp("no. of bound states = "+string(count))
disp("Binding energy in MeV")
disp(E1)
[val,index]=max(psi(:,1).*psi(:,1))
disp("Size of the nucleus in m")
disp(x(index))

//plots
//Potential plot
scf(0)
clf(0)
plot(x/10^-15,V1/10^6,'LineWidth',4)
xtitle('Woods-Saxon potential','r (fm)','Energy in MeV')
aaa=gca
aaa.data_bounds=[0,-40;5,10]
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
aaa.font_size = 4 

//Wavefunction plot
scf(1)
clf(1)
plot(x/1D-15,psi(:,1)/(3.162D+7),'LineWidth',2)
xlabel("$\text{radial coordinate r (fm)}$", "fontsize", 4)
ylabel("$\text{u(r)}\ (fm )^{(-1/2)}$", "fontsize", 4)
title("u for Woods-Saxon potential")
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
title("Probability Density for Woods-Saxon potential")
aaa=gca
aaa.x_label.font_size = 4
aaa.y_label.font_size = 4
aaa.title.font_size = 4
aaa.font_size=3

