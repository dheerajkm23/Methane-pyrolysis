function dydt = cca(t,n)  %output=input

N=6.02*10^23*10^-6; %Avogadro constant multiplied with 10^-6 for cm^3 to m^3
pch4=0.45*101325; %pressure of methane in pascal
v0=1*10^-3; %m^3, initial volume
nch4=pch4*v0/(8.314*973); %initial moles of methane
par=0.45*101325;%pressure of argon in pascal
nar=par*v0/(8.314*973);%initial moles of argon
pna=0.1*101325;%pressure of sodium in pascal
nna=pna*v0/(8.314*973);%initial moles of sodium
ni=nch4+nar+nna %total initial moles
nt=n(1)+n(2)+n(3)+n(4)+n(5)+n(6)+n(7)+n(8)+n(9)+n(10)+n(11)+n(12)+n(13)+n(14)+n(15)+n(16)+nar;
% nt is total moles present at any time t.
vt=nt*v0/ni % total volume at any time t since pressure is constant

kf=[1.1*10^-8 2.2*10^9 1.5*10^11 N*2.2*10^-21 N*3.2*10^-22 N*10^-23 N*9.4*10^-27 N*2.3*10^-23 N*3.1*10^-21 N*6.9*10^-13 N*9.9*10^-24 9.5*10^9 8.7*10^10 6.8*10^3 8.4*10^5 7.1*10^3 5*10^9 5.8*10^8 4.4*10^-1 1.4*10^8 2.4*10^1 4.7*10^6 N*1.9*10^-10 N*2.8*10^-10];

kb=[N*4.1*10^-10 N*1.5*10^-9  N*4.4*10^-9 N*1.9*10^-10 4.1*10^9 N*1.9*10^-12 N*6.6*10^-11 3.9*10^2 4.7*10^4 N*5.7*10^-15 N*3.8*10^-19 N*3.3*10^-9 N*7.3*10^-9 N*3*10^-9  N*2.3*10^-10 N*1.4*10^-9 N*1.5*10^-9 N*3.1*10^-9 N*3.7*10^-12 N*5.1*10^-10 N*5*10^-10 N*3.9*10^-11 1.8*10^-10 2.7*10^-4];


%{
ri for each species is defined as d(mole of species)/dt 
e1=ch4 & r1
e2=ch3* & r2
e3=h* & r3
e4=Na2 & r4
e5=Na &r5
e6=Na3   & r6
e7=NaH & r7
e8=HNaCh3 & r8
e9=Na2H & r9
e10=Na2Ch3 & r10
e11=HNa2Ch3 & r11
e12=HNa3Ch3 & r12
e13=h2 & r13
e14=C2H6 & r14
e15=NaCh3 & r15
e16=Na3H & r16
%}

r1=vt*((-kf(1)*(n(1)/vt) + kb(1)*(n(2)/vt)*(n(3)/vt)) + (-kf(4)*(n(1)/vt)*(n(5)/vt) + kb(4)*(n(7)/vt)*(n(2)/vt))  + (-kf(5)*(n(1)/vt)*(n(5)/vt) + kb(5)*(n(8)/vt)) + (-kf(6)*(n(1)/vt)*(n(4)/vt) + kb(6)*(n(9)/vt)*(n(2)/vt)) + (-kf(7)*(n(1)/vt)*(n(4)/vt) + kb(7)*(n(10)/vt)*(n(3)/vt)) + (-kf(8)*(n(1)/vt)*(n(4)/vt) + kb(8)*(n(11)/vt)) + (-kf(9)*(n(1)/vt)*(n(6)/vt) + kb(9)*(n(12)/vt)) + (-kf(10)*(n(1)/vt)*(n(3)/vt) + kb(10)*(n(2)/vt)*(n(13)/vt)) + (-kf(11)*(n(1)/vt)*(n(2)/vt) + kb(11)*(n(14)/vt)*(n(3)/vt)))

r2= vt*((kf(1)*(n(1)/vt) + -kb(1)*(n(2)/vt)*(n(3)/vt)) + (kf(4)*(n(1)/vt)*(n(5)/vt) + -kb(4)*(n(2)/vt)*(n(7)/vt)) + (kf(6)*(n(1)/vt)*(n(4)/vt) + -kb(6)*(n(2)/vt)*(n(9)/vt))+ (kf(10)*(n(1)/vt)*(n(3)/vt) + -kb(10)*(n(2)/vt)*(n(13)/vt))+ (-kf(11)*(n(1)/vt)*(n(2)/vt) + kb(11)*(n(14)/vt)*(n(3)/vt))+ (kf(13)*(n(8)/vt) + -kb(13)*(n(2)/vt)*(n(7)/vt))+ (kf(15)*(n(15)/vt) + -kb(15)*(n(2)/vt)*(n(5)/vt))+ (kf(18)*(n(10)/vt) + -kb(18)*(n(2)/vt)*(n(4)/vt))+ 2*(-kf(24)*(n(2)/vt)*(n(2)/vt) + kb(24)*(n(14)/vt)))

r3= vt*((kf(1)*(n(1)/vt) + -kb(1)*(n(2)/vt)*(n(3)/vt))+ (kf(7)*(n(1)/vt)*(n(4)/vt) + -kb(7)*(n(3)/vt)*(n(10)/vt))+ (-kf(10)*(n(3)/vt)*(n(1)/vt) + kb(10)*(n(2)/vt)*(n(13)/vt))+ (kf(11)*(n(1)/vt)*(n(2)/vt) + -kb(11)*(n(3)/vt)*(n(14)/vt))+ (kf(12)*(n(8)/vt) + -kb(12)*(n(3)/vt)*(n(15)/vt))+ (kf(14)*(n(7)/vt) + -kb(14)*(n(3)/vt)*(n(5)/vt))+ (kf(16)*(n(9)/vt) + -kb(16)*(n(3)/vt)*(n(4)/vt))+ (kf(21)*(n(16)/vt) + -kb(21)*(n(3)/vt)*(n(6)/vt))+ 2*(-kf(23)*(n(3)/vt)*(n(3)/vt) + kb(23)*(n(13)/vt)))

r4= vt*((-kf(2)*(n(4)/vt) + kb(2)*(n(5)/vt)*(n(5)/vt))+ (kf(3)*(n(6)/vt) + -kb(3)*(n(4)/vt)*(n(5)/vt))+ (-kf(6)*(n(1)/vt)*(n(4)/vt) + kb(6)*(n(9)/vt)*(n(2)/vt))+ (-kf(7)*(n(1)/vt)*(n(4)/vt) + kb(7)*(n(10)/vt)*(n(3)/vt))+ (-kf(8)*(n(1)/vt)*(n(4)/vt) + kb(8)*(n(11)/vt))+ (kf(16)*(n(9)/vt) + -kb(16)*(n(4)/vt)*(n(3)/vt))+ (kf(18)*(n(10)/vt) + -kb(18)*(n(4)/vt)*(n(2)/vt))+ (kf(20)*(n(16)/vt) + -kb(20)*(n(4)/vt)*(n(7)/vt)))

r5= vt*(2*(kf(2)*(n(4)/vt) + -kb(2)*(n(5)/vt)*(n(5)/vt))+ (kf(3)*(n(6)/vt) + -kb(3)*(n(4)/vt)*(n(5)/vt))+ (-kf(4)*(n(5)/vt)*(n(1)/vt) + kb(4)*(n(2)/vt)*(n(7)/vt))+ (-kf(5)*(n(5)/vt)*(n(1)/vt) + kb(5)*(n(8)/vt))+ (kf(14)*(n(7)/vt) + -kb(14)*(n(5)/vt)*(n(3)/vt))+ (kf(15)*(n(15)/vt) + -kb(15)*(n(5)/vt)*(n(2)/vt))+ (kf(17)*(n(9)/vt) + -kb(17)*(n(5)/vt)*(n(7)/vt))+ (kf(19)*(n(10)/vt) + -kb(19)*(n(5)/vt)*(n(15)/vt))+ (kf(22)*(n(16)/vt) + -kb(22)*(n(5)/vt)*(n(9)/vt)))

r6= vt*((-kf(3)*(n(6)/vt) + kb(3)*(n(4)/vt)*(n(5)/vt))+ (-kf(9)*(n(6)/vt)*(n(1)/vt) + kb(9)*(n(12)/vt))+ (kf(21)*(n(16)/vt) + -kb(21)*(n(6)/vt)*(n(3)/vt)))

r7= vt*((kf(4)*(n(1)/vt)*(n(5)/vt) + -kb(4)*(n(2)/vt)*(n(7)/vt))+ (kf(13)*(n(8)/vt) + -kb(13)*(n(2)/vt)*(n(7)/vt))+ (-kf(14)*(n(7)/vt) + kb(14)*(n(3)/vt)*(n(5)/vt))+ (kf(17)*(n(9)/vt) + -kb(17)*(n(5)/vt)*(n(7)/vt))+ (kf(20)*(n(16)/vt) + -kb(20)*(n(4)/vt)*(n(7)/vt)))

r8= vt*((kf(5)*(n(5)/vt)*(n(1)/vt) + -kb(5)*(n(8)/vt))+ (-kf(12)*(n(8)/vt) + kb(12)*(n(15)/vt)*(n(3)/vt))+ (-kf(13)*(n(8)/vt) + kb(13)*(n(2)/vt)*(n(7)/vt)))

r9= vt*((kf(6)*(n(4)/vt)*(n(1)/vt) + -kb(6)*(n(9)/vt)*(n(2)/vt))+ (-kf(16)*(n(9)/vt) + kb(16)*(n(4)/vt)*(n(3)/vt))+ (-kf(17)*(n(9)/vt) + kb(17)*(n(5)/vt)*(n(7)/vt))+ (kf(22)*(n(16)/vt) + -kb(22)*(n(9)/vt)*(n(5)/vt)))

r10=  vt*((kf(7)*(n(1)/vt)*(n(4)/vt) + -kb(7)*(n(10)/vt)*(n(3)/vt))+ (-kf(18)*(n(10)/vt) + kb(18)*(n(2)/vt)*(n(4)/vt))+ (-kf(19)*(n(10)/vt) + kb(19)*(n(5)/vt)*(n(15)/vt)))

r11=  vt*((kf(8)*(n(1)/vt)*(n(4)/vt) + -kb(8)*(n(11)/vt)))

r12= vt*((kf(9)*(n(1)/vt)*(n(6)/vt) + -kb(9)*(n(12)/vt)))

r13= vt*((kf(10)*(n(1)/vt)*(n(3)/vt) + -kb(10)*(n(2)/vt)*(n(13)/vt))+ (kf(23)*(n(3)/vt)*(n(3)/vt) + -kb(23)*(n(13)/vt)))

r14=  vt*((kf(11)*(n(1)/vt)*(n(2)/vt) + -kb(11)*(n(14)/vt)*(n(3)/vt))+ (kf(24)*(n(2)/vt)*(n(2)/vt) + -kb(24)*(n(14)/vt)))

r15=  vt*((kf(12)*(n(8)/vt) + -kb(12)*(n(15)/vt)*(n(3)/vt))+ (-kf(15)*(n(15)/vt) + kb(15)*(n(2)/vt)*(n(5)/vt))+ (kf(19)*(n(10)/vt) + -kb(19)*(n(5)/vt)*(n(15)/vt)))

r16=  vt*((-kf(20)*(n(16)/vt) + kb(20)*(n(4)/vt)*(n(7)/vt))+ (-kf(21)*(n(16)/vt) + kb(21)*(n(6)/vt)*(n(3)/vt))+ (-kf(22)*(n(16)/vt) + kb(22)*(n(9)/vt)*(n(5)/vt)))


dydt=[r1; r2; r3; r4; r5; r6; r7; r8; r9; r10; r11; r12; r13; r14; r15; r16]

end








