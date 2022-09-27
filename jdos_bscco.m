n_q = 401;
E = 0.024;
n_E = 16;
[qx,qy] = meshgrid(linspace(-1,1,n_q));
Epoints = linspace(-E,E,n_E);
I = eye(2);
Vs = 0.1;
Vm = 0.0;


t1 = 0.1305; 
t2 = -0.5951;
t3 = 0.1636;
t4 = -0.0519;
t5 = -0.1117;
t6 = 0.0510;
E_tb = t1 + t2*(cos(pi*qx)+cos(pi*qy))/2 ...
    +t3*cos(pi*qx).*cos(pi*qy)+ t4*(cos(2*pi*qx)+cos(2*pi*qy))/2 ...
    +t5*(cos(2*pi*qx).*cos(pi*qy)+cos(pi*qx).*cos(2*pi*qy))/2 ...
    + t6*cos(2*pi*qx).*cos(2*pi*qy);

D0 = 0.025;
D = D0/2*(cos(pi*qx)-cos(pi*qy));

V = [Vs + Vm,0;0,-Vs + Vm];

d = 0.0015;

G11 = zeros(n_q,n_q,n_E);
G12 = zeros(n_q,n_q,n_E);
G21 = zeros(n_q,n_q,n_E);
G22 = zeros(n_q,n_q,n_E);
G11r = zeros(n_q,n_q,n_E);
G12r = zeros(n_q,n_q,n_E);
G21r = zeros(n_q,n_q,n_E);
G22r= zeros(n_q,n_q,n_E);
for k=1:n_E
for i=1:n_q
for j=1:n_q
   H = [E_tb(i,j),D(i,j);D(i,j),-E_tb(i,j)];
   G0 = I/((Epoints(k) + 1i*d)*I - H);
   G11(i,j,k) = G0(1,1);
   G12(i,j,k) = G0(1,2);
   G21(i,j,k) = G0(2,1);
   G22(i,j,k) = G0(2,2);
end
end
end

Ak = -imag(G11+G22)/pi;

dnq = zeros(n_q,n_q,n_E);
for k=1:n_E
for m=1:n_q
for l=1:n_q
    sum1 = 0;
for i=1:n_q
for j=1:n_q
    ii = i + m;
    jj = j + l;
    if i + m > n_q
        ii = i + m - n_q ;
    end
    if j + l > n_q
        jj = j + l - n_q ;
    end
    if i + m < 1
        ii = i + m + n_q ;
    end
    if j + l < 1
        jj = j + l + n_q ;
    end
    sum1 = sum1 + Ak(i,j,k)*Ak(ii,jj,k);
end
end
    dnq(m,l,k) = sum1;
end
end
end
dnqx=zeros(n_q,n_q,n_E);
for k=1:n_E
    dnqx(:,:,k)=fftshift(dnq(:,:,k));
end