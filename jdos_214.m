n_q = 401;
E = 0.001;
n_E = 11;
[qx,qy] = meshgrid(linspace(-1,1,n_q));
Epoints = linspace(-E,E,n_E);
I = eye(4);


mu_0 = 1;
t_h = 0.1;
V_m = 0.1;
mu_z = 0.7;
t_z = 0.55;
t_z1 = 0.2;
t = 1;

epsilon_xz = - mu_0 - 2*t*cos(pi*qx) - 2*t_h*cos(pi*qy);
epsilon_yz = - mu_0 - 2*t*cos(pi*qy) - 2*t_h*cos(pi*qx);
V_hyb = - 2*V_m*sin(pi*qx).*sin(pi*qy);
epsilon_alpha = 1/2*((epsilon_xz + epsilon_yz) ...
    - sqrt((epsilon_xz - epsilon_yz).^2 + 4*V_hyb.^2));
epsilon_beta = 1/2*((epsilon_xz + epsilon_yz) ...
    + sqrt((epsilon_xz - epsilon_yz).^2 + 4*V_hyb.^2));

D0 = 0.00035;
D = D0*(cos(pi*qx) - cos(pi*qy));

V_intra = 1;
V_inter = 0.1;
Vs = [V_intra,0,V_inter,0;...
    0,-V_intra,0,-V_inter;...
    V_inter,0,V_intra,0;...
    0,-V_inter,0,-V_intra];

d = 0.0015;
G11 = zeros(n_q,n_q,n_E);
G12 = zeros(n_q,n_q,n_E);
G13 = zeros(n_q,n_q,n_E);
G14 = zeros(n_q,n_q,n_E);
G21 = zeros(n_q,n_q,n_E);
G22 = zeros(n_q,n_q,n_E);
G23 = zeros(n_q,n_q,n_E);
G24 = zeros(n_q,n_q,n_E);
G31 = zeros(n_q,n_q,n_E);
G32 = zeros(n_q,n_q,n_E);
G33 = zeros(n_q,n_q,n_E);
G34 = zeros(n_q,n_q,n_E);
G41 = zeros(n_q,n_q,n_E);
G42 = zeros(n_q,n_q,n_E);
G43 = zeros(n_q,n_q,n_E);
G44 = zeros(n_q,n_q,n_E);
G11r = zeros(n_q,n_q,n_E);
G12r = zeros(n_q,n_q,n_E);
G13r = zeros(n_q,n_q,n_E);
G14r = zeros(n_q,n_q,n_E);
G21r = zeros(n_q,n_q,n_E);
G22r = zeros(n_q,n_q,n_E);
G23r = zeros(n_q,n_q,n_E);
G24r = zeros(n_q,n_q,n_E);
G31r = zeros(n_q,n_q,n_E);
G32r = zeros(n_q,n_q,n_E);
G33r = zeros(n_q,n_q,n_E);
G34r = zeros(n_q,n_q,n_E);
G41r = zeros(n_q,n_q,n_E);
G42r = zeros(n_q,n_q,n_E);
G43r = zeros(n_q,n_q,n_E);
G44r = zeros(n_q,n_q,n_E);
for k=1:n_E
for i=1:n_q
for j=1:n_q
   H = [epsilon_alpha(i,j),D(i,j),0,0;...
     D(i,j),-epsilon_alpha(i,j),0,0;...
     0,0,epsilon_beta(i,j),D(i,j);...
     0,0,D(i,j),-epsilon_beta(i,j)];
   G0 = I/((Epoints(k) + 1i*d)*I - H);
   G11(i,j,k) = G0(1,1);
   G12(i,j,k) = G0(1,2);
   G13(i,j,k) = G0(1,3);
   G14(i,j,k) = G0(1,4);
   G21(i,j,k) = G0(2,1);
   G22(i,j,k) = G0(2,2);
   G23(i,j,k) = G0(2,3);
   G24(i,j,k) = G0(2,4);
   G31(i,j,k) = G0(3,1);
   G32(i,j,k) = G0(3,2);
   G33(i,j,k) = G0(3,3);
   G34(i,j,k) = G0(3,4);
   G41(i,j,k) = G0(4,1);
   G42(i,j,k) = G0(4,2);
   G43(i,j,k) = G0(4,3);
   G44(i,j,k) = G0(4,4);
end
end
   G11r(:,:,k) = fftshift(fft2(G11(:,:,k)))/(n_q*n_q);
   G12r(:,:,k) = fftshift(fft2(G12(:,:,k)))/(n_q*n_q);
   G13r(:,:,k) = fftshift(fft2(G13(:,:,k)))/(n_q*n_q);
   G14r(:,:,k) = fftshift(fft2(G14(:,:,k)))/(n_q*n_q);
   G21r(:,:,k) = fftshift(fft2(G21(:,:,k)))/(n_q*n_q);
   G22r(:,:,k) = fftshift(fft2(G22(:,:,k)))/(n_q*n_q);
   G23r(:,:,k) = fftshift(fft2(G23(:,:,k)))/(n_q*n_q);
   G24r(:,:,k) = fftshift(fft2(G24(:,:,k)))/(n_q*n_q);
   G31r(:,:,k) = fftshift(fft2(G31(:,:,k)))/(n_q*n_q);
   G32r(:,:,k) = fftshift(fft2(G32(:,:,k)))/(n_q*n_q);
   G33r(:,:,k) = fftshift(fft2(G33(:,:,k)))/(n_q*n_q);
   G34r(:,:,k) = fftshift(fft2(G34(:,:,k)))/(n_q*n_q);
   G41r(:,:,k) = fftshift(fft2(G41(:,:,k)))/(n_q*n_q);
   G42r(:,:,k) = fftshift(fft2(G42(:,:,k)))/(n_q*n_q);
   G43r(:,:,k) = fftshift(fft2(G43(:,:,k)))/(n_q*n_q);
   G44r(:,:,k) = fftshift(fft2(G44(:,:,k)))/(n_q*n_q);
end
Ak = -imag(G11+G22+G33+G44)/pi;
imagesc(Ak(:,:,n_E))

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
for k=1:1
    dnqx(:,:,k)=fftshift(dnq(:,:,k));
end

dnq1 = zeros(2*n_q+1,2*n_q+1,n_E);
for k=1:1
dnq_tem=zeros(3*n_q,3*n_q);
dnq_tem(1:n_q,1:n_q)=abs(dnqx(:,:,k));
dnq_tem(1:n_q,n_q+1:2*n_q)=abs(dnqx(:,:,k));
dnq_tem(1:n_q,2*n_q+1:3*n_q)=abs(dnqx(:,:,k));
dnq_tem(n_q+1:2*n_q,1:n_q)=abs(dnqx(:,:,k));
dnq_tem(n_q+1:2*n_q,n_q+1:2*n_q)=abs(dnqx(:,:,k));
dnq_tem(n_q+1:2*n_q,2*n_q+1:3*n_q)=abs(dnqx(:,:,k));
dnq_tem(2*n_q+1:3*n_q,1:n_q)=abs(dnqx(:,:,k));
dnq_tem(2*n_q+1:3*n_q,n_q+1:2*n_q)=abs(dnqx(:,:,k));
dnq_tem(2*n_q+1:3*n_q,2*n_q+1:3*n_q)=abs(dnqx(:,:,k));
dnq1(:,:,k)=dnq_tem((n_q+1)/2:(n_q+1)/2+2*n_q,(n_q+1)/2:(n_q+1)/2+2*n_q);
end
