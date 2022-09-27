n_q = 201; % number of points in k-space
E = 0.2; % energy range in eV(from -1meV to 1meV)
n_E = 9; % number of energy layers
[qx,qy] = meshgrid(linspace(-1,1,n_q));
Epoints = linspace(-E,E,n_E);
I = eye(4);

%the parameters
mu = 1.6;%%
t1 = -1;
t2 = 1.3;
t3 = -0.85;
t4 = -0.85;

%the tight binding model to fit superconducting STM spectra
epsilon_x = -2*t1*cos(pi*qx) - 2*t2*cos(pi*qy) - 4*t3*cos(pi*qx).*cos(pi*qy);
epsilon_y = -2*t2*cos(pi*qx) - 2*t1*cos(pi*qy) - 4*t3*cos(pi*qx).*cos(pi*qy);
epsilon_xy = -4*t4*sin(pi*qx).*sin(pi*qy);

%the energy gap on α : β Fermi surface
D0 = 0.1;
D1 = D0*sin(pi*qx).*sin(pi*qy);
D2 = D0*sin(pi*qx).*sin(pi*qy);

%the scattering potential
V0 = 4*D0;
V_intra = [V0,0,0,0;...
    0,-V0,0,0;...
    0,0,V0,0;...
    0,0,0,-V0];

d = 0.005;
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

   H = [epsilon_x(i,j)-mu,D1(i,j),epsilon_xy(i,j),0;...
     conj(D1(i,j)),-epsilon_x(i,j)+mu,0,-epsilon_xy(i,j);...
     epsilon_xy(i,j),0,epsilon_y(i,j)-mu,D2(i,j);...
     0,-epsilon_xy(i,j),conj(D2(i,j)),-epsilon_y(i,j)+mu];
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
   G11r(:,:,k) = fftshift(ifft2(G11(:,:,k)))/(n_q*n_q);
   G12r(:,:,k) = fftshift(ifft2(G12(:,:,k)))/(n_q*n_q);
   G13r(:,:,k) = fftshift(ifft2(G13(:,:,k)))/(n_q*n_q);
   G14r(:,:,k) = fftshift(ifft2(G14(:,:,k)))/(n_q*n_q);
   G21r(:,:,k) = fftshift(ifft2(G21(:,:,k)))/(n_q*n_q);
   G22r(:,:,k) = fftshift(ifft2(G22(:,:,k)))/(n_q*n_q);
   G23r(:,:,k) = fftshift(ifft2(G23(:,:,k)))/(n_q*n_q);
   G24r(:,:,k) = fftshift(ifft2(G24(:,:,k)))/(n_q*n_q);
   G31r(:,:,k) = fftshift(ifft2(G31(:,:,k)))/(n_q*n_q);
   G32r(:,:,k) = fftshift(ifft2(G32(:,:,k)))/(n_q*n_q);
   G33r(:,:,k) = fftshift(ifft2(G33(:,:,k)))/(n_q*n_q);
   G34r(:,:,k) = fftshift(ifft2(G34(:,:,k)))/(n_q*n_q);
   G41r(:,:,k) = fftshift(ifft2(G41(:,:,k)))/(n_q*n_q);
   G42r(:,:,k) = fftshift(ifft2(G42(:,:,k)))/(n_q*n_q);
   G43r(:,:,k) = fftshift(ifft2(G43(:,:,k)))/(n_q*n_q);
   G44r(:,:,k) = fftshift(ifft2(G44(:,:,k)))/(n_q*n_q);
end
Ak = -imag(G11+G33);
imagesc(Ak(:,:,n_E))
colormap('gray');
axis equal
T11 = zeros(1,n_E);
T12 = zeros(1,n_E);
T13 = zeros(1,n_E);
T14 = zeros(1,n_E);
T21 = zeros(1,n_E);
T22 = zeros(1,n_E);
T23 = zeros(1,n_E);
T24 = zeros(1,n_E);
T31 = zeros(1,n_E);
T32 = zeros(1,n_E);
T33 = zeros(1,n_E);
T34 = zeros(1,n_E);
T41 = zeros(1,n_E);
T42 = zeros(1,n_E);
T43 = zeros(1,n_E);
T44 = zeros(1,n_E);
sumG_ind = (n_q+1)/2;
for k=1:n_E
    A = V_intra;
    B = [G11r(sumG_ind,sumG_ind,k),G12r(sumG_ind,sumG_ind,k),...
        G13r(sumG_ind,sumG_ind,k),G14r(sumG_ind,sumG_ind,k);...
        G21r(sumG_ind,sumG_ind,k),G22r(sumG_ind,sumG_ind,k),...
        G23r(sumG_ind,sumG_ind,k),G24r(sumG_ind,sumG_ind,k);...
        G31r(sumG_ind,sumG_ind,k),G32r(sumG_ind,sumG_ind,k),...
        G33r(sumG_ind,sumG_ind,k),G34r(sumG_ind,sumG_ind,k);...
        G41r(sumG_ind,sumG_ind,k),G42r(sumG_ind,sumG_ind,k),...
        G43r(sumG_ind,sumG_ind,k),G44r(sumG_ind,sumG_ind,k)];
    C = I - A*B;
    T = I/C*A;
    T11(k) = T(1,1);
    T12(k) = T(1,2);
    T13(k) = T(1,3);
    T14(k) = T(1,4);
    T21(k) = T(2,1);
    T22(k) = T(2,2);
    T23(k) = T(2,3);
    T24(k) = T(2,4);
    T31(k) = T(3,1);
    T32(k) = T(3,2);
    T33(k) = T(3,3);
    T34(k) = T(3,4); 
    T41(k) = T(4,1);
    T42(k) = T(4,2);
    T43(k) = T(4,3);
    T44(k) = T(4,4); 
end

A11r = zeros(n_q,n_q,n_E);
A22r = zeros(n_q,n_q,n_E);
A33r = zeros(n_q,n_q,n_E);
A44r = zeros(n_q,n_q,n_E);
for k=1:n_E
    k
for i=1:n_q
for j=1:n_q
G1 = [G11r(i,j,k),G12r(i,j,k),G13r(i,j,k),G14r(i,j,k);...
    G21r(i,j,k),G22r(i,j,k),G23r(i,j,k),G24r(i,j,k);...
    G31r(i,j,k),G32r(i,j,k),G33r(i,j,k),G34r(i,j,k);...
    G41r(i,j,k),G42r(i,j,k),G43r(i,j,k),G44r(i,j,k)];
G2 = [G11r(n_q-i+1,n_q-j+1,k),G12r(n_q-i+1,n_q-j+1,k),...
    G13r(n_q-i+1,n_q-j+1,k),G14r(n_q-i+1,n_q-j+1,k);...
    G21r(n_q-i+1,n_q-j+1,k),G22r(n_q-i+1,n_q-j+1,k),...
    G23r(n_q-i+1,n_q-j+1,k),G24r(n_q-i+1,n_q-j+1,k);...
    G31r(n_q-i+1,n_q-j+1,k),G32r(n_q-i+1,n_q-j+1,k),...
    G33r(n_q-i+1,n_q-j+1,k),G34r(n_q-i+1,n_q-j+1,k);...
    G41r(n_q-i+1,n_q-j+1,k),G42r(n_q-i+1,n_q-j+1,k),...
    G43r(n_q-i+1,n_q-j+1,k),G44r(n_q-i+1,n_q-j+1,k)];
T = [T11(k),T12(k),T13(k),T14(k);...
    T21(k),T22(k),T23(k),T24(k);...
    T31(k),T32(k),T33(k),T34(k);...
    T41(k),T42(k),T43(k),T44(k)];
A = G2.*(T*G1);
A11r(i,j,k) = A(1,1);
A22r(i,j,k) = A(2,2);
A33r(i,j,k) = A(3,3);
A44r(i,j,k) = A(4,4);
end
end
end
% dnr = zeros(n_q,n_q,n_E);
% dnq = zeros(n_q,n_q,n_E);
% for k=1:n_E
%     dnq(:,:,k) = -imag(fftshift(fft2(A11r(1:n_q,1:n_q,k))))/pi;
%     dnr(:,:,k) = ifft2(ifftshift(dnq(:,:,k)));
% end
dnr = zeros(n_q,n_q,n_E);
dnq = zeros(n_q,n_q,n_E);
for k=1:n_E
    dnq(:,:,k) = -imag(fftshift(fft2(ifftshift(A11r(:,:,k)+A33r(:,:,k)))))/pi;
    dnr(:,:,k) = fftshift(fft2(dnq(:,:,k)));
end


imagesc(dnq(:,:,n_E));
colormap('gray')
axis equal;

% figure('name','dnq');
% colorscale = 'gray';
% for k=1:n_E
% subplot(ceil(n_E/4),4,k);
% imagesc((abs(dnq(:,:,k))));
% colormap(colorscale);
% axis equal
% xticks([1 n_q-1]);
% xticklabels({'-\pi','\pi'});
% yticks([1 n_q-1]);
% yticklabels({'-\pi','\pi'});
% xlabel('qx');
% ylabel('qy');
% title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
% end
