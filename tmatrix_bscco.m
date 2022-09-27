%the cauculation for bscco:

n_q = 401;% number of points in k-space
E = 0.024;% energy range in eV(from -1meV to 1meV)
n_E = 16; % number of energy layers
[qx,qy] = meshgrid(linspace(-1,1,n_q));
Epoints = linspace(-E,E,n_E);
I = eye(2);
Vs = 0.1;
Vm = 0.0;

%the parameters
t1 = 0.1305; 
t2 = -0.5951;
t3 = 0.1636;
t4 = -0.0519;
t5 = -0.1117;
t6 = 0.0510;

%the tight binding model to fit superconducting STM spectra
E_tb = t1 + t2*(cos(pi*qx)+cos(pi*qy))/2 ...
    +t3*cos(pi*qx).*cos(pi*qy)+ t4*(cos(2*pi*qx)+cos(2*pi*qy))/2 ...
    +t5*(cos(2*pi*qx).*cos(pi*qy)+cos(pi*qx).*cos(2*pi*qy))/2 ...
    + t6*cos(2*pi*qx).*cos(2*pi*qy);

%the energy gap
D0 = 0.025;
D = D0/2*(cos(pi*qx)-cos(pi*qy));

%the scattering potential
V = [Vs + Vm,0;0,-Vs + Vm];

%bias modulation
d = 0.0015;

%calculate the Green's Function
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
   H = [E_tb(i,j),D(i,j);conj(D(i,j)),-E_tb(i,j)];% the Hamiltonian
   G0 = I/((Epoints(k) + 1i*d)*I - H);%the unperturbed Green’s function
   G11(i,j,k) = G0(1,1);
   G12(i,j,k) = G0(1,2);
   G21(i,j,k) = G0(2,1);
   G22(i,j,k) = G0(2,2);
end
end
   G11r(:,:,k) = fftshift(ifft2(G11(:,:,k)))/(n_q*n_q);
   G12r(:,:,k) = fftshift(ifft2(G12(:,:,k)))/(n_q*n_q);
   G21r(:,:,k) = fftshift(ifft2(G21(:,:,k)))/(n_q*n_q);
   G22r(:,:,k) = fftshift(ifft2(G22(:,:,k)))/(n_q*n_q);
end

%plot the LDOS for bscco
Ak = -imag(G11);
figure('name','LDOS')
imagesc(Ak(:,:,12));%choose one energy layer(layer=12 and E=0.024eV) to calculate 

%calculate the T matrix:
T11 = zeros(1,n_E);
T12 = zeros(1,n_E);
T21 = zeros(1,n_E);
T22 = zeros(1,n_E);
sumG_ind = (n_q+1)/2;
for k=1:n_E
    A = [Vs + Vm,0;0,-Vs + Vm];
    B = [G11r(sumG_ind,sumG_ind,k),G12r(sumG_ind,sumG_ind,k);...
        G21r(sumG_ind,sumG_ind,k),G22r(sumG_ind,sumG_ind,k)];
    C = I - A*B;
    T = I/C*A;
    T11(k) = T(1,1);
    T12(k) = T(1,2);
    T21(k) = T(2,1);
    T22(k) = T(2,2);
end

%calculate the change in density of states
A11r = zeros(n_q,n_q,n_E);
A22r = zeros(n_q,n_q,n_E);
for k=1:n_E
for i=1:n_q
for j=1:n_q
G1r = [G11r(i,j,k),G12r(i,j,k);...
    G21r(i,j,k),G22r(i,j,k)];
G2r = [G11r(n_q-i+1,n_q-j+1,k),G12r(n_q-i+1,n_q-j+1,k);...
    G21r(n_q-i+1,n_q-j+1,k),G22r(n_q-i+1,n_q-j+1,k)];
T = [T11(k),T12(k);...
    T21(k),T22(k)];
A = G2r*T*G1r;
A11r(i,j,k) = A(1,1);
A22r(i,j,k) = A(2,2);
end
end
end

dnr = zeros(n_q-1,n_q-1,n_E);
dnq = zeros(n_q-1,n_q-1,n_E);
for k=1:n_E
    dnq(:,:,k) = -imag(fftshift(fft2(A11r(1:n_q-1,1:n_q-1,k))))/pi;
    dnr(:,:,k) = ifft2(ifftshift(dnq(:,:,k)));
end

figure('name','dnq');
colorscale = 'gray';
for k=1:n_E
subplot(ceil(n_E/4),4,k);
imagesc((abs(dnq(:,:,k))));
colormap(colorscale);
axis equal
xticks([1 n_q-1]);
xticklabels({'-\pi','\pi'});
yticks([1 n_q-1]);
yticklabels({'-\pi','\pi'});
xlabel('qx');
ylabel('qy');
title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
end