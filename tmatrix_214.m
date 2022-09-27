n_q = 501; % number of points in k-space
E = 0.001; % energy range in eV(from -1meV to 1meV)
n_E = 11; % number of energy layers
[qx,qy] = meshgrid(linspace(-1,1,n_q));
Epoints = linspace(-E,E,n_E);
I = eye(4);

%the parameters
mu_0 = 1;
t_h = 0.1;
V_m = 0.1;
mu_z = 0.7;
t_z = 0.55;
t_z1 = 0.2;
t = 1;

%the tight binding model to fit superconducting STM spectra
epsilon_xz = - mu_0 - 2*t*cos(pi*qx) - 2*t_h*cos(pi*qy);
epsilon_yz = - mu_0 - 2*t*cos(pi*qy) - 2*t_h*cos(pi*qx);
V_hyb = - 2*V_m*sin(pi*qx).*sin(pi*qy);
epsilon_alpha = 1/2*((epsilon_xz + epsilon_yz) ...
    - sqrt((epsilon_xz - epsilon_yz).^2 + 4*V_hyb.^2));
epsilon_beta = 1/2*((epsilon_xz + epsilon_yz) ...
    + sqrt((epsilon_xz - epsilon_yz).^2 + 4*V_hyb.^2));

%the energy gap on α : β Fermi surface
D0 = 0.00035;
D = D0*(cos(pi*qx) - cos(pi*qy));

%the scattering potential
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
%    H = [epsilon_alpha(i,j),D(i,j),V_hyb(i,j),0;...
%      D(i,j),-epsilon_alpha(i,j),0,-V_hyb(i,j);...
%      V_hyb(i,j),0,epsilon_beta(i,j),D(i,j);...
%      0,-V_hyb(i,j),D(i,j),-epsilon_beta(i,j)];
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
    A = [V_intra,0,V_inter,0;...
        0,-V_intra,0,-V_inter;...
        V_inter,0,V_intra,0;...
        0,-V_inter,0,-V_intra];
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
dnr = zeros(n_q,n_q,n_E);
dnq = zeros(n_q,n_q,n_E);
for k=1:n_E
%     dnr(:,:,k) = -imag(A11r(:,:,k)+A22r(:,:,n_E-k+1)...
%     +A33r(:,:,k)+A44r(:,:,n_E-k+1))/pi/2 ;
%     dnq(:,:,k) = fftshift(fft2(dnr(:,:,k)));
    dnq(:,:,k) = -imag(fftshift(fft2(ifftshift(A11r(:,:,k)+A33r(:,:,k)))))/pi;
    dnr(:,:,k) = fftshift(fft2(dnq(:,:,k)));
end

dnq1 = zeros(2*n_q+1,2*n_q+1,n_E);
for k=1:n_E
dnq_tem=zeros(3*n_q,3*n_q);
dnq_tem(1:n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(1:n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(1:n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq1(:,:,k)=dnq_tem((n_q+1)/2:(n_q+1)/2+2*n_q,(n_q+1)/2:(n_q+1)/2+2*n_q);
end

%%


clim = [0 54];
C = max(max(abs(dnq1(:,:,11))));
imagesc(imadjust(abs(dnq1(:,:,11))/C));
% imagesc(epsilon_alpha);
axis equal
colormap('gray');
xticks([1 n_q 2*n_q]);
xticklabels({'-2\pi','0','2\pi'});
yticks([1 n_q 2*n_q]);
yticklabels({'-2\pi','0','2\pi'});
xlabel('kx');
ylabel('ky');


colorbar();
% caxis(clim);
colorbar('Ticks',[0 1],...
         'TickLabels',{'Low density ','High density '})
%%
top_margin = 0.03; % top margin
btm_margin = 0.03; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.15;% right margin
 
fig_margin = 0; % margin beween figures(sub) 

row = 2; % rows
col = 5; % cols
 
clim = [0 54];
% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;
 
for i = 1 : row
    for j = 1 : col
        % figure position: you can refer to 'help axes' to review the
        % parameter meaning, note that original point is lower-left
        position = [left_margin + (j-1)*(fig_margin+fig_w), ...
           1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
           fig_w, fig_h]
       axes('position', position)
       % draw colorful pictures... 
       k = (i - 1) * col +j;
       imagesc(Ak(:,:,k));
       axis equal
       colormap('gray');
       % title, labels
       xticks([1 n_q]);
       xticklabels({'-\pi','\pi'});
       yticks([1 n_q]);
       yticklabels({'-\pi','\pi'});
       xlabel('kx');
       ylabel('ky');
       title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
    end   
end
% draw colorbar
axes('position', [1-right_margin-fig_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar();caxis(clim);
colorbar('Ticks',[0 50],...
         'TickLabels',{'Low density ','High density '})

%%
top_margin = 0.03; % top margin
btm_margin = 0.03; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.15;% right margin
 
fig_margin = 0.0; % margin beween figures(sub) 

row = 2; % rows
col = 5; % cols
 
clim = [0 54];
% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;
 
for i = 1 : row
    for j = 1 : col
        % figure position: you can refer to 'help axes' to review the
        % parameter meaning, note that original point is lower-left
        position = [left_margin + (j-1)*(fig_margin+fig_w), ...
           1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
           fig_w, fig_h]
       axes('position', position)
       % draw colorful pictures... 
       k = (i - 1) * col +j;
       imagesc(abs(dnr(240:260,240:260,k)));
       axis equal
       colormap('gray');
       % title, labels
       set(gca,'xtick',[],'xticklabel',[])
       set(gca,'ytick',[],'yticklabel',[])
       xlabel('x');
       ylabel('y');
       title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
       
    end
end
% draw colorbar
axes('position', [1-right_margin-fig_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar();caxis(clim);
colorbar('Ticks',[0 50],...
         'TickLabels',{'Low density ','High density '})

%%
top_margin = 0.03; % top margin
btm_margin = 0.03; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.15;% right margin
 
fig_margin = 0.0; % margin beween figures(sub) 

row = 2; % rows
col = 5; % cols
 
clim = [0 54];
% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;
 
for i = 1 : row
    for j = 1 : col
        % figure position: you can refer to 'help axes' to review the
        % parameter meaning, note that original point is lower-left
        position = [left_margin + (j-1)*(fig_margin+fig_w), ...
           1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
           fig_w, fig_h]
       axes('position', position)
       % draw colorful pictures... 
       k = (i - 1) * col +j;
       C = max(max(abs(dnq1(:,:,11))));
       imagesc(imadjust(abs(dnq1(:,:,11))/C));
%        imagesc(abs(dnq1(:,:,k)));
       axis equal
       colormap('gray');
       % title, labels
       xticks([1 n_q 2*n_q]);
       xticklabels({'-2\pi','0','2\pi'});
       yticks([1 n_q 2*n_q]);
       yticklabels({'-2\pi','0','2\pi'});
       xlabel('qx');
       ylabel('qy');
       title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
       
    end
end
% draw colorbar
axes('position', [1-right_margin-fig_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar();
caxis(clim);
colorbar('Ticks',[0 50],...
         'TickLabels',{'Low density ','High density '})
%%

figure('name','dnq');
colorscale = 'gray';
for k=1:n_E
dnq_tem=zeros(3*n_q,3*n_q);
dnq_tem(1:n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(1:n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(1:n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(n_q+1:2*n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,1:n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,n_q+1:2*n_q)=abs(dnq(:,:,k));
dnq_tem(2*n_q+1:3*n_q,2*n_q+1:3*n_q)=abs(dnq(:,:,k));
dnq1=dnq_tem((n_q+1)/2:(n_q+1)/2+2*n_q,(n_q+1)/2:(n_q+1)/2+2*n_q);
subplot(ceil(n_E/4),4,k);
imagesc(dnq1);
colormap(colorscale);
axis equal
xticks([1 n_q 2*n_q]);
xticklabels({'-2\pi','0','2\pi'});
yticks([1 n_q 2*n_q]);
yticklabels({'-2\pi','0','2\pi'});
xlabel('qx');
ylabel('qy');
end
