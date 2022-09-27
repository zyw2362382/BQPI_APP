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
k
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
%%
top_margin = 0.03; % top margin
btm_margin = 0.03; % bottom margin
left_margin = 0.03;% left margin
right_margin = 0.15;% right margin
 
fig_margin = 0.08; % margin beween figures(sub) 

row = ceil(n_E/4); % rows
col = 4; % cols
 
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
       C = max(max(abs(dnqx(:,:,k))));
       imagesc(imadjust(abs(dnqx(:,:,k))/C));
       axis equal
       colormap('gray');
       % title, labels
       xticks([1 n_q]);
       xticklabels({'-\pi','\pi'});
       yticks([1 n_q]);
       yticklabels({'-\pi','\pi'});
       xlabel('qx');
       ylabel('qy');
       title([num2str(-E+(k-1)*2*E/(n_E-1)) 'eV']);
    end   
end
% draw colorbar
axes('position', [1-right_margin-fig_margin, btm_margin, 0.2, 1-(top_margin+btm_margin)]);
axis off;
colorbar();
% caxis(clim);
colorbar('Ticks',[0 1],...
         'TickLabels',{'Low density ','High density '})
%%

% 0.1305 + -0.5951*(cos(pi*qx)+cos(pi*qy))/2+0.1636*cos(pi*qx).*cos(pi*qy)-0.0519*(cos(2*pi*qx)+cos(2*pi*qy))/2 -0.1117*(cos(2*pi*qx).*cos(pi*qy)+cos(pi*qx).*cos(2*pi*qy))/2 + 0.0510*cos(2*pi*qx).*cos(2*pi*qy)
%0.025/2*(cos(pi*qx)-cos(pi*qy))
% 0.1305 -0.5951*(cos(kx)+cos(ky))/2+0.1636*cos(kx).*cos(ky)-0.0519*(cos(2*kx)+cos(2*ky))/2 -0.1117*(cos(2*kx).*cos(ky)+cos(kx).*cos(2*ky))/2 + 0.0510*cos(2*kx).*cos(2*ky)
% epsilon_alpha = 1/2*(((- 1 - 2*cos(kx) - 2*0.1*cos(ky)) + (- 1 - 2*cos(ky) - 2*0.1*cos(kx))) - sqrt(((- 1 - 2*cos(kx) - 2*0.1*cos(ky)) - (- 1 - 2*cos(ky) - 2*0.1*cos(kx))).^2 + 4*(- 2*0.1*sin(kx).*sin(ky)).^2))
% epsilon_beta = 1/2*(((- 1 - 2*cos(kx) - 2*0.1*cos(ky)) + (- 1 - 2*cos(ky) - 2*0.1*cos(kx))) + sqrt(((- 1 - 2*cos(kx) - 2*0.1*cos(ky)) - (- 1 - 2*cos(ky) - 2*0.1*cos(kx))).^2 + 4*(- 2*0.1*sin(kx).*sin(ky)).^2))
%D = 0.00035*(cos(kx) - cos(ky))