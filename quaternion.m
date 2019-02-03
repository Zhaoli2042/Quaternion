% rotation quaternion
% define theta, phi and thi in degrees
theta = 0;
phi = 0; % rotation within xy plane, around the z axis
thi = 0;
theta = deg2rad(theta);
phi = deg2rad(phi); % rotation within xy plane, around the z axis
thi = deg2rad(thi);
q0 = cos(theta/2)*cos(0.5*(phi+thi));
q1 = sin(theta/2)*cos(0.5*(phi-thi));
q2 = sin(theta/2)*sin(0.5*(phi-thi));
q3 = cos(theta/2)*sin(0.5*(phi+thi));

rot_mat=[1-2*(q2*q2+q3*q3) 2*(q1*q2-q0*q3) 2*(q1*q3+q2*q2);...
         2*(q1*q2+q0*q3) 1-2*(q1*q1+q3*q3) 2*(q2*q3-q0*q1);...
         2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) 1-2*(q1*q1+q2*q2)];
rot_mat_T=rot_mat';
norm=q0*q0+q1*q1+q2*q2+q3*q3;

dist=1; % bond length set to be 0.5

sigma=1;

box=10;
n=5;

pos=zeros(n,3);
rota=zeros(n,3); % 3 rotation angles of center of mass

pos(1,:)=[0.5*box,0,0];
% set the molecule to be benzene

for i=2:1:n
    isafe = 0;
    num=0;
    fprintf('Placing Molecule Number %.0f\n',i);
    while (isafe == 0)
        num=num+1;
        the_new=randi(360);phi_new=randi(360);thi_new=randi(360);
        rota(i,:)=[the_new phi_new thi_new];
        pos(i,:)=[0.5*box*(2*rand-1) 0.5*box*(2*rand-1) 0.5*box*(2*rand-1)];
        for j = 1:1:6 % six atoms in this molecule
            posj=pos(i,:)+[dist dist dist]/sqrt(3)*quat([theta,phi+(60)*(j-1),thi]);
                for k = 1:1:i-1
                    for l = 1:1:6
                        posl=pos(k,:)+[dist dist dist]/sqrt(3)*quat([theta,phi+(60)*(k-1),thi]);
                        rxjl=posj(1)-posl(1);
                        ryjl=posj(2)-posl(2);
                        rzjl=posj(3)-posl(3);
                        % minimum image convention
                        rxjl=rxjl-box*round(rxjl/box);
                        ryjl=ryjl-box*round(ryjl/box);
                        rzjl=rzjl-box*round(rzjl/box);
                        rjlsq=rxjl*rxjl+ryjl*ryjl+rzjl*rzjl;
                        if (rjlsq < sigma)
                            isafe = 0;
                        else
                            isafe = 1;
                        end
                    end
                end
        end
    end
end


% dump the positions to a real coordinate matrix
matrix = zeros(1,5);
for i = 1:1:n
    for k=1:1:6
        posk=pos(i,:)+[dist dist dist]/sqrt(3)*quat([theta,phi+(60)*(k-1),thi]);
        aa=size(matrix);aa=aa(1);
        if (matrix(1,1) == 0)
            matrix(i,:)=[aa 1 posk];
        else
            matrix=[matrix;aa 1 posk];
        end
    end
end
aa=size(matrix);aa=aa(1);
fid=fopen('pos_benzene.txt','w');
fprintf(fid,'ITEM: TIMESTEP\n');
fprintf(fid,'%.0f\n',1000);
fprintf(fid,'ITEM: NUMBER OF ATOMS\n');
fprintf(fid,'%.0f\n',aa);
fprintf(fid,'ITEM: BOX BOUNDS pp pp pp\n');
% print box information to txt file
for x=1:1:3
    fprintf(fid,'%.0f %.0f\n',-0.5*box,0.5*box);
end
fprintf(fid,'ITEM: ATOMS id type xu yu zu \n');
% print position information to txt file

for x=1:1:aa
    fprintf(fid,'%.0f %.0f %.6f %.6f %.6f\n',matrix(x,1),matrix(x,2),matrix(x,3),matrix(x,4),matrix(x,5));
end