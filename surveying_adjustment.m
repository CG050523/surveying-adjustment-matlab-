clear;
clc;
[file_name,file_path] = uigetfile('*.txt','请选择角度文件','角度文件内容.txt');
ang = load(strcat([file_path,file_name]));
A=ang(:,1:5);%A矩阵为角度文件内容，测站点/左照准点/右照准点/角度/测角中误差
[r,c]=size(A);
[file_name,file_path] = uigetfile('*.txt','请选择点坐标文件','点坐标文件内容.txt');
point = load(strcat([file_path,file_name]));
P=point(:,1:3);%P矩阵为点坐标文件内容，点名/X坐标/Y坐标
[r1,c1]=size(P);
point_known = P(r1,1);%已知点个数
[file_name,file_path] = uigetfile('*.txt','请选择边长观测值文件','边长观测值文件内容.txt');
edge = load(strcat([file_path,file_name]));
E=edge(:,1:4);%E矩阵为边长文件内容，起点/终点/边长观测值/测边中误差
[r2,c2]=size(E);
B = zeros(r+r2,2*(r1-1));%B矩阵为未知数改正数的系数阵
l = zeros(r+r2,1);%l矩阵为常数项列矩阵 V=B*X-l
length_back = zeros(r2,1);%length_back矩阵为边长的近似坐标反算值
%计算角度系数阵
for i=1:r
    for m=1:r1-1
        if A(i,1)==P(m,1)
            Coor(1,1)=P(m,2);
            Coor(1,2)=P(m,3);
        end
    end
    for n=1:r1-1
        if A(i,2)==P(n,1)
            Coor(2,1)=P(n,2);
            Coor(2,2)=P(n,3);
        end
    end
    for o=1:r1-1
        if A(i,3)==P(o,1)
            Coor(3,1)=P(o,2);
            Coor(3,2)=P(o,3);
        end
    end
    ajk =206265*(Coor(3,2)-Coor(1,2))/((Coor(3,1)-Coor(1,1))*(Coor(3,1)-Coor(1,1))+(Coor(3,2)-Coor(1,2))*(Coor(3,2)-Coor(1,2)));
    ajh =206265*(Coor(2,2)-Coor(1,2))/((Coor(2,1)-Coor(1,1))*(Coor(2,1)-Coor(1,1))+(Coor(2,2)-Coor(1,2))*(Coor(2,2)-Coor(1,2)));
    bjk =206265*(Coor(3,1)-Coor(1,1))/((Coor(3,1)-Coor(1,1))*(Coor(3,1)-Coor(1,1))+(Coor(3,2)-Coor(1,2))*(Coor(3,2)-Coor(1,2)));
    bjh =206265*(Coor(2,1)-Coor(1,1))/((Coor(2,1)-Coor(1,1))*(Coor(2,1)-Coor(1,1))+(Coor(2,2)-Coor(1,2))*(Coor(2,2)-Coor(1,2)));
    
    B(i,A(i,1)*2-1)=ajk-ajh;B(i,A(i,1)*2)=-bjk+bjh;
    B(i,A(i,3)*2-1)=-ajk;B(i,A(i,3)*2)=bjk;
    B(i,A(i,2)*2-1)=ajh;B(i,A(i,2)*2)=-bjh;
    
end
%计算边长系数阵
for i=r+1:r+r2
    for m=1:r1-1
        if E(i-r,1)==P(m,1)
            Coor(1,1)=P(m,2);
            Coor(1,2)=P(m,3);
        end
    end
    for n=1:r1-1
        if E(i-r,2)==P(n,1)
            Coor(2,1)=P(n,2);
            Coor(2,2)=P(n,3);
        end
    end
    cjk = (Coor(2,1)-Coor(1,1))/sqrt((Coor(2,1)-Coor(1,1))*(Coor(2,1)-Coor(1,1))+(Coor(2,2)-Coor(1,2))*(Coor(2,2)-Coor(1,2)));
    djk = (Coor(2,2)-Coor(1,2))/sqrt((Coor(2,1)-Coor(1,1))*(Coor(2,1)-Coor(1,1))+(Coor(2,2)-Coor(1,2))*(Coor(2,2)-Coor(1,2)));
    B(i,E(i-r,1)*2-1)=-cjk;B(i,E(i-r,1)*2)=-djk;
    B(i,E(i-r,2)*2-1)=cjk;B(i,E(i-r,2)*2)=djk;
end
%已知点系数赋0
for i=1:r+r2
    for j = 1:point_known
        B(i,j*2-1)=0;
        B(i,j*2)=0;
    end
end
Coor=zeros(3,2);
q = zeros(r+r2,r+r2);%定义权阵
%计算角度权值
for i=1:r
    q(i,i)=1;
end
%计算边长权值
for i = r+1:r+r2
    q(i,i)=100/E(i-r,3);
end
%计算角度常系数阵
for i=1:r
    for j=1:r1-1
        if A(i,1)==P(j,1)
            Coor(1,1)=P(j,2);Coor(1,2)=P(j,3);
        end
        if A(i,2)==P(j,1)
            Coor(2,1)=P(j,2);Coor(2,2)=P(j,3);
        end
        if A(i,3)==P(j,1)
            Coor(3,1)=P(j,2);Coor(3,2)=P(j,3);
        end
    end
    detx=Coor(3,1)-Coor(1,1);dety=Coor(3,2)-Coor(1,2);
    arc_R = rad2deg(atan(abs(dety)/abs(detx)));
    if (detx>0)&&(dety>0)
        arc_R=arc_R;
    end
    if (detx<0)&&(dety>0)
        arc_R=180-arc_R;
    end
    if (detx<0)&&(dety<0)
        arc_R=180+arc_R;
    end
    if (detx>0)&&(dety<0)
        arc_R=2*180-arc_R;
    end
    if (detx==0)&&(dety<0)
        arc_R=180+arc_R;
    end
    detx=Coor(2,1)-Coor(1,1);dety=Coor(2,2)-Coor(1,2);
    arc_L = rad2deg(atan(abs(dety)/abs(detx)));
    if (detx>0)&&(dety>0)
        arc_L=arc_L;
    end
    if (detx<0)&&(dety>0)
        arc_L=180-arc_L;
    end
    if (detx<0)&&(dety<0)
        arc_L=180+arc_L;
    end
    if (detx>0)&&(dety<0)
        arc_L=2*180-arc_L;
    end
    if arc_R>arc_L
        arc(i)=arc_R-arc_L;
    else
        arc(i)=arc_R-arc_L+360;
    end
    l(i) =(A(i,4)-arc(i))*3600;
    
end
%近似坐标反算边长值
for i=r+1:r+r2
    for m=1:r1-1
        if E(i-r,1)==P(m,1)
            Coor(1,1)=P(m,2);
            Coor(1,2)=P(m,3);
        end
    end
    for n=1:r1-1
        if E(i-r,2)==P(n,1)
            Coor(2,1)=P(n,2);
            Coor(2,2)=P(n,3);
        end
    end
    detx=Coor(2,1)-Coor(1,1);dety=Coor(2,2)-Coor(1,2);
    length_back(i-r)=sqrt(detx*detx+dety*dety);
end
%计算边长常系数阵
for i=r+1:r+r2
    l(i,1)=E(i-r,3)-length_back(i-r);
end
NBB=B'*q*B;
W=B'*q*l;
X=pinv(NBB)*W;
V=B*X-l;
Q = zeros(r2,r2);%Q矩阵为边长观测值的权阵
for i=r+1:r+r2
    Q(i-r,i-r)=q(i,i);
end
%计算测边中误差
for i=1:r2
    D=100*inv(Q);
    sigma_edge(i)=sqrt(D(i,i));
end
sigma_edge=sigma_edge';
result_X=zeros(r1-1,3);%result_X矩阵为坐标平差值矩阵，点名/X坐标/Y坐标
%计算坐标平差值
for i=1:r1-1
    result_X(i,1)=i;
    result_X(i,2)=P(i,2)+X(2*i-1);
    result_X(i,3)=P(i,3)+X(2*i);
end
result_A=zeros(r,2);%result_A矩阵为角度平差值矩阵，角名/角度平差值（度）
%计算角度平差值
for i=1:r
    result_A(i,1)=i;
    result_A(i,2)=V(i,1)/3600+A(i,4);
end
result_final_A=zeros(r,4);%result_final_A矩阵为角度平差值矩阵，角名/度/分/秒
%角度转换单位度为度分秒
for i=1:r
    degrees=result_A(i,2);
    result_final_A(i,1)=i;
    result_final_A(i,2) = fix(degrees);
    result_final_A(i,3) = fix((degrees - result_final_A(i,2)) * 60);
    result_final_A(i,4) = (degrees - result_final_A(i,2) - result_final_A(i,3)/60)*3600;
    result_final_A(i,4)=round(result_final_A(i,4));
end
jdgcz=zeros(r,4);%jdgcz矩阵为角度观测值矩阵，角名/度/分/秒
for i=1:r
    degrees=A(i,4);
    jdgcz(i,1)=i;
    jdgcz(i,2) = fix(degrees);
    jdgcz(i,3) = fix((degrees - jdgcz(i,2)) * 60);
    jdgcz(i,4) = (degrees - jdgcz(i,2) - jdgcz(i,3)/60)*3600;
    jdgcz(i,4)=round(jdgcz(i,4));
end
result_E = zeros(r2,3);%result_E矩阵为边长平差值矩阵，起点/终点/边长平差值
%计算边长平差值
for i=r+1:r2+r
    result_E(i-r,1)=E(i-r,1);
    result_E(i-r,2)=E(i-r,2);
    result_E(i-r,3)=E(i-r,3)+V(i,1);
    result_E(i-r,3) = roundn(result_E(i-r,3),-4);
end
%数据输出
[file_name,file_path]=uigetfile('*.txt','请选择结果输出位置');
fn=fopen(strcat([file_path,file_name]),'w');
fprintf(fn,'%s\n\n','点号对照表');
for i = 1:r1-1
    if i <= point_known
        fprintf(fn,'%c',char(64+i));
        fprintf(fn,'%14s','  -----------  ');
        fprintf(fn,'%d\n',i);
    else
        fprintf(fn,'%s','P');
        fprintf(fn,'%d',i-point_known);
        fprintf(fn,'%14s','  -----------  ');
        fprintf(fn,'%d\n',i);
    end
end
fprintf(fn,'\n');
fprintf(fn,'%s\n\n','角度观测值为：');
fprintf(fn,'%s\n\n','测站点	左照准点    右照准点    角度观测值');
for i=1:r
    fprintf(fn,'%d',A(i,1));
    fprintf(fn,'%15d',A(i,2));
    fprintf(fn,'%15d',A(i,3));
    fprintf(fn,'%15.4f\n',A(i,4));
end
fprintf(fn,'%s\n\n','边长观测值为：');
fprintf(fn,'%s\n\n','起点 终点	边长观测值');
for i=1:r2
    fprintf(fn,'%d',E(i,1));
    fprintf(fn,'%10d',E(i,2));
    fprintf(fn,'%20.4f\n',E(i,3));
end
fprintf(fn,'%s\n\n','坐标改正数和坐标平差值：');
fprintf(fn,'%s\n\n','点名     X坐标     Y坐标     X改正数（米)     Y改正数（米）     X平差值        Y平差值');
for i=1:r1-1
    fprintf(fn,'%d',P(i,1));
    fprintf(fn,'%15.4f',P(i,2));
    fprintf(fn,'%15.4f',P(i,3));
    fprintf(fn,'%10.4f',X(2*i-1,1));
    fprintf(fn,'%10.4f',X(2*i,1));
    fprintf(fn,'%15.4f',result_X(i,2));
    fprintf(fn,'%15.4f\n',result_X(i,3));
end
fprintf(fn,'%s\n\n','角度改正数和平差值：');
fprintf(fn,'%s\n\n','角序号    测站点   左照准点   右照准点   角度观测值(°′″)      角度改正数（″）     角度平差值(°′″)');
for i=1:r
    fprintf(fn,'%s','β');
    fprintf(fn,'%d',i);
    fprintf(fn,'%14d',A(i,1));
    fprintf(fn,'%14d',A(i,2));
    fprintf(fn,'%14d',A(i,3));
    fprintf(fn,'%12d',jdgcz(i,2:4));
    fprintf(fn,'%15.4f',V(i,1));
    fprintf(fn,'%10d',result_final_A(i,2:4));
    fprintf(fn,'\n');
end
fprintf(fn,'%s\n\n','边长改正数和平差值：');
fprintf(fn,'%s\n\n','起点   终点   边长观测值    边长改正数     边长平差值     边长中误差（毫米)');
for i=1:r2
    fprintf(fn,'%d',E(i,1));
    fprintf(fn,'%10d',E(i,2));
    fprintf(fn,'%20.4f',E(i,3));
    fprintf(fn,'%14.4f',V(r+i,1));
    fprintf(fn,'%14.4f',result_E(i,3));
    fprintf(fn,'%14.4f\n',sigma_edge(i,1));
end
Coor = zeros(2,2);
figure(1);
for i=1:r2
    for m=1:r1-1
        if edge(i,1)==point(m,1)
            Coor(1,1)=point(m,2);
            Coor(1,2)=point(m,3);
        end
    end
    for n=1:r1-1
        if edge(i,2)==point(n,1)
            Coor(2,1)=point(n,2);
            Coor(2,2)=point(n,3);
        end
    end
    line([Coor(1,2),Coor(2,2)],[Coor(1,1),Coor(2,1)]);
    hold on;
end
for i=1:r1-1
    x(i,1) = P(i,2);
    y(i,1) = P(i,3);
end
scatter(y,x,'fill');
QXX = pinv(NBB);
for i=point_known + 1:r1-1
    K_error = sqrt((QXX(i*2-1,i*2-1)-QXX(i*2,i*2))^2+4*(QXX(i*2-1,i*2)^2));
    E_error(i-point_known,1) = 10*sqrt(5*(QXX(i*2-1,i*2-1)+QXX(i*2,i*2)+K_error));
    F_error(i-point_known,1) = 10*sqrt(5*(QXX(i*2-1,i*2-1)+QXX(i*2,i*2)-K_error));
    direction_error(i-point_known,1) = rad2deg(atan(QXX(2*i-1,2*i)/((QXX(2*i-1,2*i-1)+QXX(2*i,2*i)+K_error)/2-QXX(2*i,2*i))));
end
sita = 0:0.1:2*pi;
for i=1 + point_known:r1-1
    plot(P(i,3)+cos(direction_error(i-point_known,1)*pi/180)*E_error(i-point_known,1)*cos(sita)-sin(direction_error(i-point_known,1)*pi/180)*F_error(i-point_known,1)*sin(sita),P(i,2)+cos(direction_error(i-point_known,1)*pi/180)*F_error(i-point_known,1)*sin(sita)+sin(direction_error(i-point_known,1)*pi/180)*E_error(i-point_known,1)*cos(sita),'m');
end
hold on;
title('导线网及误差椭圆');
grid on;