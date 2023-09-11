clear;
clc;
[file_name,file_path] = uigetfile('*.txt','��ѡ��Ƕ��ļ�','�Ƕ��ļ�����.txt');
ang = load(strcat([file_path,file_name]));
A=ang(:,1:5);%A����Ϊ�Ƕ��ļ����ݣ���վ��/����׼��/����׼��/�Ƕ�/��������
[r,c]=size(A);
[file_name,file_path] = uigetfile('*.txt','��ѡ��������ļ�','�������ļ�����.txt');
point = load(strcat([file_path,file_name]));
P=point(:,1:3);%P����Ϊ�������ļ����ݣ�����/X����/Y����
[r1,c1]=size(P);
point_known = P(r1,1);%��֪�����
[file_name,file_path] = uigetfile('*.txt','��ѡ��߳��۲�ֵ�ļ�','�߳��۲�ֵ�ļ�����.txt');
edge = load(strcat([file_path,file_name]));
E=edge(:,1:4);%E����Ϊ�߳��ļ����ݣ����/�յ�/�߳��۲�ֵ/��������
[r2,c2]=size(E);
B = zeros(r+r2,2*(r1-1));%B����Ϊδ֪����������ϵ����
l = zeros(r+r2,1);%l����Ϊ�������о��� V=B*X-l
length_back = zeros(r2,1);%length_back����Ϊ�߳��Ľ������귴��ֵ
%����Ƕ�ϵ����
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
%����߳�ϵ����
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
%��֪��ϵ����0
for i=1:r+r2
    for j = 1:point_known
        B(i,j*2-1)=0;
        B(i,j*2)=0;
    end
end
Coor=zeros(3,2);
q = zeros(r+r2,r+r2);%����Ȩ��
%����Ƕ�Ȩֵ
for i=1:r
    q(i,i)=1;
end
%����߳�Ȩֵ
for i = r+1:r+r2
    q(i,i)=100/E(i-r,3);
end
%����Ƕȳ�ϵ����
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
%�������귴��߳�ֵ
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
%����߳���ϵ����
for i=r+1:r+r2
    l(i,1)=E(i-r,3)-length_back(i-r);
end
NBB=B'*q*B;
W=B'*q*l;
X=pinv(NBB)*W;
V=B*X-l;
Q = zeros(r2,r2);%Q����Ϊ�߳��۲�ֵ��Ȩ��
for i=r+1:r+r2
    Q(i-r,i-r)=q(i,i);
end
%�����������
for i=1:r2
    D=100*inv(Q);
    sigma_edge(i)=sqrt(D(i,i));
end
sigma_edge=sigma_edge';
result_X=zeros(r1-1,3);%result_X����Ϊ����ƽ��ֵ���󣬵���/X����/Y����
%��������ƽ��ֵ
for i=1:r1-1
    result_X(i,1)=i;
    result_X(i,2)=P(i,2)+X(2*i-1);
    result_X(i,3)=P(i,3)+X(2*i);
end
result_A=zeros(r,2);%result_A����Ϊ�Ƕ�ƽ��ֵ���󣬽���/�Ƕ�ƽ��ֵ���ȣ�
%����Ƕ�ƽ��ֵ
for i=1:r
    result_A(i,1)=i;
    result_A(i,2)=V(i,1)/3600+A(i,4);
end
result_final_A=zeros(r,4);%result_final_A����Ϊ�Ƕ�ƽ��ֵ���󣬽���/��/��/��
%�Ƕ�ת����λ��Ϊ�ȷ���
for i=1:r
    degrees=result_A(i,2);
    result_final_A(i,1)=i;
    result_final_A(i,2) = fix(degrees);
    result_final_A(i,3) = fix((degrees - result_final_A(i,2)) * 60);
    result_final_A(i,4) = (degrees - result_final_A(i,2) - result_final_A(i,3)/60)*3600;
    result_final_A(i,4)=round(result_final_A(i,4));
end
jdgcz=zeros(r,4);%jdgcz����Ϊ�Ƕȹ۲�ֵ���󣬽���/��/��/��
for i=1:r
    degrees=A(i,4);
    jdgcz(i,1)=i;
    jdgcz(i,2) = fix(degrees);
    jdgcz(i,3) = fix((degrees - jdgcz(i,2)) * 60);
    jdgcz(i,4) = (degrees - jdgcz(i,2) - jdgcz(i,3)/60)*3600;
    jdgcz(i,4)=round(jdgcz(i,4));
end
result_E = zeros(r2,3);%result_E����Ϊ�߳�ƽ��ֵ�������/�յ�/�߳�ƽ��ֵ
%����߳�ƽ��ֵ
for i=r+1:r2+r
    result_E(i-r,1)=E(i-r,1);
    result_E(i-r,2)=E(i-r,2);
    result_E(i-r,3)=E(i-r,3)+V(i,1);
    result_E(i-r,3) = roundn(result_E(i-r,3),-4);
end
%�������
[file_name,file_path]=uigetfile('*.txt','��ѡ�������λ��');
fn=fopen(strcat([file_path,file_name]),'w');
fprintf(fn,'%s\n\n','��Ŷ��ձ�');
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
fprintf(fn,'%s\n\n','�Ƕȹ۲�ֵΪ��');
fprintf(fn,'%s\n\n','��վ��	����׼��    ����׼��    �Ƕȹ۲�ֵ');
for i=1:r
    fprintf(fn,'%d',A(i,1));
    fprintf(fn,'%15d',A(i,2));
    fprintf(fn,'%15d',A(i,3));
    fprintf(fn,'%15.4f\n',A(i,4));
end
fprintf(fn,'%s\n\n','�߳��۲�ֵΪ��');
fprintf(fn,'%s\n\n','��� �յ�	�߳��۲�ֵ');
for i=1:r2
    fprintf(fn,'%d',E(i,1));
    fprintf(fn,'%10d',E(i,2));
    fprintf(fn,'%20.4f\n',E(i,3));
end
fprintf(fn,'%s\n\n','���������������ƽ��ֵ��');
fprintf(fn,'%s\n\n','����     X����     Y����     X����������)     Y���������ף�     Xƽ��ֵ        Yƽ��ֵ');
for i=1:r1-1
    fprintf(fn,'%d',P(i,1));
    fprintf(fn,'%15.4f',P(i,2));
    fprintf(fn,'%15.4f',P(i,3));
    fprintf(fn,'%10.4f',X(2*i-1,1));
    fprintf(fn,'%10.4f',X(2*i,1));
    fprintf(fn,'%15.4f',result_X(i,2));
    fprintf(fn,'%15.4f\n',result_X(i,3));
end
fprintf(fn,'%s\n\n','�Ƕȸ�������ƽ��ֵ��');
fprintf(fn,'%s\n\n','�����    ��վ��   ����׼��   ����׼��   �Ƕȹ۲�ֵ(����)      �Ƕȸ��������壩     �Ƕ�ƽ��ֵ(����)');
for i=1:r
    fprintf(fn,'%s','��');
    fprintf(fn,'%d',i);
    fprintf(fn,'%14d',A(i,1));
    fprintf(fn,'%14d',A(i,2));
    fprintf(fn,'%14d',A(i,3));
    fprintf(fn,'%12d',jdgcz(i,2:4));
    fprintf(fn,'%15.4f',V(i,1));
    fprintf(fn,'%10d',result_final_A(i,2:4));
    fprintf(fn,'\n');
end
fprintf(fn,'%s\n\n','�߳���������ƽ��ֵ��');
fprintf(fn,'%s\n\n','���   �յ�   �߳��۲�ֵ    �߳�������     �߳�ƽ��ֵ     �߳���������)');
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
title('�������������Բ');
grid on;