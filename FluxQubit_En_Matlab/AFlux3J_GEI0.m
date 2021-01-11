function [En0, Ef0, Ia0, Ia0ss, Ia0bs]=AFlux3J_GEI0(fArray, mac)


N = 8;  % ���ڼ���Ļ�ʸ������, ȡֵԽ��Խ��ȷ
mu=mac(1);    % mu = 40;     %Ej=mu*Ec
alpha=mac(2); % alpha =0.8;   %small size junction

% fArray�ĳ���, ȷ��EI�����С
fNmax = length(fArray);

% �������� Hm�����չ�ʽ�����и�Hm��ֵ
Hm = zeros((2*N+1)^2); %289*289

%��Hm����f�ı�Ĳ��ָ�ֵ
for k = -N:N
    for l = -N:N
        row = (k+N)*(2*N+1)+l+N+1;
        
        %�Խ�Ԫ (k,l)��,(k,l)��
        col = row;
        Hm(row,col) = Hm(row,col)+4/mu/(1+2*alpha)*(k^2 +l^2 +alpha*(k+l)^2) +2+alpha;
        %Hm(row,col) = Hm(row,col)+2/mu *(k+l)^2 + 2/mu/(1+2*alpha) * (k-l)^2 +2+alpha;
        
        %�ǶԽ�Ԫ (k,l)��,(k-1,l)��
        k1 = k-1;    
        l1 = l;
        if k1>=-N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %�ǶԽ�Ԫ (k,l)��,(k+1,l)��
        k1 = k+1; 
        l1 = l;
        if k1<=N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %�ǶԽ�Ԫ (k,l)��,(k,l-1)��
        k1 = k;     
        l1 = l-1;
        if l1>=-N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %�ǶԽ�Ԫ (k,l)��,(k,l+1)��
        k1 = k;     
        l1 = l+1;
        if l1<=N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
      
    end
end

%��ÿ��fֵ����Hm��f�仯�Ĳ��ָ�ֵ
%�õ�������Hm�����Hm�ı���ֵ�뱾��ʸ�������������Ԫt��I
englev=zeros(6,fNmax); %�����ܼ�
Isp=zeros(6,fNmax); %����I��������
Idp=zeros(6,fNmax); %����Iλ�Ƶ�������

Ibs=zeros(6,fNmax); %��ڳ�������

for fN = 1:fNmax
    f = fArray(fN); 
    H = Hm;
    
    %��Hm��f�仯�Ĳ��ָ�ֵ
    for k = -N:N
        for l = -N:N
            
            row = (k+N)*(2*N+1)+l+N+1;
            
            %�ǶԽ�Ԫ (k,l)��,(k-1,l+1)��
            k1 = k-1;
            l1 = l+1;
            if k1>=-N && l1<=N
                col = (k1+N)*(2*N+1)+l1+N+1;
                H(row,col) = H(row,col)-alpha/2*exp(-1i*2*pi*f);
            end
            
            %�ǶԽ�Ԫ (k,l)��,(k+1,l-1)��
            k1 = k+1;
            l1 = l-1;
            if k1<=N && l1>=-N
                col = (k1+N)*(2*N+1)+l1+N+1;
                H(row,col) = H(row,col)-alpha/2*exp(1i*2*pi*f);
            end
            
        end
    end

    %��Hm�ı���ֵ�뱾��ʸ��
    [V_eig, d_eig] = eig(H); % d_eigΪ�Խ���, H * V_eig = V_eig * d_eig
    d_eig = d_eig*ones(length(d_eig),1); % d_eig��Ϊ������
    
    % sort������d_eig��Ԫ�ش�С�������У��������d_sort��d_eig������
    % d_sortIָʾ���к��˳��, ��d_sort=d(d_sortI)
    [d_sort, d_sortI]=sort(d_eig);   
    
    % ������ͼ����ܼ��� [����ֵ] �� [����ʸ��]
    englev(:,fN) = d_sort(1:6);  %lowest 6 levels
    ef1 = V_eig(:,d_sortI(1)); % ef1: ��E1��Ӧ�ı�������
    ef2 = V_eig(:,d_sortI(2));
    ef3 = V_eig(:,d_sortI(3));
    
    %[���������Im]
    
    %��������Imsp
    Imsp = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %��������Ԫ (k,l)��,(k-1,l+1)��
            k1 = k-1;
            l1 = l+1;
            if k1 >= -N && l1 <= N
                col = (k1+N)*(2*N+1)+l1+N+1; 
                Imsp(row,col) = Imsp(row,col) -1i/2 *exp(-1i*2*pi*f);
            end
            
            %��������Ԫ (k,l)��,(k+1,l-1)��
            k1 = k+1;
            l1 = l-1;
            if l1 >= -N && k1 <= N
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imsp(row,col) = Imsp(row,col) +1i/2 *exp(1i*2*pi*f);
            end

        end
    end
    
    %λ�Ƶ�������Imdp
    Imdp = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %��������Ԫ (k,l)��,(k-1,l)��
            k1 = k-1;
            l1 = l;
            if k1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) +1i/2;
            end
 
            %��������Ԫ (k,l)��,(k+1,l)��
            k1 = k+1;
            l1 = l;
            if k1 <= N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) -1i/2;
            end
            
            %��������Ԫ (k,l)��,(k,l-1)��
            k1 = k;
            l1 = l-1;
            if l1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) -1i/2;
            end            
            
            %��������Ԫ (k,l)��,(k,l+1)��
            k1 = k;
            l1 = l+1;
            if l1 <= N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) +1i/2;
            end

        end
    end
    
    %��ڳ���
    Imbs = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %��������Ԫ (k,l)��,(k-1,l)��
            k1 = k-1;
            l1 = l;
            if k1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imbs(row,col) = Imbs(row,col) +1i/2;
            end
 
            %��������Ԫ (k,l)��,(k+1,l)��
            k1 = k+1;
            l1 = l;
            if k1 <= N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imbs(row,col) = Imbs(row,col) -1i/2;
            end
            
        end
    end
    
    Isp(1,fN)=ef1'*Imsp*ef2;
    Isp(2,fN)=ef1'*Imsp*ef3;
    Isp(3,fN)=ef2'*Imsp*ef3;
    Isp(4,fN)=real(ef1'*Imsp*ef1);
    Isp(5,fN)=real(ef2'*Imsp*ef2);
    Isp(6,fN)=real(ef3'*Imsp*ef3);
    
    Idp(1,fN)=ef1'*Imdp*ef2;
    Idp(2,fN)=ef1'*Imdp*ef3;
    Idp(3,fN)=ef2'*Imdp*ef3;
    Idp(4,fN)=real(ef1'*Imdp*ef1);
    Idp(5,fN)=real(ef2'*Imdp*ef2);
    Idp(6,fN)=real(ef3'*Imdp*ef3);
    
    Ibs(1,fN)=ef1'*Imbs*ef2;
    Ibs(2,fN)=ef1'*Imbs*ef3;
    Ibs(3,fN)=ef2'*Imbs*ef3;
    Ibs(4,fN)=real(ef1'*Imbs*ef1);
    Ibs(5,fN)=real(ef2'*Imbs*ef2);
    Ibs(6,fN)=real(ef3'*Imbs*ef3);
    %f
end

Ia0=(Isp+Idp)*alpha/(1+2*alpha); %�����ٵĵ���Ԫ [I01 I02 I12 I00 I11 I22]
Ia0ss=Isp*alpha; %С�ڳ���
Ia0bs=Ibs; %��ڳ���


En0 = real(englev); %Hm�Ƕ��׵ģ�englevΪʵ�����������ݸ�ʽΪ�������˴�ȡʵ��ʹ���ʽ��Ϊʵ��
Ef(1,:) = englev(2,:)-englev(1,:);     %E01
Ef(2,:) = englev(3,:)-englev(1,:);     %E02
Ef(3,:) = englev(3,:)-englev(2,:);     %E12

Ef0=Ef; %�����ٵ��ܼ����










