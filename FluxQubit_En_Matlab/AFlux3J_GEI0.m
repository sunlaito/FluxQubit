function [En0, Ef0, Ia0, Ia0ss, Ia0bs]=AFlux3J_GEI0(fArray, mac)


N = 8;  % 用于计算的基矢的数量, 取值越大越精确
mu=mac(1);    % mu = 40;     %Ej=mu*Ec
alpha=mac(2); % alpha =0.8;   %small size junction

% fArray的长度, 确定EI矩阵大小
fNmax = length(fArray);

% 给出矩阵 Hm，依照公式，逐行给Hm赋值
Hm = zeros((2*N+1)^2); %289*289

%给Hm不随f改变的部分赋值
for k = -N:N
    for l = -N:N
        row = (k+N)*(2*N+1)+l+N+1;
        
        %对角元 (k,l)行,(k,l)列
        col = row;
        Hm(row,col) = Hm(row,col)+4/mu/(1+2*alpha)*(k^2 +l^2 +alpha*(k+l)^2) +2+alpha;
        %Hm(row,col) = Hm(row,col)+2/mu *(k+l)^2 + 2/mu/(1+2*alpha) * (k-l)^2 +2+alpha;
        
        %非对角元 (k,l)行,(k-1,l)列
        k1 = k-1;    
        l1 = l;
        if k1>=-N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %非对角元 (k,l)行,(k+1,l)列
        k1 = k+1; 
        l1 = l;
        if k1<=N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %非对角元 (k,l)行,(k,l-1)列
        k1 = k;     
        l1 = l-1;
        if l1>=-N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
        
        %非对角元 (k,l)行,(k,l+1)列
        k1 = k;     
        l1 = l+1;
        if l1<=N
            col = (k1+N)*(2*N+1)+l1+N+1;
            Hm(row,col) = Hm(row,col)-0.5;
        end
      
    end
end

%对每个f值，给Hm随f变化的部分赋值
%得到完整的Hm后，算出Hm的本征值与本征矢量，并计算电流元t与I
englev=zeros(6,fNmax); %本征能级
Isp=zeros(6,fNmax); %电流I超流部分
Idp=zeros(6,fNmax); %电流I位移电流部分

Ibs=zeros(6,fNmax); %大节超流部分

for fN = 1:fNmax
    f = fArray(fN); 
    H = Hm;
    
    %给Hm随f变化的部分赋值
    for k = -N:N
        for l = -N:N
            
            row = (k+N)*(2*N+1)+l+N+1;
            
            %非对角元 (k,l)行,(k-1,l+1)列
            k1 = k-1;
            l1 = l+1;
            if k1>=-N && l1<=N
                col = (k1+N)*(2*N+1)+l1+N+1;
                H(row,col) = H(row,col)-alpha/2*exp(-1i*2*pi*f);
            end
            
            %非对角元 (k,l)行,(k+1,l-1)列
            k1 = k+1;
            l1 = l-1;
            if k1<=N && l1>=-N
                col = (k1+N)*(2*N+1)+l1+N+1;
                H(row,col) = H(row,col)-alpha/2*exp(1i*2*pi*f);
            end
            
        end
    end

    %求Hm的本征值与本征矢量
    [V_eig, d_eig] = eig(H); % d_eig为对角阵, H * V_eig = V_eig * d_eig
    d_eig = d_eig*ones(length(d_eig),1); % d_eig变为列向量
    
    % sort函数令d_eig的元素从小到大排列，结果赋给d_sort，d_eig本身不变
    % d_sortI指示排列后的顺序, 有d_sort=d(d_sortI)
    [d_sort, d_sortI]=sort(d_eig);   
    
    % 给出最低几个能级的 [本征值] 与 [本征矢量]
    englev(:,fN) = d_sort(1:6);  %lowest 6 levels
    ef1 = V_eig(:,d_sortI(1)); % ef1: 与E1对应的本征向量
    ef2 = V_eig(:,d_sortI(2));
    ef3 = V_eig(:,d_sortI(3));
    
    %[求电流矩阵Im]
    
    %超流部分Imsp
    Imsp = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %电流矩阵元 (k,l)行,(k-1,l+1)列
            k1 = k-1;
            l1 = l+1;
            if k1 >= -N && l1 <= N
                col = (k1+N)*(2*N+1)+l1+N+1; 
                Imsp(row,col) = Imsp(row,col) -1i/2 *exp(-1i*2*pi*f);
            end
            
            %电流矩阵元 (k,l)行,(k+1,l-1)列
            k1 = k+1;
            l1 = l-1;
            if l1 >= -N && k1 <= N
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imsp(row,col) = Imsp(row,col) +1i/2 *exp(1i*2*pi*f);
            end

        end
    end
    
    %位移电流部分Imdp
    Imdp = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %电流矩阵元 (k,l)行,(k-1,l)列
            k1 = k-1;
            l1 = l;
            if k1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) +1i/2;
            end
 
            %电流矩阵元 (k,l)行,(k+1,l)列
            k1 = k+1;
            l1 = l;
            if k1 <= N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) -1i/2;
            end
            
            %电流矩阵元 (k,l)行,(k,l-1)列
            k1 = k;
            l1 = l-1;
            if l1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) -1i/2;
            end            
            
            %电流矩阵元 (k,l)行,(k,l+1)列
            k1 = k;
            l1 = l+1;
            if l1 <= N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imdp(row,col) = Imdp(row,col) +1i/2;
            end

        end
    end
    
    %大节超流
    Imbs = zeros((2*N+1)^2);
    for k = -N:N
        for l = -N:N
            row = (k+N)*(2*N+1)+l+N+1;
            
            %电流矩阵元 (k,l)行,(k-1,l)列
            k1 = k-1;
            l1 = l;
            if k1 >= -N 
                col = (k1+N)*(2*N+1)+l1+N+1;
                Imbs(row,col) = Imbs(row,col) +1i/2;
            end
 
            %电流矩阵元 (k,l)行,(k+1,l)列
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

Ia0=(Isp+Idp)*alpha/(1+2*alpha); %无量纲的电流元 [I01 I02 I12 I00 I11 I22]
Ia0ss=Isp*alpha; %小节超流
Ia0bs=Ibs; %大节超流


En0 = real(englev); %Hm是厄米的，englev为实数，但其数据格式为复数，此处取实部使其格式变为实数
Ef(1,:) = englev(2,:)-englev(1,:);     %E01
Ef(2,:) = englev(3,:)-englev(1,:);     %E02
Ef(3,:) = englev(3,:)-englev(2,:);     %E12

Ef0=Ef; %无量纲的能级间距










