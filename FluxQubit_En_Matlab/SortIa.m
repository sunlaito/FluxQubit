function Ia_sort=SortIa(fArray, Ia)
%Ia [I01 I02 I12 I00 I11 I22]

% [����ǶԽ�Ԫ]
Ic=Ia(1:3,:);  %[I01 I02 I12]

% �����ǶԽ�Ԫ��λ�����Ľ����Ϊ������ʵ��,��I01, I12ȡΪ��ʵ��, I02ȡΪ��ʵ��
sigIc = sum(sign(real(Ic)));
for i=1:length(fArray)
    if sigIc(i)==3 || sigIc(i)==-1
        Ic(:,i) = abs(Ic(:,i));
    
    elseif sigIc(i)==-3 || sigIc(i)==1
        Ic(:,i) = abs(Ic(:,i));
        Ic(2,i) = -Ic(2,i);
    end
end

% [���� f=0.5 ʱ�ĵ���Ԫ]

% f_c=round((fmax_g-fmin_g)/fstep_g/2+1)
f_c = find(fArray);

% ����f=0.5ʱ�ĵ����Խ�Ԫ
Ia(4:6,f_c) = 0;

% ����f=0.5ʱ�ĵ����ǶԽ�Ԫ
Ic(2,f_c) = 0;



%�����ֵ
Ia_sort=Ia;
Ia_sort(1:3,:)=Ic;

% figure();
% plot(fmin_g:fstep_g:fmax_g,Ia_sort(1:3,:));
% xlim([fmin_g,fmax_g]);

end