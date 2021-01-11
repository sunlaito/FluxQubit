function Ia_sort=SortIa(fArray, Ia)
%Ia [I01 I02 I12 I00 I11 I22]

% [处理非对角元]
Ic=Ia(1:3,:);  %[I01 I02 I12]

% 电流非对角元相位调整的结果或为三项正实数,或I01, I12取为正实数, I02取为负实数
sigIc = sum(sign(real(Ic)));
for i=1:length(fArray)
    if sigIc(i)==3 || sigIc(i)==-1
        Ic(:,i) = abs(Ic(:,i));
    
    elseif sigIc(i)==-3 || sigIc(i)==1
        Ic(:,i) = abs(Ic(:,i));
        Ic(2,i) = -Ic(2,i);
    end
end

% [处理 f=0.5 时的电流元]

% f_c=round((fmax_g-fmin_g)/fstep_g/2+1)
f_c = find(fArray);

% 处理f=0.5时的电流对角元
Ia(4:6,f_c) = 0;

% 处理f=0.5时的电流非对角元
Ic(2,f_c) = 0;



%结果赋值
Ia_sort=Ia;
Ia_sort(1:3,:)=Ic;

% figure();
% plot(fmin_g:fstep_g:fmax_g,Ia_sort(1:3,:));
% xlim([fmin_g,fmax_g]);

end