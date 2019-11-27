clc;
clear all;
close all;

load('multiPath.txt');
% read data
length_all = multiPath(end,1);% 1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  get the first data (ref data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_ref = [];
AMP_ref = [];
ANG_ref = [];
P_ref = [];
time_ref = [];
for r = 1:size(multiPath,1)  % rows
    % get the ref (1)
    if multiPath(r,1) == 1
        T_ref =[T_ref multiPath(r,2)];
        AMP_ref =[AMP_ref multiPath(r,3)];
        ANG_ref =[ANG_ref multiPath(r,4)];
        P_ref =[P_ref multiPath(r,5)];
        time_ref =[time_ref multiPath(r,6)];
    else
        path_num = length(T_ref);  % path number
        now_row = multiPath(r,1);
        num = 0;
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the space for "save data"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_path = cell(1,path_num);
for i = 1:path_num
    final_path{i} = zeros(length_all,5);
end

% save first data
for p = 1:path_num % first path
    final_path{p}(1,1) = T_ref(p);
    final_path{p}(1,2) = AMP_ref(p);
    final_path{p}(1,3) = ANG_ref(p);
    final_path{p}(1,4) = P_ref(p);
    final_path{p}(1,5) = time_ref(p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search all data and 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = path_num+1:size(multiPath,1)  % rows
    % get the ref (1)
    if multiPath(r,1)== now_row
        T(num+1) = multiPath(r,2);
        AMP(num+1) = multiPath(r,3);
        ANG(num+1) = multiPath(r,4);
        P(num+1) = multiPath(r,5);
        Time(num+1) = multiPath(r,6);
        num = num +1;
        if now_row == length_all
            ind_T0 = 0;
            for p = 1:path_num % first path
                for k = 1:length(T)
                    T_diff(k) = T_ref(p) - T(k);
                    Amp_diff(k) = AMP_ref(p) - AMP(k);
                    Ang_diff(k) = ANG_ref(p) - ANG(k);
                end
                [val_T,ind_T]=min(abs(T_diff));
                if ind_T == ind_T0
                    T_diff(ind_T) = 100;
                    [val_T,ind_T]=min(abs(T_diff));
                else
                    ind_T = ind_T;
                end
                if abs(T(ind_T)-T_ref(p))<1 %&& abs(ANG_ref(p) - ANG(ind_T))<2
                    final_path{p}(now_row,1) = T(ind_T);
                    final_path{p}(now_row,2) = AMP(ind_T);
                    final_path{p}(now_row,3) = ANG(ind_T);
                    final_path{p}(now_row,4) = P(ind_T);
                    final_path{p}(now_row,5) = Time(ind_T);
                    ind_T0 = ind_T;
                else
                    final_path{p}(now_row,1) = 0;
                    final_path{p}(now_row,2) = 0;
                    final_path{p}(now_row,3) = 0;
                    final_path{p}(now_row,4) = 0;
                    final_path{p}(now_row,5) = 0;
                    ind_T0 = ind_T;
                end
            end
        end
    else
        num=0;
        ind_T0 = 0;
        for p = 1:path_num % first path
            for k = 1:length(T)
                T_diff(k) = T_ref(p) - T(k);
            end
            [val_T,ind_T]=min(abs(T_diff));
            if ind_T == ind_T0
                T_diff(ind_T) = 100;
                [val_T,ind_T]=min(abs(T_diff));
            else
                ind_T=ind_T;
            end
            
            if abs(T(ind_T)-T_ref(p))<1 %&& abs(ANG_ref(p) - ANG(ind_T))<2
                final_path{p}(now_row,1) = T(ind_T);
                final_path{p}(now_row,2) = AMP(ind_T);
                final_path{p}(now_row,3) = ANG(ind_T);
                final_path{p}(now_row,4) = P(ind_T);
                final_path{p}(now_row,5) = Time(ind_T);
                ind_T0 = ind_T;
            else
                final_path{p}(now_row,1) = 0;
                final_path{p}(now_row,2) = 0;
                final_path{p}(now_row,3) = 0;
                final_path{p}(now_row,4) = 0;
                final_path{p}(now_row,5) = 0;
                ind_T0 = ind_T;
            end

        end
        
        now_row = now_row +1; 
        
        T(num+1) = multiPath(r,2);
        AMP(num+1) = multiPath(r,3);
        ANG(num+1) = multiPath(r,4);
        P(num+1) = multiPath(r,5);
        Time(num+1) = multiPath(r,6);
        num = num +1;
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
px = 1:length_all;
pn = length(px);
orign_delay = 0;
py1 = orign_delay + final_path{1}(:,1);
py2 = orign_delay + final_path{2}(:,1);
py3 = orign_delay + final_path{3}(:,1);

f990=figure(990)
whitebg(f990,[0 .5 .6]);
hold on
for i = 1:pn
   plot(px(i),py1(i),'r-p',px(i),py2(i),'g-p',px(i),py3(i),'y-p');
%    stem(px(i),py1(i),'r-p');
%    stem(px(i),py2(i),'g-p');
%    stem(px(i),py3(i),'b-p');
%    axis([1 length_all 0 8])
   M(i) = getframe;
   pause(0.001);
end
% movie(M);


f991 = figure(991)
whitebg(f991,[0 .5 .6]);
subplot(3,2,1)
plot(final_path{1, 1}(:,1),'r')
hold on
plot(final_path{1, 2}(:,1),'g')
hold on
plot(final_path{1, 3}(:,1),'y')
title("Tau")
subplot(3,2,2)
plot(diff(final_path{1, 1}(:,1)),'r')
hold on
plot(diff(final_path{1, 2}(:,1)),'g')
hold on
plot(diff(final_path{1, 3}(:,1)),'y')
title("Tau diff")


subplot(3,2,3)
plot(unwrap(final_path{1, 1}(:,3)),'r')
hold on 
plot(unwrap(final_path{1, 2}(:,3)),'g')
plot(unwrap(final_path{1, 3}(:,3)),'y')
title("Angle")
subplot(3,2,4)
plot(unwrap(diff(final_path{1, 1}(:,3))),'r')
hold on
plot(unwrap(diff(final_path{1, 2}(:,3))),'g')
plot(unwrap(diff(final_path{1, 3}(:,3))),'y')
title("Angle diff")

subplot(3,2,5)
plot(final_path{1, 1}(:,2),'r')
hold on 
plot(final_path{1, 2}(:,2),'g')
hold on 
plot(final_path{1, 3}(:,2),'y')
title("Amp")
subplot(3,2,6)
plot(diff(final_path{1, 1}(:,2)),'r')
hold on
plot(diff(final_path{1, 2}(:,2)),'g')
hold on
plot(diff(final_path{1, 3}(:,2)),'y')
title("Amp diff")