
%% 模拟退火算法实现
clc;clear;
global x_a x_b
x_a = 0;
x_b = 0;
%% 参数初始化
T = 1000;   % 初始温度
T_min = 1; % 温度阈值
Markov  = 100;  % 每个温度下的迭代次数
alfa = 0.98;  % 温度衰减系数
X_max = 100; % 自变量的上界
X_min = -100; % 自变量的下界
x = rand(1,2) .* (X_max - X_min) + X_min; 
% x(1,1) = X_min + rand * (X_max - X_min); % 随机生成一个解
% x(1,2) = X_min + rand * (X_max - X_min); % 随机生成一个解
fitt = fit(x); % 计算其适应度
fitt_min = fitt; % 将初始适应度设为最小
x_min = x; % 将初始位置记为最小位置
FITT = []; % 记录中间变量
X = []; % 记录中间变量
%% 迭代开始
tic
while T > T_min
    for i = 1:Markov
        x_new(1,1) = X_min + rand * (X_max - X_min); % 生成新解
        x_new(1,2) = X_min + rand * (X_max - X_min); % 生成新解
        fitt_new = fit(x_new); % 计算新解的适应度
        if fitt_new < fitt % 因为是找最小值，如果新解小于最小解，则更新解
            fitt = fitt_new;
            x = x_new;
        else % 如果不小于，则随机选择
            p = exp(-(fitt_new - fitt)/T); % 接受概率
            if rand < p % 如果随机数小于随机概率，则接受新解
                fitt = fitt_new;
                x = x_new;
            end
        end
        if fitt < fitt_min
            fitt_min = fitt;
            x_min = x_new;
        end
    end
    % 绘制过程图
    scatter(x_min(1,1),x_min(1,2),'.','r'); % 绘制散点图
    hold on
    scatter(x_a,x_b,'.','b')
    hold off
    title(['当前温度为：',num2str(T)])
    axis([-100,100,-100,100])
    pause(0.1)
    FITT = [FITT fitt_min];
    X = [X;x_min];
    T = T * alfa; % 温度下降
end
toc
%% 结果可视化
[fit_best,ind_best] = min(FITT);
x_best = X(ind_best,:);
disp(['模拟退火算法求解的最优解为：',num2str(fit_best)])
disp(['模拟退火算法求解的最优位置为：',num2str(x_best)])

figure(2)
plot(FITT);
xlabel('迭代次数');
ylabel('最小时间');
