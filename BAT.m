%% 蝙蝠算法实现
% 目标函数还是老样子，下面一部分是图形的绘制，不用管
% x = -100:100;
% [X,Y] = meshgrid(x);
% Z = (X - 50).^2 + (Y - 20).^2;
% mesh(X,Y,Z)
clc;
clear;
global x_a x_b
x_a = 0;
x_b = 0;
%% 交互式界面初始化
figure(1)
plotbutton = uicontrol('style','pushbutton','string','运行','fontsize',12,'position',[10,400,50,20],'callback','run = 1;');
stop = 0;run = 0; % 初始化界面标志
while stop == 0
    if run == 1
        %% 参数初始化
        %         wmax = 2.8;%惯性权重最大值
        %         wmin = 0.4;%惯性权重最小值
        X_max = 100; % 定义自变量的上界
        X_min = -100; % 定义自变量的下界
        alph = 0.90; % 脉冲响度衰减系数
        gamma = 0.05; % 脉冲频率增强系数
        Q_min = -1; % 频率下界
        Q_max = 1; % 频率上界
        A_min = 1; % 响度下界
        A_max = 2; % 响度上界
        MAX_Iter = 100; % 迭代总次数
        Iter = 1; % 初始迭代次数
        N = 20; % 种群大小
        % 初始化种群
        x =  zeros(N,2); % 初始化种群
        fitt = zeros(N,1); % 初始化适应度
        v = zeros(N,2); % 初始化速度
        Q = zeros(N,1); % 初始化频率
        X = zeros(N,2);
        A = 1+1*ones(N,1);
        %         R = 0.5*ones(N,1);
        R = 0.5;
        FIT_BEST = []; % 定义中间变量
        BEST = []; % 定义中间变量
        for i = 1:N
            x(i,1) = X_min + rand()*(X_max - X_min);
            x(i,2) = X_min + rand()*(X_max - X_min);
            fitt(i) = fit(x(i,:));
        end
        [fit_best,ind] = min(fitt); % 找到最佳适应度和其位置下标
        X_best = x(ind,:); %读取最佳位置
        %% 迭代过程
        while Iter <= MAX_Iter
            Iter = Iter + 1;
            for i = 1:N
                Q = Q_min + rand * (Q_max - Q_min); % 更新频率
                %                 w=(wmax-wmin)*exp(-2*(Iter / MAX_Iter)^2)+wmin;%惯性权重因子
                v(i,:) =  v(i ,:) - (x(i,:) - X_best) * Q; % 更新速度
                X(i,:) = x(i,:) + v(i,:); % 更新位置
                %判断位置是否越界
                if X(i,1) < X_min || X(i,1) > X_max % 判断是否超出边界
                    X(i,1) = X_min + rand(1) * (X_max - X_min);
                end
                if X(i,2) < X_min || X(i,2) > X_max
                    X(i,2) = X_min + rand(1) * (X_max - X_min);
                end
                if  rand > R % 如果大于R则进行随机搜索
                    rand_ind = unidrnd(N);
                    %                     X(i,:) = x(rand_ind,:) + unifrnd(-1,1) * A(i);
                    %                     X(i,:) = X_best + unifrnd(-1,1) * mean(A);
                    %                     X(i,:) = X_best + unifrnd(-1,1) *  randn(1,2);
                    X(i,:) = x(rand_ind,:) + unifrnd(-1,1) * randn(1,2);
                    if X(i,1) < X_min || X(i,1) > X_max % 判断是否超出边界
                        X(i,1) = X_min + rand(1) * ( X_max - X_min );
                    end
                    if X(i,2) < X_min || X(i,2) > X_max
                        X(i,2) = X_min + rand(1) * ( X_max - X_min );
                    end
                end
                Fnew = fit(X(i,:)); % 计算当前位置的适应度
                if (Fnew <= fitt(i))
                    x(i,:) = X(i,:);
                    fitt(i) = Fnew;  % 如果满足条件，就将位置更新
                end
                % 判断当前适应度与最优适应度之间的关系
                if (Fnew <= fit_best)
                    X_best = X(i,:);
                    fit_best = Fnew;
                    A(i) = alph * A(i);
                    R = R * ( 1 - exp( - (gamma * Iter)) ); % R在不断减小，这样有利于后期更多的游动
                end
                
            end
            scatter(x(:,1),x(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(Iter)])
            pause(0.05)
            [fit_b,ind_b] = min(fitt); % 找到一次迭代后的最优解
            x_B = x(ind_b,:);
            FIT_BEST = [FIT_BEST;fit_b]; % 记录最优解
            BEST = [BEST;x_B]; % 记录最优位置
        end
        [Fit_Best,ind_be] = min(FIT_BEST);
        X_BEST = BEST(ind_be,:);
        disp(['找到的最优解为：',num2str(Fit_Best)])
        disp(['最优解对应的X：',num2str(X_BEST)])
        figure(2)
        plot(1:MAX_Iter,FIT_BEST)
        title('种群迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        stop = 1;
    end
    drawnow
end
% 改进方向：
% 1、适应性权重
% 2、随机移动的参照物选取
% 3、
X = [1,2,3,4,5,6];
Y = [4,5,6,7,1,9];
corrcoef(X',Y')