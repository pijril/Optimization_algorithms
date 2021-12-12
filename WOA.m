%% 鲸鱼算法
clc;clear
a = 0;
b = 0;
y = (x(1) - a)^2 + (x(2) - b)^2;
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
        a = 2; % 初始游动因子的值，后续会不断缩减
        MAX_iter = 100; % 迭代次数
        iter = 1; % 初始化迭代次数
        b = 3; % 螺旋系数
        number = 20; % 定义鲸鱼的数量
        fitt = zeros(number,1); % 初始化适应度
        x = zeros(number,2); % 初始化种群
        l_min = -100; %自变量下界
        l_max = 100; % 自变量上界
        for i = 1:number % 初始化种群和适应度
            x(i,1) = l_min + rand(1) * (l_max - l_min);
            x(i,2) = l_min + rand(1) * (l_max - l_min);
            fitt(i) = fit(x(i,:));
        end
        fit_best = min(fitt); % 记录最优适应度
        ind = find(fitt == min(fitt),1);
        X_best = x(ind,:); % 记录最优位置
        X_BEST = []; % 记录全局
        FITT = []; % 记录全局
        %% 迭代开始
        while iter <= MAX_iter
            for i = 1:number
                if rand() < 0.5 % 采用游动捕食
                    a = 2 * (1 - iter / MAX_iter); % 更新游动因子
                    r = rand(1,2); % 随机数
                    C = 2 * r; % 计算摆动因子
                    A =  2 * r * a - a; % 计算收敛因子
                    if sqrt(A(1)^2 + A(2)^2) < 1 % 采取包围捕食
                        D = abs(C .* X_best - x(i,:)); % 计算最优距离
                        x(i,:) = X_best - A .* D; % 这里是位置更新
                        if x(i,1) < l_min || x(i,1) > l_max % 判断是否超出边界
                            x(i,1) = l_min + rand(1) * (l_max - l_min);
                        end
                        if x(i,2) < l_min || x(i,2) > l_max
                            x(i,2) = l_min + rand(1) * (l_max - l_min);
                        end
                    else % 采取随机游动
                        rand_ind = unidrnd(20); % 随机抽取鲸鱼
                        D = abs(C .* x(rand_ind,:) - x(i,:));
                        x(i,:) = x(rand_ind,:) - A .* D;
                        if x(i,1) < l_min || x(i,1) > l_max % 判断是否超出边界
                            x(i,1) = l_min + rand(1) * (l_max - l_min);
                        end
                        if x(i,2) < l_min || x(i,2) > l_max
                            x(i,2) = l_min + rand(1) * (l_max - l_min);
                        end
                    end
                else % 采用气泡捕食
                    D_1 = abs(X_best - x(i,:));
                    l = unifrnd(-1,1); % 随机数
                    x(i,:) = D_1 * exp(b * l) * cos(2 * pi * l) + X_best;
                    if x(i,1) < l_min || x(i,1) > l_max % 判断是否超出边界
                        x(i,1) = l_min + rand(1) * (l_max - l_min);
                    end
                    if x(i,2) < l_min || x(i,2) > l_max
                        x(i,2) = l_min + rand(1) * (l_max - l_min);
                    end
                end
                fitt(i) = fit(x(i,:)); % 每更新一次位置，求解一次适应度
            end
            % 绘制迭代过程图
            scatter(x(:,1),x(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(iter)])
            pause(0.07)
            %     disp(['当前迭代次数为：',num2str(iter)])
            % 迭代完一次更新数据
            fit_best = min(fitt); % 更新当前代数中种群的最优解
            ind = find(fitt == min(fitt),1);
            X_best = x(ind,:); % 更新最优位置
            X_BEST = [X_BEST;X_best];
            FITT = [FITT;fit_best];
            iter = iter + 1;
        end
        Fit_Best = min(FITT);
        ind_best = find(FITT == min(FITT),1);
        X_BE = X_BEST(ind_best,:);
        disp(['找到的最优解为：',num2str(Fit_Best)])
        disp(['最优解对应的X：',num2str(X_BE)])
        figure(2)
        plot(1:MAX_iter,FITT)
        title('种群迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        stop = 1;
    end
    drawnow
end