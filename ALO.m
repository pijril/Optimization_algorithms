clc;clear
global x_a x_b
x_a = 48;
x_b = -48;
%% 交互式界面初始化
figure(1)
plotbutton = uicontrol('style','pushbutton','string','运行','fontsize',12,'position',[10,400,50,20],'callback','run = 1;');
stop = 0;run = 0; % 初始化界面标志
while stop == 0
    if run == 1
        %% 参数初始化
        ant_N = 20; % 蚂蚁的数量
        antlion_N = 20; % 蚁狮的数量
        D = 2; % 问题的维数
        X_max = [100,100]; % 上界
        X_min = [-100,-100]; % 下界
        Iter = 0; % 初始迭代次数
        Iter_max = 100; % 最大迭代次数
        ant_X = rand(ant_N,D) .* (X_max - X_min) + X_min; % 初始化蚂蚁种群
        antlion_X = rand(ant_N,D) .* (X_max - X_min) + X_min; % 初始化蚁狮种群
        ant_fitness = zeros(ant_N,1); % 初始化蚂蚁适应度
        antlion_fitness = zeros(antlion_N,1); % 初始化蚁狮适应度
        % 计算适应度
        for i = 1:ant_N
            ant_fitness(i) = fit(ant_X(i,:)); % 计算蚂蚁适应度
            antlion_fitness(i) = fit(antlion_X(i,:)); % 计算蚁狮适应度
        end
        % 存储最优的蚁狮
        antlion_best_fitness = min(antlion_fitness);
        antlion_best_X = antlion_X(antlion_fitness == antlion_best_fitness,:);
        
        TEMP = []; % 记录中间过程
        TEMP_X = [];
        % 迭代过程
        while Iter < Iter_max
            TEMP = [TEMP;antlion_best_fitness];
            TEMP_X = [TEMP_X;antlion_best_X];
            Iter = Iter + 1;
            for i = 1:ant_N % 遍历所有蚂蚁
                % 为蚂蚁挑选合适的蚁狮,采用轮盘法
%                 antlion_probability = antlion_fitness./sum(antlion_fitness);
%                 antlion_probability_sort = cumsum(antlion_probability,1);
%                 % antlion_probability_X = [antlion_probability,antlion_X];
%                 r = rand;
%                 tar = find(antlion_probability_sort >= r);
%                 chose = tar(1);
%                 antlion_current = antlion_X(chose,:); % 假设挑选最好的蚁狮
                % % 随机为蚂蚁选择合适的蚁狮
                chose = unidrnd(ant_N);
                antlion_current = antlion_X(chose,:);
%                 更新c和d
%                 先计算陷阱的范围
                if Iter < 0.1 * Iter_max
                    w = 2;
                elseif (0.1 * Iter_max <= Iter) && (Iter < 0.5 * Iter_max)
                    w = 2;
                elseif (0.5 * Iter_max <= Iter) && (Iter < 0.75 * Iter_max)
                    w = 3;
                elseif (0.75 * Iter_max <= Iter) && (Iter < 0.9 * Iter_max)
                    w = 4;
                elseif (0.9 * Iter_max <= Iter) && (Iter < 0.95 * Iter_max)
                    w = 5;
                elseif (0.95 * Iter_max <= Iter) && (Iter <= Iter_max)
                    w = 6;
                end
                I = 10^w * Iter / Iter_max;
                % 创建一个随机游走
                % 对每一维数进行更新,注意这里需要运行D次
                for q = 1:D
                    c = X_min(q) / I + antlion_current(q); % 计算此维数的陷阱的下界
                    d = X_max(q) / I + antlion_current(q); % 计算此维数的陷阱的上界
                    Random_walk = zeros(Iter_max,1);
                    Random_walk(i) = 1;
                    for j = 2:Iter_max % 遍历所有迭代次数
                        for k = 1:j
                            if rand < 0.5
                                temp = -1;
                            else
                                temp = 1;
                            end
                            Random_walk(j) = Random_walk(j - 1) + temp;
                        end
                    end
                    % 更新蚂蚁位置
                    ARW = (Random_walk(i) - min(Random_walk) + 1) / (max(Random_walk) - min(Random_walk) + 1);
                    ant_temp(1,q) =  ARW * (d - c) + c; % 计算该维数的变化值
                end
                % 蚂蚁精英化
                ant_X(i,:) = (ant_temp + antlion_best_X) / 2;
                % 判断蚁狮是否捕食蚂蚁
                if fit(ant_X(i,:)) < fit(antlion_current) % 如果蚂蚁比蚁狮要好，则更新
                    antlion_X(chose,:) = ant_X(i,:); % 更新对应这个蚂蚁的蚁狮
                    antlion_fitness(chose) = fit(ant_X(i,:));
                end
            end
            % 更新最优
            antlion_best_fitness = min(antlion_fitness);
            antlion_best_X = antlion_X(find(antlion_fitness == antlion_best_fitness,1),:);
            scatter(antlion_X(:,1),antlion_X(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(Iter)])
            pause(0.05)
        end
        fit_best = min(TEMP);
        
        fit_best_x = TEMP_X(TEMP == fit_best,:);
        disp(['最优的结果为：' num2str(fit_best)])
        disp(['最优的结果为：' num2str(fit_best_x(1,:))])
        plot(TEMP)
        title('算法迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')    
        stop = 1;
    end
    drawnow
end