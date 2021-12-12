% 实现人工鱼群算法
clc,clear;
global x_a x_b % 全局变量
x_a = 0;x_b = 0; % 控制最优解的位置

%% 交互式界面初始化
figure(1)
plotbutton = uicontrol('style','pushbutton','string','运行','fontsize',12,'position',[10,400,50,20],'callback','run = 1;');
stop = 0;run = 0; % 初始化界面标志
while stop == 0
    if run == 1
        %% 参数初始化
        tic % 计算算法运行时间
        Visual = 30; % 视野范围
        Step = 3; % 移动步长
        N = 50; % 种群大小
        D = 2; % 解空间的维数
        Try_number =10; % 觅食时的最大迭代次数
        delta = 7; % 拥挤因子
        x_min = -100; % 自变量下界
        x_max = 100; % 自变量上界
        Iter = 0; % 初始迭代次数
        Iter_max = 100; % 最大迭代次数
        
        %% 种群初始化
        x = rand(N,D) * (x_max - x_min) + x_min;
        % 计算初始种群的适应度
        for i = 1:N
            fitness(i) = fit(x(i,:));
        end
        % 挑选最优适应度并记录全程
        [best_fitness,I] = min(fitness);
        best_x = x(I,:);
        BEST_Iter_X = []; % 定义中间变量
        BEST_Iter_fitness = [];
        %% 迭代开始
        while Iter < Iter_max
            Iter = Iter + 1; % 迭代次数加1
            Evolution_flag = 0; % 增加一个进化条件
            % 遍历每一只人工鱼
            for i = 1:N
                %% 聚群行为(Cluster behavior)
                nf = 0; % 初始化第i只鱼的聚群数量为0
                x_Cluster = [0,0]; % 定义一个累加的变量，用于储存聚群鱼的位置
                for j=1:N % 遍历每一只人工鱼
                    if i~=j
                        if (norm(x(i,:) - x(j,:)) < Visual) % 如果这只鱼在第i只鱼的视野内
                            nf = nf + 1; % 聚群鱼加1
                            x_Cluster = x_Cluster + x(j,:); % 将视野内的鱼进行累加
                        end
                    end
                end
%                 x_Cluster = x_Cluster - x(i,:); % 上面遍历j时，将j=i的鱼也计入其中了，所以需要减去i鱼
%                 nf = nf - 1; % 同样数量减1
                x_c = x_Cluster./nf; % 计算聚群的中心位置
                if (fit(x_c)*nf) < delta*fit(x(i,:)) % 如果满足条件
                    X_next1 = x(i,:) + rand() * Step * ((x_c-x(i,:))/norm(x_c - x(i,:))); % 发生移动
                    % 判断是否越界
                    for k=1:D
                        if X_next1(k) > x_max || X_next1(k) < x_min
                            X_next1(k) = rand() * (x_max - x_min) + x_min;
                        end
                    end
                    % 计算移动后的适应度
                    X_next1_fitness = fit(X_next1); % 计算适应度
                else % 进行觅食
                    flag = 0; % 如果在觅食次数内发生更新的标志
                    for j = 1:Try_number
                        x_new = x(i,:) + unifrnd(-1,1) * Visual; % 在第i只鱼的视野内随机选择一个鱼
                        % 判断是否越界
                        for k=1:D
                            if x_new(k) > x_max || x_new(k) < x_min
                                x_new(k) = rand() * (x_max - x_min) + x_min;
                            end
                        end
                        if fit(x(i,:)) > fit(x_new) % 如果选择的鱼比第i只鱼要好，则第i只鱼向选择的鱼移动
                            X_next1 = x(i,:) + rand() * Step * ((x_new - x(i,:))/norm(x_new - x(i,:)));
                            flag = 1; % 更新标志
                            % 判断是否越界
                            for k=1:D
                                if X_next1(k) > x_max || X_next1(k) < x_min
                                    X_next1(k) = rand() * (x_max - x_min) + x_min;
                                end
                            end
                            X_next1_fitness = fit(X_next1); % 计算觅食后的适应度
                            break; % 已经更新位置，跳出循环
                        end
                    end
                    if flag ==0 % 如果没有觅食成功
                        X_next1 = x(i,:) + unifrnd(-1,1) * Step; % 如果没有更新，则采取随机移动的方式更新
                        % 判断是否越界
                        for k=1:D
                            if X_next1(k) > x_max || X_next1(k) < x_min
                                X_next1(k) = rand() * (x_max - x_min) + x_min;
                            end
                        end
                        X_next1_fitness = fit(X_next1); % 计算适应度
                    end
                end
                
                %% 追尾行为（This behavior）
                fitness_follow = inf; % 定义初始的适应度
                best_flag = 0; % 判断视野内是否有鱼伙伴的标志，如果有，则将其置为1
                best_follow = []; % 定义视野内最优的鱼
                for j = 1:N
                    if (i~=j) % 防止与自己追尾
                        if (norm(x(j,:) - x(i,:)) < Visual) && (fit(x(j,:)) < fitness_follow) % 判断视野内最优解
                            best_follow = x(j,:); % 更新视野内最优鱼
                            fitness_follow = fit(x(j,:));
                            best_flag = 1; % 如果有一只鱼在其视野内，则其值变为1
                        end
                    end
                end
                nf_follow = 0; % 视野内最优鱼的周围鱼的个数
                for j=1:N % 遍历所有鱼
                    if i~=j % 视野内最优鱼的周围不包括其本身这只鱼
                        if best_flag == 1 % 判断是否拥挤的前提是视野内存在鱼
                            if (norm(x(j,:) - best_follow) < Visual) % 判断是否拥挤
                                nf_follow = nf_follow + 1;
                            end
                        end
                    end
                end
                % 进行判断是觅食移动还是随机移动
                if (best_flag == 1)&&fit(best_follow)/nf_follow < delta*fit(x(i,:)) % 如果满足条件，则该鱼向追尾鱼移动
                    X_next2 = x(i,:) + rand() * Step * (best_follow - x(i,:))/norm(best_follow - x(i,:));
                    % 判断是否越界
                    for k=1:D
                        if X_next2(k) > x_max || X_next2(k) < x_min
                            X_next2(k) = rand() * (x_max - x_min) + x_min;
                        end
                    end
                    X_next2_fitness = fit(X_next2); % 记录适应度
                else % 如果不满足条件，则发生觅食行为
                    flag = 0; % 如果在觅食次数内发生更新的标志
                    for j = 1:Try_number
                        x_new = x(i,:) + unifrnd(-1,1) * Visual; % 在第i只鱼的视野内随机选择一个鱼
                        % 判断是否越界
                        for k=1:D
                            if x_new(k) > x_max || x_new(k) < x_min
                                x_new(k) = rand() * (x_max - x_min) + x_min;
                            end
                        end
                        if fit(x(i,:)) > fit(x_new) % 如果选择的鱼比第i只鱼要好，则第i只鱼向选择的鱼移动
                            X_next2 = x(i,:) + rand() * Step * ((x_new - x(i,:))/norm(x_new - x(i,:)));
                            flag = 1; % 更新标志
                            % 判断是否越界
                            for k=1:D
                                if X_next2(k) > x_max || X_next2(k) < x_min
                                    X_next2(k) = rand() * (x_max - x_min) + x_min;
                                end
                            end
                            X_next2_fitness = fit(X_next2); % 计算觅食后的适应度
                            break; % 已经更新位置，跳出循环
                        end
                    end
                    if flag == 0 % 如果没有觅食成功
                        X_next2 = x(i,:) + unifrnd(-1,1) * Step; % 如果没有更新，则采取随机移动的方式更新
                        % 判断是否越界
                        for k=1:D
                            if X_next2(k) > x_max || X_next2(k) < x_min
                                X_next2(k) = rand() * (x_max - x_min) + x_min;
                            end
                        end
                        X_next2_fitness = fit(X_next2); % 计算适应度
                    end
                end
                % 第i只鱼更新情况
                if X_next1_fitness < X_next2_fitness
                    x(i,:) = X_next1; % 如果聚群最优，则鱼更新
                else
                    x(i,:) = X_next2; % 如果追尾最优，则鱼更新
                end
            end
            % 一轮迭代之后，更新信息
            for i =1:N
                if fit(x(i,:)) < best_fitness
                    best_x = x(i,:); % 更新种群中的最优鱼
                    best_fitness = fit(x(i,:)); % 更新最优适应度
                end
            end
            % 绘制迭代过程图
            scatter(x(:,1),x(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(Iter)])
            pause(0.01)
            BEST_Iter_X = [ BEST_Iter_X;best_x]; % 定义中间变量
            BEST_Iter_fitness = [BEST_Iter_fitness;best_fitness];
        end
        % 输出结果
        [Fit_Best,II] = min(BEST_Iter_fitness);
        X_BE = BEST_Iter_X(II,:);
        disp(['找到的最优解为：',num2str(Fit_Best)])
        disp(['最优解对应的X：',num2str(X_BE)])
        figure(2)
        plot(1:Iter_max,BEST_Iter_fitness)
        title('种群迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        toc
        stop = 1;
    end
    drawnow
end