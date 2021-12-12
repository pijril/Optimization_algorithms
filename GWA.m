%% 灰狼算法实现
clc;
clear;
rng('shuffle')
global x_a x_b
x_a = 5;
x_b = -5;
%% 交互式界面初始化
figure(1)
plotbutton = uicontrol('style','pushbutton','string','运行','fontsize',12,'position',[10,400,50,20],'callback','run = 1;');
stop = 0;run = 0; % 初始化界面标志
while stop == 0
    if run == 1
        %% 参数初始化
        N = 20; % 种群的个数
        D = 2; % 维数
        X_max = 100;
        X_min = -100;
        a = 2; % 收敛因子
        MAX_Iter = 100; % 定义迭代次数
        x = zeros(N,D); % 分配内存
        fitt = zeros(N,1); % 分配内存
        Iter = 0; % 迭代次数
        for i = 1:N
            x(i,1) = X_min + rand()*(X_max - X_min);
            x(i,2) = X_min + rand()*(X_max - X_min);
            fitt(i) = fit(x(i,:));
        end
        % 选出三头狼王，用X_First、X_Second、X_Third表示；
        temp = [fitt x];
        temp_sort = sortrows(temp,1); % 对重组矩阵按第一列的进行升序
        X_First = temp_sort(1,2:end); % 选出狼王
        X_Second = temp_sort(2,2:end); % 选出次狼
        X_Third = temp_sort(3,2:end);  % 选出次次狼
        % [fit_best,ind] = min(fitt); % 找到最佳适应度和其位置下标
        % X_best = x(ind,:); %读取最佳位置
        FIT_BEST = []; % 记录每次迭代过程中的最优适应度
        X_BEST = []; % 记录一下每次迭代过程中的最优适应度对应的位置
        
        %% 迭代开始
        while Iter < MAX_Iter
            Iter = Iter + 1;
            a = 2 - 2 * (Iter/MAX_Iter); % 更新收敛因子（这里是线性收敛因子)
            %             a = 2 - 2 * ((1 / (exp(1) - 1)) * (exp(Iter / MAX_Iter) - 1)); % 非线性收敛因子
            for i = 1:N
                % 分别计算第i只狼与三只头狼的协同系数C与A
                %                 A_1 = 2 * a * rand() - a; % 计算协同系数
                A_1 = unifrnd(-a,a,1,2);
                C_1 = 2 * rand(1,2);  % 计算协同系数
                %                 A_2 = 2 * a * rand() - a; % 计算协同系数
                A_2 = unifrnd(-a,a,1,2);
                
                C_2 = 2 * rand(1,2);  % 计算协同系数
                %                 A_3 = 2 * a * rand() - a; % 计算协同系数
                A_3 = unifrnd(-a,a,1,2);
                
                C_3 = 2 * rand(1,2);  % 计算协同系数
                % 计算第i只狼与三只头狼的直接距离
                D_1 = abs(C_1 .* X_First - x(i,:));
                D_2 = abs(C_2 .* X_Second - x(i,:));
                D_3 = abs(C_3 .* X_Third - x(i,:));
                % 更新位置公式
                X_1 = X_First - A_1 .* D_1;
                X_2 = X_Second - A_2 .* D_2;
                X_3 = X_Third - A_3 .* D_3;
                %更新位置，这里给三只狼王加上一个贪心准则
                if (sum((x(i,:) ~= X_First)) == 2)&&(sum((x(i,:) ~= X_Second)) == 2) && (sum((x(i,:) ~= X_Third)) == 2)
                    x(i,:) = (X_1 + X_2 + X_3) / 3;
                else
                    X = (X_1 + X_2 + X_3) / 3;
                    Fem = fit(X(:)); % 计算三只狼王如果更新位置他们的适应度
                    if Fem < fit(x(i,:)) % 如果更新后的适应度小于原适应度
                        x(i,:) = X(:); % 则头狼位置也更新
                    end
                end
                if x(i,1) < X_min || x(i,1) > X_max % 判断是否超出边界
                    x(i,1) = X_min + rand(1) * (X_max - X_min);
                end
                if x(i,2) < X_min || x(i,2) > X_max
                    x(i,2) = X_min + rand(1) * (X_max - X_min);
                end
                fitt(i) = fit(x(i,:)); % 计算更新后的适应度
            end
            % 绘制狼群分布图
            scatter(x(:,1),x(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            scatter(X_First(1,1),X_First(1,2),'.','g')
            scatter(X_Second(1,1),X_Second(1,2),'.','g')
            scatter(X_Third(1,1),X_Third(1,2),'.','g')
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(Iter)])
            hold off
            pause(0.1)
            % 每一次迭代后重新选出三头狼王，用X_First、X_Second、X_Third表示；
            temp = [fitt x];
            temp_sort = sortrows(temp,1); % 对重组矩阵按第一列的进行升序
            X_First = temp_sort(1,2:end); % 选出狼王
            X_Second = temp_sort(2,2:end); % 选出次狼
            X_Third = temp_sort(3,2:end);  % 选出次次狼
            X_BEST = [X_BEST;X_First]; % 保存中间变量
            FIT_BEST = [FIT_BEST;temp_sort(1,1)]; % 保存中间变量
        end
        
        %% 结果可视化
        figure(2)
        plot(1:MAX_Iter,FIT_BEST)
        title('种群迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        [fit_B,ind_Best] = min(FIT_BEST);
        x_b = X_BEST(ind_Best,:);
        disp(['灰狼算法求解的最优解为：' ,num2str(fit_B)])
        disp(['灰狼算法求解的最优位置为：' ,num2str(x_b)])
        stop = 1;
    end
    drawnow
end

