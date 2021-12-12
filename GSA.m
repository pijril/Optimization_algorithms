%% 万有引力算法实现
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
        N = 20; % 种群的大小
        D = 2; % 种群的维数
        alph = 5; % G中常数取值 经过我多次手动测试，最好的alph值为8
        sigo = 00.1; % 防止出现分母为0的数
        Iter_max = 100; % 迭代总次数
        Iter = 0; % 初始迭代次数
        X_MAX = 100; % 粒子的上界
        X_MIN= -100; % 粒子的下界
        m = zeros(N,1); % 初始化粒子的惯性质量
        M = zeros(N,1); % 初始化粒子的质量
        A = zeros(N,D); % 初始化粒子的加速度
        V = zeros(N,D); % 初始化粒子的速度
        F = zeros(N,D); % 初始化粒子的万有引力
        f = zeros(N,D); % 初始化中间变量
        G_0 = 100; % 定义万有引力常量初始值
        fitt = zeros(N,1); % 初始化种群适应度
        X = zeros(N,D);
        % 初始化种群及适应度
        x = X_MIN + (X_MAX-X_MIN) .* rand(N,D);
        for i =1:N
            fitt(i) = fit(x(i,:)); % 计算第i个物体的适应度
        end
        X_BEST = [];
        FITT = [];
        
        %% 迭代开始
        while Iter < Iter_max
            Iter = Iter + 1; % 每运行一次，迭代次数加一
            % 惯性质量的计算
            worst = max(fitt); % 找到种群中最差的适应度
            best = min(fitt); % 找到种群中最优的适应度
            for i =1:N
                m(i) =  (fitt(i) - worst) / (best - worst); % 计算惯性质量
            end
            for i = 1:N
                M(i) = m(i) / sum(m); % 计算质量
            end
                G = G_0 * exp(-( alph * Iter) / Iter_max); % 计算第Iter时的G的值
%             G = G_0 * (1-Iter/Iter_max); % 这里是对G的一种改进，将非线性改为线性下降，但其结果并不理想。
            % 下面计算物体万有引力的大小
            for i =1:N % 遍历所有物体
                for j = 1:N % 当第i个物体时，需要计算第i个物体与所有N个物体之间的引力大小
                    R = sqrt((x(i,1) - x(j,1))^2 + (x(i,2) - x(j,2))^2);
                    for q = 1:D
                        % 计算在Iter时刻，物体i在D维上收到第j个粒子的引力大小
                        f(j,q) = G * ((M(i) * M(j)) / (R + sigo)) * (x(j,q) - x(i,q));
                    end
                end
                for q = 1:D
                    F(i,q) = rand(1,N) * f(:,q); % 计算第i个粒子的第q维的总引力
                end
            end
            % 计算物体加速度
            for i =1:N
                if M(i) ~= 0
                    for q = 1:D
                        A(i,q) = F(i,q) / M(i); % 计算第i个物体第q维的加速度
                    end
                else
                    A(i,:) = 0; % 如果第i个物体的惯性质量为0，则说明其位置差，则给其加速度为0
                end
            end
            % 计算物体速度
            for i = 1:N
                for q = 1:D
                    V(i,q) = rand() * V(i,q) + A(i,q); % 计算第i个粒子第q维的粒子速度
                end
            end
            % 更新位置
            for i = 1:N
                for q = 1:D
                    X(i,q) = x(i,q) + V(i,q);% 由速度更新位置
                    if X(i,q) > X_MAX
                        X(i,q) = X_MIN + (X_MAX-X_MIN) * rand();
                        %                 X(i,q) = X_MAX;
                    end
                    if X(i,q) < X_MIN
                        X(i,q) = X_MIN + (X_MAX-X_MIN) * rand();
                        %                 X(i,q) = X_MIN;
                    end
                end
                if fitt(i) > fit(X(i,:))
                    x(i,:) = X(i,:);
                    fitt(i) = fit(X(i,:));
                end
            end
            
            
            % 绘制一次迭代后的粒子分布情况
            figure(1)
            scatter(x(:,1),x(:,2),'.','r'); % 绘制散点图
            hold on
            scatter(x_a,x_b,'.','b')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(Iter)])
            pause(0.02)
            % 记录中间变量
            [fit_min,index_best] = min(fitt); % 找到最优位置
            X_BEST = [X_BEST;x(index_best,:)];
            FITT = [FITT;fit_min];
        end
        [fit_best,index_BEST] = min(FITT); % 找到全局最优的解
        disp(['最优的结果为：' num2str(fit_best)])
        disp(['最优的结果为：' num2str(X_BEST(index_BEST,:))])
        figure(2)
        plot(1:Iter_max,FITT)
        title('万有引力算法迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        stop =1;
    end
    drawnow
end


