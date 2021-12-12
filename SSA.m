% 松鼠算法实现
clc;clear
tic
global x_a x_b
x_a = 0;x_b=0;
SS_t =[];
figure(1)
plotbutton = uicontrol('style','pushbutton','string','运行','fontsize',12,'position',[10,400,50,20],'callback','run = 1;');
stop = 0;run = 0; % 初始化界面标志
while stop == 0
    if run == 1
        % 参数初始化
        t_max = 100; % 最大迭代次数
        t = 1; % 初始迭代次数
        n = 50; % 松鼠的个数
        d = 2; % 带搜索空间维数
        P_dp = 0.1; % 捕食者出现概率
        x_min = -100; % 自变量下界
        x_max = 100; % 自变量上界
        G_c = 1.9; % 滑动常数
        beta = 1.5; % levy中的常数
        times = 0;
        
        % 种群初始化
        X = x_min + rand(n,d).* (x_max - x_min);
        fitt = zeros(1,n); % 初始化内存空间
        % 计算适应度
        for i = 1:n
            fitt(i) = fit(X(i,:));
        end
        % 下面进行种群分类
        fitt_X = [fitt' X];
        fitt_X_sort = sortrows(fitt_X,1); % 先排序，方便后面分类
        FS_ht = fitt_X_sort(1,2:3); % 适应度最小的处在核桃树上
        FS_at = fitt_X_sort(2:4,2:3); % 适应度较小的处在橡树上
        FS_nt = fitt_X_sort(5:end,2:3); % 适应度一般的处在普通树上
        
        
        TT = [];
        XX = [];
        while t<t_max
            t = t + 1;
            % 更新处在橡树的松鼠——》核桃树上
            for i = 1:length(FS_at)
                % 进行位置更新
                if rand > P_dp
                    d_g = cl_dg();
                    FS_at(i,:) = FS_at(i,:) + d_g * G_c * (FS_ht - FS_at(i,:));
                else
                    FS_at(i,:) = x_min + rand(1,d).* (x_max - x_min);
                end
            end
            % 更新普通树上的松鼠——》橡树上
            for i = 1:length(FS_nt)
                % 进行位置更新
                d_g = cl_dg();
                FS_at_c = mean(FS_at);
                if rand > P_dp
                    FS_nt(i,:) = FS_nt(i,:) + d_g * G_c * (FS_at_c - FS_nt(i,:));
                else
                    FS_nt(i,:) = x_min + rand(1,d).* (x_max - x_min);
                end
            end
            % 更新在普通树上的松鼠——》核桃树上
            for i =1:length(FS_nt)
                % 进行位置更新
                if rand > P_dp
                    d_g = cl_dg();
                    FS_nt(i,:) = FS_nt(i,:) + d_g * G_c * (FS_ht - FS_nt(i,:));
                else
                    FS_nt(i,:) = x_min + rand(1,d).* (x_max - x_min);
                end
            end
            % 计算当前季节最小常数值S_min
            S_min = 10e-2/((365)^(t/(t_max/2.5)));
            % 计算季节常量
            temp = FS_ht - FS_at;
            s = 0;
            for i =1:3
                s = s + temp(i,1)^2 + temp(i,2)^2;
            end
            S_t = sqrt(s);
            SS_t = [SS_t;S_t];
            % 判断季节
            
            if S_t < S_min
                times = times+1;
                sigma_u  = (gamma(1+beta)*sin(pi*beta*1/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
                % 遍历每一个松鼠，然后levy飞行
                for i = 1:length(FS_nt)
                    for j = 1:d
                        levy = 0.1*(normrnd(0,1)*sigma_u)/(abs(normrnd(0,1))^(1/beta));
                        FS_nt(i,j) = x_min + levy * (x_max - x_min);
                    end
                end
            end
            % 绘制迭代过程图
            % 绘制普通松鼠（蓝色点）
            scatter(FS_nt(:,1),FS_nt(:,2),'.','b'); % 绘制散点图
            hold on
            % 绘制橡胶树的松鼠（绿色）
            scatter(FS_at(:,1),FS_at(:,2),'.','g')
            % 绘制核桃树的松鼠（紫色）
            scatter(FS_ht(:,1),FS_ht(:,2),'*','m')
            % 绘制最优点位置
            scatter(x_a,x_b,'.','r')
            hold off
            axis([-100,100,-100,100])
            title(['当前迭代次数为：',num2str(t)])
            pause(0.07)
            % 更新新一代的种群情况
            X = [FS_nt;FS_ht;FS_at];
            for i = 1:n
                fitt(i) = fit(X(i,:));
            end
            TT = [TT;min(fitt)]; % 记录核桃树上的松鼠位置
            ind = find(fitt == min(fitt),1);
            X_best = X(ind,:); % 更新最优位置
            XX =[XX;X_best];
            % 下面进行种群分类
            fitt_X = [fitt' X];
            fitt_X_sort = sortrows(fitt_X,1); % 先排序，方便后面分类
            FS_ht = fitt_X_sort(1,2:3); % 适应度最小的处在核桃树上
            FS_at = fitt_X_sort(2:4,2:3); % 适应度较小的处在橡树上
            FS_nt = fitt_X_sort(5:end,2:3); % 适应度一般的处在普通树上
            
        end
        % 结果输出
        Fit_Best = min(TT);
        ind_best = find(TT == min(TT),1);
        X_BE = XX(ind_best,:);
        disp(['找到的最优解为：',num2str(Fit_Best)])
        disp(['最优解对应的X：',num2str(X_BE)])
        disp(['季节更新次数：',num2str(times)])
        figure(2)
        plot(1:length(TT),TT)
        title('种群迭代图')
        xlabel('迭代次数')
        ylabel('适应度值')
        stop = 1;
    end
    drawnow
end