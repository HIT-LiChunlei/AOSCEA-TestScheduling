classdef InstanceOP10 < PROBLEM
    % <multi> <real> <constrained>
    % case1 with 100 dies
    methods
        %% Initialization
        function Setting(obj)
            nDie         = 100;     % 晶片的数量
            nChainMax    = 100;     % 所允许的最多的测试链
            obj.M        = 2;
            obj.D        = nChainMax+nChainMax*nDie;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            load('Instance_10.mat');
            data.nDie      = 100;                                % 晶片的数量
            data.nChainMax = 100;                                % 所允许的最多的测试链
            data.P         = Instance_10.power;                  % 每个晶片对应的测试功率
            data.T         = Instance_10.time;                   % 每个晶片对应的测试时间
            data.Pmax      = 60000;                             % 每个测试组所允许的最大测试功率
            data.Ca        = 0.005;                             % 单位时间的测试成本
            data.Cb        = 993*4.24*1e-8;                     % 每个BIST测试结构的硬件成本
            data.Tmax      = sum(data.T)/2;                   % 所允许的最大测试时间

            pop = varargin{1};
            popsize = size(pop,1);
            Fit = zeros(popsize,obj.M);
            Cons = zeros(popsize,2);
            for np = 1:popsize
                %% 选择出激活的测试链
                x = pop(np,:);
                active = x(1:data.nChainMax); x(1:data.nChainMax) = []; % 取出前nChainMax数据，用作激活测试链指示
                x = reshape(x,data.nChainMax,data.nDie);                % 剩余数据按照nChainMax*nDie矩阵排列
                p = find(active>=0.5);                                  % 只有大于0.5的测试链才能被激活
                if isempty(p)
                    n_active  = 1;
                    active(1) = 1;
                else
                    n_active  = length(p);
                end
                F1 = data.Cb*n_active*1e6; % 测试链硬件成本，假设有one million(1e6)个chips需要测试

                %% 晶片分配测试链
                flag       = zeros(data.nDie,1);      % 指示晶片是否被分配到测试链
                ChainGroup = cell(data.nChainMax,1);  % 测试链上晶片分组情况记录
                ChainPsum  = zeros(data.nChainMax,1); % 每条测试链上的测试功率总和
                ChainTsum  = zeros(data.nChainMax,1); % 每条测试链上的测试时间总和
                ChainDnum  = zeros(data.nChainMax,1); % 每条测试链上的晶片数量
                recording  = [];                      % 分配过程数据记录
                while 1
                    p = find(flag==0);                % p表示还未分配测试链的晶片序号
                    if isempty(p); break; end         % p为空则表示已分配完毕，退出循环
                    matchMatric = zeros(length(p),data.nChainMax); % 初始化分配矩阵
                    PMatric     = matchMatric;        % 累计测试功率矩阵
                    TMatric     = matchMatric;        % 累计测试时间矩阵
                    for i = 1:length(p)
                        for noC = 1:data.nChainMax
                            noD = p(i);                           % 未分配测试链的晶片序号
                            if active(noC)<0.5
                                matchMatric(i,noC) = inf;         % 如果第noc条测试链未被激活，那么相应的匹配程度值设置为inf
                            else
                                matchMatric(i,noC) = x(noC,noD);  % 如果已经激活，第i个晶片对第noC条测试链的匹配程度设置为x中相应的值
                            end
                            PMatric(i,noC) = ChainPsum(noC);      % 当前测试链的累计测试功率
                            TMatric(i,noC) = ChainTsum(noC);      % 当前测试链的累计测试时间
                        end
                    end
                    temp1   = myMapminmax(matchMatric);           % 归一化匹配程度矩阵
                    temp2   = myMapminmax(PMatric);               % 归一化测试功率矩阵
                    temp3   = myMapminmax(TMatric);               % 归一化测试时间矩阵
                    temp    = temp1+temp2+temp3;                  % 最终的匹配程度由三方面体现
                    [p1,p2] = find(temp==min(min(temp)));         % p1表示temp中最小值的行序号,p2表示最小值的列序号
                    noD     = p(p1(1));                           % 选择出匹配程度最高的晶片与测试链
                    noC     = p2(1);
                    ChainPsum(noC)  = ChainPsum(noC)+data.P(noD); % 更新测试链上的测试功耗
                    ChainTsum(noC)  = ChainTsum(noC)+data.T(noD); % 更新测试链上的测试时间
                    ChainDnum(noC)  = ChainDnum(noC)+1;           % 更新测试链上晶片的个数
                    ChainGroup{noC} = [ChainGroup{noC},noD];      % 更新分组用的数据（此时未分组）
                    recording       = [recording;noC,noD,data.P(noD),data.T(noD)]; % 更新分配过程：测试链，晶片，功率，时间
                    flag(noD)       = 1;                          % 已分配的晶片标识置为1
                end
                %% 按顺序将测试链中的晶片分组
                dieGroup   = []; % 晶片分组情况
                recording1 = [];
                recording3 = [];
                noG        = 1;  % 测试组编号从1开始
                ST         = 0;  % 起始时间
                Pcons      = 0;  % 功率约束违反程度
                Tcons      = 0;  % 时间约束违反程度
                while 1
                    newG = [];
                    newT = zeros(data.nChainMax,1)+ST;
                    ST0  = ST;
                    for noC = 1:data.nChainMax
                        if ~isempty(ChainGroup{noC}) % 如果测试链中分配的有晶片
                            newG           = [newG,ChainGroup{noC}(1)];
                            ChainPsum(noC) = ChainPsum(noC)-data.P(ChainGroup{noC}(1)); % 更新测试链上的测试功耗
                            ChainTsum(noC) = ChainTsum(noC)-data.T(ChainGroup{noC}(1)); % 更新测试链上的测试时间

                            ST         = newT(noC);                            % 起始测试时间
                            ET         = newT(noC)+data.T(ChainGroup{noC}(1)); % 结束测试时间
                            newT(noC)  = newT(noC)+data.T(ChainGroup{noC}(1));
                            recording1 = [recording1;noC,noG,ChainGroup{noC}(1),data.P(ChainGroup{noC}(1)),data.T(ChainGroup{noC}(1)),ST,ET];
                            ChainGroup{noC}(1) = [];
                        end
                    end
                    % 测试链1 Die2 Die4
                    % 测试链3 Die3 Die1
                    % 经过上述操作之后，会将Die2和Die3放在一个测试组里面 newG
                    if isempty(newG); break; end

                    while 1
                        temp1 = myMapminmax(ChainPsum);
                        temp2 = myMapminmax(ChainTsum);
                        temp3 = myMapminmax(newT);
                        temp  = temp1+temp2-temp3;
                        [p1]  = find(temp==max(max(temp)));
                        noC   = p1(1); % 找出累计测试时间最短的测试链编号

                        if isempty(ChainGroup{noC}); break; end
                        if sum(data.P(newG))+data.P(ChainGroup{noC}(1))>data.Pmax; break; end %如果再往测试组里加一个功率超限，则退出循环

                        newG  = [newG,ChainGroup{noC}(1)];
                        ChainPsum(noC) = ChainPsum(noC)-data.P(ChainGroup{noC}(1));
                        ChainTsum(noC) = ChainTsum(noC)-data.T(ChainGroup{noC}(1));

                        ST         = newT(noC);
                        ET         = newT(noC)+data.T(ChainGroup{noC}(1));
                        recording1 = [recording1;noC,noG,ChainGroup{noC}(1),data.P(ChainGroup{noC}(1)),data.T(ChainGroup{noC}(1)),ST,ET]; % 1测试链编号 2测试组编号 3芯片编号 4功率 5时间 6起始时间 7截止时间
                        newT(noC)  = newT(noC)+data.T(ChainGroup{noC}(1));
                        ChainGroup{noC}(1) = [];
                    end
                    recording3 = [recording3;noG,ST0,max(newT),sum(data.P(newG))];
                    if sum(data.P(newG))>data.Pmax
                        Pcons = Pcons+sum(data.P(newG))-data.Pmax;
                    end
                    ST  = max(newT);
                    noG = noG+1;
                    dieGroup = [dieGroup,{newG}];
                end
                TestTime = max(recording1(:,7))*1/50/1e6;
                F2 = data.Ca*TestTime*1e6;

                if data.Tmax<max(recording1(:,7))
                    Tcons = Tcons+(max(recording1(:,7))-data.Tmax);
                end
                Fit(np,:) = [F1,F2];
                Cons(np,:) = [Pcons,Tcons];
            end
            Population = SOLUTION(varargin{1},Fit,Cons,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [4.2103e+03 19209.9]*10; % HV
            % R = [589.4448,5431.1173]; % IGD
        end
    end
end