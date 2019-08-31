    
****ACKNOWLEDGMENTS****
   A large chunk of the code used (~50%) is written by Professor Alexey Rubtsov and the original code can be found here https://sites.google.com/site/alexeyvrubtsov/quantitative-finance-materials within the "Portfolio Optimization" zip file. I could not have finished this project without his guidance.




%% Attempt to read data from stockData.xlsx and create dynamic backtesting
%script that allows you to change the start date and the length of the
%backtesting. The script will record portfolio value and rebalance
%portfolio at the beginning of each new month

%reading data from source: yahoo finance),
%all prices in $USD

clear
clc

format long

Short = 0; % if Short=1, we allow for short positions; if Short = 0, we don't.
numScenarios = 15;  %number of scenarios generated when creating Rs linspace... determines the # of portfolios along the efficient frontier

%% Reading in and cleaning the data

%import data
[num,text,raw] = xlsread('energyStocks.xlsx');

%num is the array of closing prices will be used to calculate returns matrices and portfolio values in
%backtesting process

%calculate log returns since data is daily close price for each stock
returns = diff(log(num));

%remove NaN rows, 87 found in the original data set of 4567 elements 
NaNarray = isnan(returns); %create logical array where 1 denotes a NaN and 0 for anything else
for i = length(returns):-1:1 %delete from end of array to avoid dealing with updated indices while deleting rows in data set
    if NaNarray(i,1) == 1
        returns(i,:) = [];
        text(i,:) = [];
        raw(i,:) = [];
        num(i,:) = [];        
    end
end

%convert the date character vector into serial date number
t = datenum(text(:,1));

startDate = '01/04/2013'; %the start date for the backtesting script, assumed to be at the beginnning of a month
endDate = '06/01/2018'; %end date

startDateSerial = datenum(startDate);   %enter date to start backtesting, datenum converts date to serial number
startIndex = find(t == startDateSerial); %finds index in which the startDate is located in serial date arary
initialMonthIndex = startIndex; %record first day of first month index for purpose of calculating portfolio returns

endDateSerial = datenum(endDate);  %end date for backtesting,
endIndex = find(t == endDateSerial);  %find index in which date is located in t array (array of dates in serial form)

%calculate numMonths which will be the number of months between the
%start and end date

numYears = year(endDate) - year(startDate);
numDays = abs(datenum(startDate) - datenum(endDate));
numMonths = numYears*12 + month(endDate) - month(startDate);

%load all hisorical data into array 
rowLength = size(returns);
currentData = returns(1:startIndex,1:rowLength(2)); %current data, this is the historical data to calculate for the 5 stocks

%% Calculate efficient frontier using historical data    
% Expected returns and covariance matrix
[T,N] = size(currentData); % N is the number of equities
Sigma = cov(currentData)*22; %Scaling daily covariance and mean by number of days in a trading month = 22
mu = mean(currentData)*22;
% Plotting the distribution of the equally weighted portfolio
wn = ones(N,1)/N; % equal weighted portfolio weights
figure(1)
alpha1 = 0.01;
CVaR1 = cvarhist(wn,currentData,alpha1);
alpha2 = 0.05;
CVaR2 = cvarhist(wn,currentData,alpha2);
h1 = hist(currentData*wn,T/10);
hist(currentData*wn,T/10)
set(gca,'Box', 'On', 'LineWidth',1.5, 'FontSize', 12)
hold on
plot([CVaR1 CVaR1],[0 max(h1)],'r', 'LineWidth', 1.5)
plot([CVaR2 CVaR2],[0 max(h1)],'r--', 'LineWidth', 1.5)
legend('returns',['CVar-Level, \alpha = ' num2str(alpha1)],['CVar-Level, \alpha = ' num2str(alpha2)])

%% Optimization for calculating effiicient frontier at start date of backtesting
Rs = linspace(max(min(mu),0), max(mu),numScenarios); % target portfolio returns (no negative returns)
alpha = 0.1; % confidence level for CVaR
% Set the lower and upper bounds on portfolio weights. 
% To enforce diversification, we could also use, e.g., UB = 1/N*min(N,10);
UB = 1;
LB = -UB;

% optimization loop
% If CVaR is defined as a negative value (i.e., as a gain, not as a loss),
% we want to maximize CVaR or minimize -CVaR.
% Since MATLAB does minimization we minimize -CVaR
objfun = @(w) -(cvarhist(w,currentData,alpha));
% Markowitz problem for comparison
meanvar = @(w) w'*Sigma*w;
% Initial guess
w0 = wn; % equally weighted portfolio
ES0 = cvarhist(w0,currentData,alpha);
for i = 1:length(Rs) % for each target expected return
    R0 = Rs(i);    % target expected return i
    A = -mu; % negative of expected returns
    if Short==0
        % Constraints: 1) expected returns should not be less than target
        % returns, 2) weights should be nonnegative, 3) weights should be
        % less than upper bound UB
        A = [A;-eye(N)];
        A = [A; eye(N)];
        b = [-R0 zeros(1,N) UB*ones(1,N)];
    elseif Short==1
        A=[A;-eye(N)];
        A=[A; eye(N)]; 
        b=[-R0 -LB*ones(1,N) UB*ones(1,N)];
    elseif Short~=0||Short~=1
        error('Set Short=1 or 0')
    end
    b = b';
    Aeq = ones(1,N);
    beq = 1;
    options = optimset(optimset('fmincon'), 'Display', 'off', 'Algorithm', 'sqp');
    % CVaR optimization
    [w,~,~,~,~,~,hessian] = fmincon(objfun,w0,A,b,Aeq,beq,[],[],[],options);
    % Mean-Variance optimization
    [wm,~,~,~] = fmincon(meanvar,w0,A,b,Aeq,beq,[],[],[],options);
    % Evaluate data for plots
    CVaR(i) = -cvarhist(w,currentData,alpha);
    CVaRwm(i) = -cvarhist(wm,currentData,alpha);
    SD(i) = std(currentData*w);
    R(i) = mean(currentData*w);
    effSig(i) = std(currentData*wm);
    effMu(i) = mean(currentData*wm);
    Weight(:,i) = w;
    MarkoW(:,i) = wm;
end
weightHistoryMinVar = [Weight(:,1)]; %record weights of minimum variance portfolio 
weightHistoryMinVarMarkowitz = [MarkoW(:,1)]; %record weight of markowitz frontier minimum variance portfolioIWRZ

weightHistoryMarkowitzP10 = [MarkoW(:,10)];


%del = weightHistoryMinVar(:,1)./num(startIndex,:);
%% begin backtesting forloop walking forword 1 month at a time until end date
currentIndex = startIndex;
updatedData = currentData;
for i = 1:numMonths %put i = 1:numMonths in end, sets the number of times portfolio will rebalance
    index = currentIndex;
    
    %for-loop to find the index of the first day of the new month since
    %not all months have a uniform amount of days
    for j = index:index+29
        if month(t(j+1)) ~= month(t(j))
            newMonthStart = j+1; 
            initialMonthIndex = [initialMonthIndex newMonthStart];
            break;
        else
           continue; 
        end
    end
    %add extra month of data to updatedData set to be able to recalculate
    %new covariance and expected returns matrixd
    updatedData = [updatedData; returns(currentIndex+1:newMonthStart,:)];
    currentIndex = newMonthStart;
    [T,N] = size(updatedData); % N is the number of equities
    Sigma = cov(updatedData)*22; %Scaling daily covariance and mean by number of days in a trading month = 22
    mu = mean(updatedData)*22;
    wn = ones(N,1)/N; % equal weighted portfolio weights
    alpha1 = 0.01;
    CVaR1 = cvarhist(wn,updatedData,alpha1);
    alpha2 = 0.05;
    CVaR2 = cvarhist(wn,updatedData,alpha2);

    %% Optimization for calculating effiicient frontier at start date of backtesting
    Rs = linspace(max(min(mu),0), max(mu),numScenarios); % target portfolio returns (no negative returns)
    alpha = 0.1; % confidence level for CVaR
    % Set the lower and upper bounds on portfolio weights. 
    % To enforce diversification, we could also use, e.g., UB = 1/N*min(N,10);
    UB = 1;
    LB = -UB;

    % optimization loop
    % If CVaR is defined as a negative value (i.e., as a gain, not as a loss),
    % we want to maximize CVaR or minimize -CVaR.
    % Since MATLAB does minimization we minimize -CVaR
    objfun = @(w) -(cvarhist(w,updatedData,alpha));
    % Markowitz problem for comparison
    meanvar = @(w) w'*Sigma*w;
    % Initial guess
    w0 = wn; % equally weighted portfolio
    ES0 = cvarhist(w0,updatedData,alpha);
    for i = 1:length(Rs) % for each target expected return
        R0 = Rs(i);    % target expected return i
        A = -mu; % negative of expected returns
        if Short==0
            % Constraints: 1) expected returns should not be less than target
            % returns, 2) weights should be nonnegative, 3) weights should be
            % less than upper bound UB
         A = [A;-eye(N)];
            A = [A; eye(N)];
            b = [-R0 zeros(1,N) UB*ones(1,N)];
        elseif Short==1
            A=[A;-eye(N)];
            A=[A; eye(N)]; 
            b=[-R0 -LB*ones(1,N) UB*ones(1,N)];
        elseif Short~=0||Short~=1
            error('Set Short=1 or 0')
        end
        b = b';
        Aeq = ones(1,N);
        beq = 1;
        options = optimset(optimset('fmincon'), 'Display', 'off', 'Algorithm', 'sqp');
        % CVaR optimization
        [w,~,~,~,~,~,hessian] = fmincon(objfun,w0,A,b,Aeq,beq,[],[],[],options);
        % Mean-Variance optimization
        [wm,~,~,~] = fmincon(meanvar,w0,A,b,Aeq,beq,[],[],[],options);
        % Evaluate data for plots
        CVaR(i) = -cvarhist(w,updatedData,alpha);
        CVaRwm(i) = -cvarhist(wm,updatedData,alpha);
        SD(i) = std(updatedData*w);
        R(i) = mean(updatedData*w);
        effSig(i) = std(updatedData*wm);
        effMu(i) = mean(updatedData*wm);
        Weight(:,i) = w;
        MarkoW(:,i) = wm;        
    end
    weightHistoryMinVar = [weightHistoryMinVar Weight(:,1)]; %update weight history with each iteration to record weight history
    weightHistoryMinVarMarkowitz = [weightHistoryMinVarMarkowitz MarkoW(:,1)];
    
    weightHistoryMarkowitzP10 = [weightHistoryMarkowitzP10 MarkoW(:,10)];
end

%% Plot the results
% Plot of the CVaR efficient frontier
figure(2) 
hold on
plot(CVaR,R*100,'k','LineWidth',1.5) % CVaR when CVaR is optimal
plot(CVaRwm,effMu*100,'k--','LineWidth',1.5) % CVaR for Markowitz
set(gca,'Box','on','LineWidth', 1.5, 'FontSize',12')
grid on
ylabel('Portfolio return (in %)')
title('CVaR efficient frontier')
legend('CVaR','Markowitz')
legend('Location','NorthWest')
xlabel('CVaR')

% % Plot of the Mean-Variance efficient frontier
% figure(3) 
% hold on
% plot(SD,R*100,'k','LineWidth',1.5) % CVaR optimal
% plot(effSig,effMu*100,'k--','LineWidth',1.5) % Markowitz optimal
% set(gca,'Box','on','LineWidth', 1.5, 'FontSize',12')
% grid on
% ylabel('Portfolio return (in %)')
% legend('CVaR','Markowitz')
% legend('Location','NorthWest')
% title('Mean-variance efficient frontier')
% xlabel('Standard deviation')

% Plot of the portfolio weights
figure(3)
subplot(1,2,1)
area(Weight');
title('CVaR');
set(get(gcf,'Children'),'YLim',[0 1]);
xlabel('Portfolio')
ylabel('Asset weight')

subplot(1,2,2)
area(MarkoW');
title('Markowitz');
set(get(gcf,'Children'),'YLim',[0 1]);
xlabel('Portfolio')
ylabel('Asset weight')

% saveas(2,'FigureCVaR.png')
% saveas(3,'FigureMV.png')
% saveas(4,'FigurePortfolios.png')

%% calculate portfolio value a.k.a portfolio performance by theoretically
%investing 1 dollar and given the weights, see how much of the fraction of
%each share we buy, then multiply that fraction by the share price one
%month in advance to determine wether the portfolio gained or lost value

portfolioValueMarkowitz = []; %min varian markowitz optimal portfolio
portfolioValueCVAR = []; %min variance cvar optimal portfolio
portfolioValueMarkowitzP10 = []; %for portfolio #10
for i = 1:numMonths  %numMonths 
    del = weightHistoryMinVar(:,i)'./num(initialMonthIndex(i),:);  %invest 1 dollar, find the fraction of each share we own
    val = del.*num(initialMonthIndex(i+1),:); %multiply fraction of share we own with new share prices, sum up vector to find total return
    portfolioValueCVAR(i) = sum(val); %sum vector to find new portfolio value (compare to value of $1 at start)
    
    del2 = weightHistoryMinVarMarkowitz(:,i)'./num(initialMonthIndex(i),:);
    val2 = del2.*num(initialMonthIndex(i+1),:);
    portfolioValueMarkowitz(i) = sum(val2);
    
    del3 = weightHistoryMarkowitzP10(:,i)'./num(initialMonthIndex(i),:);
    val3 = del3.*num(initialMonthIndex(i+1),:);
    portfolioValueMarkowitzP10(i) = sum(val3);
end

portfolioReturnsCVAR = portfolioValueCVAR - 1; %subtract 1 from total return to obtain rate of return
portfolioReturnsMarkowitz = portfolioValueMarkowitz - 1;%do the same for the markowitz efficient frontier portfolio
portfolioReturnsMarkowitzP10 = portfolioValueMarkowitzP10 - 1;


riskless = zeros(length(portfolioReturnsMarkowitz));

SharpeMarkowitz = sharpe(portfolioReturnsMarkowitz,0)
SharpeCVAR = sharpe(portfolioReturnsCVAR,0)
SharpeMarkowitzP10 = sharpe(portfolioReturnsMarkowitzP10,0);

figure(5)
plot(portfolioReturnsMarkowitz)
title('Min Var Portfolio returns')
xlabel('Month')
ylabel('Rate of Return')

meanCVAR = mean(portfolioReturnsCVAR);
stCVAR =  std(portfolioReturnsCVAR);


    
