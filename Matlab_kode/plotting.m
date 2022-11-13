format long

%Reading the execution time from Python Model 1
B = readmatrix('execution_timeM1.txt');
python_times1 = B
python_times1 = python_times1(2:101).*10^(-9)

%Reading while iterations from Python Model 1
W1 = readmatrix('while_iterations1.txt');
python_iterations1 = W1;
python_iterations1 = W1(2:101);

%Reading the execution time from Python Model 2
C = readmatrix('execution_timeM2.txt');
python_times2 = C
python_times2 = python_times2(2:101).*10^(-9)

%Reading while iterations from Python Model 2
W2 = readmatrix('while_iterations2.txt');
python_iterations2 = W2;
python_iterations2 = W2(2:101);

if (Model == 1)
    
    f = figure("Name", "myfig");
    axes('XScale', 'linear', 'YScale', 'log')
    hold on;
    plot(python_times1, 'LineWidth', 2);
    hold on;
    plot(tosave1, 'LineWidth', 2);
    hold off;
    legend('Python','Matlab')
    title('Execution time model 1');
    %set(gca, 'YScale', 'log')
    ylabel('time [seconds]')
    xlabel('step')
    %ylim([0, 0.00004])
    grid
    saveas(f,"python_vs_matlabM1.png")
end

if (Model == 2)

    f = figure("Name", "myfig");
    axes('XScale', 'linear', 'YScale', 'log')
    hold on;
    plot(python_times2);
    hold on;
    plot(tosave2);
    hold off;
    legend('Python','Matlab')
    title('Execution time model 2');
    ylabel('time [seconds]')
    xlabel('step')
    %ylim([0, 0.004])
    grid
    saveas(f,"python_vs_matlabM2.png")
end


%% Plots
%Plot matlab runtimes vs python runtimes against each other for each model
fig = figure;
subplot(2,1,1)
plot(python_times1);
hold on;
plot(tosave1);
hold off;
legend('Python','Matlab')
title('Execution time model 1');
ylabel('time [seconds]')
xlabel('step')
%ylim([0, 0.004])
grid

subplot(2,1,2)
plot(python_times2);
hold on;
plot(tosave2);
hold off;
legend('Python','Matlab')
title('Execution time model 2');
ylabel('time [seconds]')
xlabel('step')
%ylim([0, 0.004])
grid

%Plot matlab while iterations vs python while iterations against each other
%for each model

fig = figure;
subplot(2,1,1)
plot(python_iterations1, 'LineWidth', 2, 'LineStyle','--');
hold on;
plot(while_iterations1, 'LineWidth', 1);
hold off;
legend('Python','Matlab')
title('While iterations each step model 1');
ylabel('iterations')
xlabel('step')

ylim([0,7])
grid

subplot(2,1,2)
plot(python_iterations2, 'LineWidth', 2, 'LineStyle','--');
hold on;
plot(while_iterations2, 'LineWidth', 1);
hold off;
legend('Python','Matlab')
title('While iterations each step model 2');
ylabel('iterations')
xlabel('step')
ylim([0,5])
grid
