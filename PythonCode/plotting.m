format long

if (Model == 1)
    B = readmatrix('execution_timeM1.txt');
    python_times1 = B
    python_times1 = python_times1(2:101)./1000000
    
    
    f = figure("Name", "myfig");
    plot(python_times1);
    hold on;
    plot(tosave1);
    hold off;
    legend('Python','Matlab')
    title('Execution time');
    ylabel('time [seconds]')
    xlabel('step')
    %ylim([0, 0.004])
    grid
    saveas(f,"python_vs_matlabM1.png")
end

if (Model == 2)
    B = readmatrix('execution_timeM2.txt');
    python_times2 = B
    python_times2 = python_times2(2:101)./1000000
    
    
    f = figure("Name", "myfig");
    plot(python_times2);
    hold on;
    plot(tosave2);
    hold off;
    legend('Python','Matlab')
    title('Execution time');
    ylabel('time [seconds]')
    xlabel('step')
    %ylim([0, 0.004])
    grid
    saveas(f,"python_vs_matlabM2.png")
end


%% Plots
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
