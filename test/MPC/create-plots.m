%%  Create plots of MPC examples.
% 5/7/11: This code may not work with the current setup.

disp(['--  There are ',num2str(n_act),' figures.  --']);

outputNames = {'sen[1]',...
    'sen[2]'};
disturbanceNames = {'addDisturbance[1].u2',...
    'addDisturbance[1].u2'};
noiseNames = {'addNoise[1].u2',...
    'addNoise[1].u2'};

Dym = dymload([SimulinkModel,'_Dym.mat']); % Dymola result data
t_Dym = dymget(Dym,'Time'); % time vector for Dymola result data

% Plot the response of each control loop
for loop=1:n_act;
    figure(loop); % Create or activate the figure
    clf(loop,'reset'); % Clear the figure
    set(loop,'name',['LQMPC - Loop ',num2str(loop)]) % Set the title of the figure

    subplot(2,2,1);
    hold on;
    plot(Sim.Controller.ref.Time,Sim.Controller.ref.Data(:,loop),'g--','DisplayName','Set-point');
    plot(t_Dym,dymget(Dym,outputNames{loop}),'b-','DisplayName','Measurement');
    %    plot(Sim.Controller.obs.Time,Sim.Controller.obs.Data(:,loop),'r-','DisplayName','Observation');

    axis tight
    ax = axis;
    y_margin = (ax(4) - ax(3))/20;
    axis([0, runtime, ax(3) - y_margin, ax(4) + y_margin]);
    title('Set-point and measurement');
    legend('Location','East');

    subplot(2,2,2);
    hold on;
    plot(Sim.Controller.act_fb.Time,-Sim.Controller.act_fb.Data(:,loop),'r-','DisplayName','-Feedback');
    plot(Sim.Controller.act_ff.Time,Sim.Controller.act_ff.Data(:,loop),'b--','DisplayName','Feedforward');
    plot(Sim.Controller.act_opt.Time,Sim.Controller.act_opt.Data(:,loop),'g:','DisplayName','Perturbation by optimizer');
    axis tight
    ax = axis;
    y_margin = (ax(4) - ax(3))/20;
    axis([0, runtime, ax(3) - y_margin, ax(4) + y_margin]);
    title('Composition of actuation');
    legend('Location','East');

    subplot(2,2,3);
    hold on;
    plot(Sim.Controller.act.Time,Sim.Controller.act.Data(:,loop),'r-','DisplayName','Actuator input');
    axis tight
    ax = axis;
    y_margin = (ax(4) - ax(3))/20;
    axis([0, runtime, ax(3) - y_margin, ax(4) + y_margin]);
    plot([0, runtime],act_min(loop)*ones(2,1),'k--','DisplayName','Lower limit');
    plot([0, runtime],act_max(loop)*ones(2,1),'k:','DisplayName','Upper limit');
    title('Total actuation');
    legend('Location','East');
    xlabel('Time /s');

    subplot(2,2,4);
    hold on;
    plot(t_Dym,dymget(Dym,disturbanceNames{loop}),'r-','DisplayName','Disturbance');
    plot(t_Dym,dymget(Dym,noiseNames{loop}),'b--','DisplayName','Noise');
    axis tight
    ax = axis;
    y_margin = (ax(4) - ax(3))/20;
    axis([0, runtime, ax(3) - y_margin, ax(4) + y_margin]);
    title('Actuator disturbance and sensor noise');
    legend('Location','East');
    xlabel('Time /s');
end

for state=1:n_x
    figNum = n_act + state;
    figure(figNum); % Create or activate the figure
    clf(figNum,'reset'); % Clear the figure
    set(figNum,'name',['LQMPC - State ',num2str(state)]) % Set the title of the figure
    DymInd = ['[',num2str(state),']'];
    hold on;
    plot(t_Dym,dymget(Dym,['plant.x',DymInd]),'r-','DisplayName','Actual state');
    plot(Sim.Controller.x_est.Time,Sim.Controller.x_est.Data(:,state),'b--','DisplayName','Estimated state');
    axis tight
    ax = axis;
    y_margin = (ax(4) - ax(3))/20;
    axis([0, runtime, ax(3) - y_margin, ax(4) + y_margin]);
    title('Actual and Estimated State');
    legend('Location','East');
    xlabel('Time /s');
end

