% parameters:
a = 20; b = 21; h = b-a;
rho=7932; mu=rho*3260^2; lbd=rho*5960^2-2*mu;	
mat = Material('noname', lbd, mu, rho);
N = 12;
k = linspace(1e-3, 0.5, 100)/h; % wavenumber-thickness (solve for frequency)
cyl = Cylinder(mat, [a, b], N);
load('data/disperse.mat')

%% plot
gew = cyl.fullyCoupled(0); % first order circumferential
dat0 = computeW(gew, k); 
gew = cyl.fullyCoupled(1); % second order circumferential
dat1 = computeW(gew, k); 

% plot against reference data:
fig = figure; hold on
xlim([0, 0.5]), ylim([0, 200])
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
plot(disperse.k/1e-3, disperse.f*1e6, 'g.')
ax = gca; ax.ColorOrderIndex = 1; % reset color order 
plot(dat0.k(:), dat0.w(:)/2/pi, '.'); 
plot(dat1.k(:), dat1.w(:)/2/pi, '.'); 
legend({'disperse', 'n = 0', 'n = 1'})

% request user to evaluate test:
assert( userTestConfirmation(fig) )

% some unused stuff:
% pass = uiconfirm(fig, 'pass?','Is plot good?')
% uiwait(fig)
% function [fig, ax] = figureConfirm()
% fig = uifigure('WindowButtonDownFcn',@(src,event)uiresume(src));
% ax = uiaxes('Parent',fig,...
%             'Units','pixels',...
%             'Position', [25, 70, 500, 340]);   
% passBtn = uibutton(fig,'push',...
%                'Text', 'pass',...
%                'Position',[200, 10, 100, 40],...
%                'UserData',false,...
%                'ButtonPushedFcn', @(btn, event) passButtonPushed(btn, fig));
% failBtn = uibutton(fig,'push',...
%                'Text', 'fail',...
%                'Position',[350, 10, 100, 40],...
%                'UserData',false,...
%                'ButtonPushedFcn', @(btn, event) passButtonPushed(btn, fig));
% end
% function passButtonPushed(btn, fig)  
%         btn.UserData = true; 
%         disp(btn.Text);
%         close(fig);
% end
% function choice = waitForButtonPress(fig)
%     uiwait(fig); % wait for button press
%     btns = findobj(fig, 'Type', 'uibutton');
% end
