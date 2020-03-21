function [h] = plotStabDiag(fn,Az,fs,stablity_status,Nmin,Nmax)
% -------------------------------------------------------------------------
% [h] = plotStabDiag(fn,Az,fs,stablity_status,Nmin,Nmax) plots the
% stabilization diagram of the identified eigen frequencies as a function
% of the model order, calculated with the SSI-COV method.
% -------------------------------------------------------------------------
% Input:
% fn: cell : eigen frequencies identified for multiple system orders.
% Az : vector: Time serie of acceleration response (illustrative purpose)
% fs: sampling frequency
% stablity_status: cell of stability status for each model order
% Nmin: scalar: minimal number of model order
% Nmax: scalar: maximal number of model order
% Output: h: handle of the figure
% -------------------------------------------------------------------------
% See also: SSICOV.m
% -------------------------------------------------------------------------
% Author: Etienne Cheynet, UIS
% Updated on: 08/03/2016
% -------------------------------------------------------------------------
Npoles =Nmin:1:Nmax;
[Saz,f]=pwelch(Az,[],[],[],fs);


h = figure;
ax1 = axes;
hold on;box on
for jj=0:4,
    y = [];
    x = [];
    for ii=1:numel(fn)
        ind = find(stablity_status{ii}==jj);
        x = [x;fn{ii}(ind)'];
        y = [y;ones(numel(ind),1).*Npoles(ii)];
    end
    x1{jj+1}=x;
    y1{jj+1}=y;
end

h1=plot(x1{1},y1{1},'k+','markersize',5);% new pole
h2=plot(x1{2},y1{2},'ko','markerfacecolor','r','markersize',5);  % stable pole
h3=plot(x1{3},y1{3},'bo','markersize',5); % pole with stable frequency and vector
h4=plot(x1{4},y1{4},'gsq','markersize',5);  % pole with stable frequency and damping
h5=plot(x1{5},y1{5},'gx','markersize',5); % pole with stable frequency
if isempty(h1),        h1=0;
elseif isempty(h2),    h2=0;
elseif isempty(h3),    h3=0;
elseif isempty(h4),    h4=0;
elseif isempty(h5),    h5=0;
end



H = [h1(1),h2(1),h3(1),h4(1),h5(1)];
legend(H,...
    'new pole',...
    'stable pole',...
    'stable freq. & MAC',...
    'stable freq. & damp.',...
    'stable freq.',...
    'location','Northoutside','orientation','horizontal');

ylabel('number of poles');
xlabel('f (Hz)')
xlim([0,max([fn{:}])*1.1])
hold off

ax2 = axes('YAxisLocation', 'Right');
linkaxes([ax1,ax2])
plot(ax2,f,Saz./max(Saz).*0.001,'k');
ax2.YLim = [0,Nmax];
ax2.XLim = [0,max([fn{:}])*1.1];
set(ax2,'yscale','log')
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set(gcf,'color','w')



end