% distribution histogram
%{
load pos_times.txt;
load neg_times.txt;

hist(pos_times(3:130), 20);
%% hold on
hist(neg_times(3:131), 20);
%}
%%

a = load(['pos_times_' num2str(1) '.txt']);
n = load(['neg_times_' num2str(1) '.txt']);
a = a(2:end-1);
n = n(2:end-1);

for i = 1:13
    b = load(['pos_times_' num2str(i) '.txt']);
    m = load(['neg_times_' num2str(i) '.txt']);

    a = cat(1, a, b(2:end-1));
    n = cat(1, n, m(2:end-1));
end

%hist(a, 50);
dur1 = a;
bins = 50;
a1 = 180;
x1 = 0:1:a1;
phat = gamfit(dur1);
y1 = gampdf(x1,phat(1),phat(2));
histogram(dur1,bins);
plot(x1, y1,'red')
%%
% TRY USING GAMFIT INSTEAD
% pdp = fitdist(a, 'Gamma');
% pdn = fitdist(n, 'Gamma');
