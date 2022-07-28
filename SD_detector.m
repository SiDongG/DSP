function [X_hat]=SD_detector(y,H,NT)
% Input parameters
% y : received signal, NRx1
% H : Channel matrix, NRxNT
% NT : number of Tx antennas
% Output parameter
% X_hat : estimated signal, NTx1
global x_list; % candidate symbols in real constellations
global x_now; % temporary x_vector elements
global x_hat; % inv(H)*y
global x_sliced; % sliced x_hat
global x_pre; % x vectors obtained in the previous stage
global real_constellation; % real constellation
global R; % R in the QR decomposition
global radius_squared; % radius^2
global x_metric; % ML metrics of previous stage candidates
global len; % NT*2
QAM_table2 = [-3-3j, -3-j, -3+3j, -3+j, -1-3j, -1-j, -1+3j, -1+j,3-3j,
3-j, 3+3j, 3+j, 1-3j, 1-j, 1+3j, 1+j]/sqrt(10); % 16-QAM
real_constellation = [-3 -1 1 3]/sqrt(10);
y =[real(y); imag(y)]; % y : complex vector -> real vector
H =[real(H) -(imag(H)) ; imag(H) real(H)];
% H : complex vector -> real vector
len = NT*2; % complex -> real
x_list = zeros(len,4); % 4 : real constellation length, 16-QAM
x_now = zeros(len,1);
x_hat = zeros(len,1);
x_pre = zeros(len,1);
x_metric = 0;
[Q,R] = qr(H); % NR x NT QR decomposition
x_hat = inv(H)*y; % zero forcing equalization
x_sliced = QAM16_real_slicer(x_hat,len)’; % slicing
radius_squared = norm(R*(x_sliced-x_hat))^2; % Radious^2
transition = 1; % meaning of transition
% 0 : radius*2, 1len : stage number
% len+1 : compare two vectors in terms of norm values
% len+2 : finish
flag = 1; % transition tracing
% 0 : stage index increases by +1
% 1 : stage index decreases by -1
% 2 : 1->len+2 or len+1->1
while(transition<len+2)
Signal Detection for Spatially Multiplexed MIMO Systems 333
if transition == 0 % radius_squared*2
[flag,transition,radius_squared,x_list]
= radius_control(radius_squared,transition);
elseif transition <= len
[flag,transition] = stage_processing(flag,transition);
elseif transition == len+1 %
[flag,transition] = compare_vector_norm(transition);
end
end
ML = x_pre;
for i=1:len/2
X_hat(i) = ML(i)+j*ML(i+len/2);
end