function v = shift_curve(v_orig,shift)

if 0 == shift 
    v = v_orig';
else
N = length(v_orig);

integer_shift = floor(shift); 
fractional_shift = shift - integer_shift; % fraction of interval

y = zeros(1,3*N);

y(1,N+1:2*N) = v_orig';

y1 = circshift(y,integer_shift*ones(1,numel(y)));
y2 = circshift(y,(integer_shift+1)*ones(1,numel(y)));
y3 = zeros(1,numel(y1));

for k = 1:3*N-1
    
    y3(k) = y1(k) + fractional_shift*( y2(k) - y1(k) );          
    
end;

v = y3(1,N+1:2*N);

v(v<0)=0;
end;

