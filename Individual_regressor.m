

syms wd1 wd2 wd3 a1 a2 a3 real
syms w1 w2 w3 v1 v2 v3 real

a = [wd1 wd2 wd3 a1 a2 a3]';
v = [w1 w2 w3 v1 v2 v3]';

vv = [v1 v2 v3]';
ww = [w1 w2 w3]';
wd = [wd1 wd2 wd3]';
aa = [a1 a2 a3]';


clear F;
for k = 1:10
 ak = zeros(10,1); ak(k) = 1;
 Ik = inertiaVecToMat(ak);
 Ik(1:3,1:3) = trace(Ik(1:3,1:3))*eye(3) - Ik(1:3,1:3);
 F(:,k) = Ik * a + crf(v)*Ik*v;
end


FF = [F(4:6,:) ; F(1:3,:)]