function [ a ] = pinertiaToVec( J )
    m = J(4,4);
    h = J(1:3,4);
    c = h/m;
    E = J(1:3,1:3);
    Ibar = trace(E)*eye(3) - E;
    
    a = inertiaMatToVec([Ibar skew(h) ; skew(h)' m*eye(3)]);
end

