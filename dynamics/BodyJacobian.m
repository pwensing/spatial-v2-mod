function  J = BodyJacobian( model, q, qd, body_num, Xend)

J = cell(1,model.NB);

X = Xend;
j = body_num;
while j > 0
    [ XJ, S ] = jcalc( model.jtype{j}, q{j} );
    J{j} = X * S;
    if model.parent(j) > 0
        X = X *  XJ * model.Xtree{ j };
    end
    j = model.parent(j);
end

for i = 1:model.NB
    if isempty(J{i})
        J{i} = zeros(size(Xend,1), length(qd{i}));
    end
end
J = cell2mat(J);



