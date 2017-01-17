function  Y = SlotineLiY( model, qddr, qdr,qd,q, factorFunction)

if nargin == 5
    factorFunction = @(I,v)(factorFunctions(I,v, 3));
end


for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  vJr= S{i}*qdr{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    w{i} = vJr;
    wd{i}= -Xup{i}*model.a_grav + S{i} * qddr{i} + crm(v{i})*vJr;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    w{i} = Xup{i}*w{model.parent(i)} + vJr;
    wd{i}= Xup{i}*wd{model.parent(i)} + S{i} * qddr{i} + crm(v{i})*vJr;
  end
  F{i} = zeros(6,10);
  for k = 1:10
     ak = zeros(10,1); ak(k) = 1;
     Ik = inertiaVecToMat(ak);
     F{i}(:,k) = Ik * wd{i} + factorFunction(Ik, v{i}) * w{i};
  end
end

Y = cell(model.NB, model.NB);

for i = model.NB:-1:1
  Y{i,i} = S{i}' * F{i};
  Fji = F{i};
  j = i;
  while model.parent(j) ~= 0
      Fji = Xup{j}' * Fji;
      j = model.parent(j);
      Y{j,i} = S{j}' * Fji;
  end
end

for i = 1:model.NB
    for j = 1:model.NB
        if size(Y{i,j},1) == 0
           Y{i,j} = zeros(length(qd{i}), 10);
        end
    end
end

Y = cell2mat(Y);


