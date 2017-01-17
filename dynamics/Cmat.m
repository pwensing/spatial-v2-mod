function  [C, H, Hd] = Cmat( model, q, qd, factorFunction)

IC = model.I;

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v, 3));
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
  end
  Sd{i} = crm(v{i}) * S{i};
  IC{i} = model.I{i};
  BC{i} = factorFunction(model.I{i}, v{i});
  ICd{i} = crf(v{i}) * model.I{i} - model.I{i} *crm(v{i});
end

C = cell(model.NB, model.NB);
H = C;
Hd = C;


for i = model.NB:-1:1
  if model.parent(i) ~= 0
     p = model.parent(i);
     IC{p} = IC{p} + Xup{i}' * IC{i} * Xup{i};
     ICd{p} = ICd{p} + Xup{i}' * ICd{i} * Xup{i};
     BC{p} = BC{p} + Xup{i}' * BC{i} * Xup{i};
  end
  F1 = IC{i} * Sd{i} + BC{i} * S{i};
  F2 = IC{i} * S{i};
  F3 = (S{i}'*BC{i})';
  F2d = ICd{i}*S{i} + IC{i} *Sd{i};
 
  C{i,i} = S{i}'*F1;
  H{i,i} = S{i}'*F2;
  Hd{i,i} = Sd{i}' * F2 + S{i}'*F2d;
  
  j = i;
  while model.parent(j) ~= 0
      F1 = Xup{j}'*F1;
      F2 = Xup{j}'*F2;
      F3 = Xup{j}'*F3;
      F2d = Xup{j}'*F2d;
      
      
      j = model.parent(j);
      if j > 0
          C{j,i} = S{j}' * F1;
          C{i,j} = (Sd{j}'*F2 + S{j}'*F3)';
          H{j,i} = S{j}' * F2;
          H{i,j} = H{j,i}';
          Hd{j,i} =Sd{j}' * F2 + S{j}'*F2d;
          Hd{i,j} = Hd{j,i}';
      end
  end
end

for i = 1:model.NB
    for j = 1:model.NB
        if size(H{i,j},1) == 0
           H{i,j} = zeros(length(qd{i}), length(qd{j}));
           Hd{i,j} = zeros(length(qd{i}), length(qd{j}));
           C{i,j} = zeros(length(qd{i}), length(qd{j}));
        end
    end
end

