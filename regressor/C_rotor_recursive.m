function  [H, C, g]  = C_rotor_recursive( model, q , qd, factorFunction)

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v, 1));
end

NJ = length(cell2mat(qd));
a_grav = get_gravity(model);

IC     = model.I;
BC     = repmat({zeros(6,6)},model.NB,1);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  [ XJrot, ~ ] = jcalc( model.jtype{i}, q{i}*model.gr{i} );
  Srot{i} = S{i} * model.gr{i};
  
  vJ    = S{i}*qd{i};
  vJrot = Srot{i}*qd{i};
  
  Xup{i} = XJ * model.Xtree{i};
  Xuprot{i} = XJrot * model.Xrotor{i};
  
  if model.parent(i) == 0
    v{i} = vJ;
    vrot{i} = vJrot;
    ag{i}   = Xup{i} * a_grav;
    agrot{i}= Xuprot{i} * a_grav;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    vrot{i} = Xuprot{i}*v{model.parent(i)} + vJrot;
    ag{i}   = Xup{i} * ag{model.parent(i)};
    agrot{i}= Xuprot{i} * ag{model.parent(i)};
  end
  
  Sd{i} = crm(v{i})*S{i};
  Sdrot{i} = crm(vrot{i})*Srot{i};
   
  BC{i} = factorFunction(model.I{i}, v{i});
  Brot{i} = factorFunction(model.Irot{i}, vrot{i});
end

H     = repmat({},model.NB,model.NB);
C     = repmat({},model.NB,model.NB);
g     = repmat({},model.NB,1);

for i = model.NB:-1:1
  
  
  g{i} = -S{i}'*IC{i}*ag{i} - Srot{i}' * model.Irot{i} * agrot{i};
  
  F1 = IC{i} * S{i};
  H{i,i} = S{i}' * F1 + Srot{i}'*model.Irot{i}*Srot{i};
  F1 = Xup{i}' * F1 + Xuprot{i}'* model.Irot{i}*Srot{i};
  
  F2 = IC{i} * Sd{i} + BC{i} * S{i};
  C{i,i} = S{i}' * F2 + Srot{i}' * (model.Irot{i}*Sdrot{i} + Brot{i} * Srot{i});
  F2 = Xup{i}'*F2 + Xuprot{i}'*(model.Irot{i}*Sdrot{i} + Brot{i} * Srot{i});
  

  
  F3 = ( S{i}'* BC{i} * Xup{i} + Srot{i}'*Brot{i}*Xuprot{i})';
  
  if model.parent(i) ~= 0
    p = model.parent(i);
    IC{p} = IC{p} + Xup{i}' * IC{i} * Xup{i};
    BC{p} = BC{p} + Xup{i}' * BC{i} * Xup{i};
    IC{p} = IC{p} + Xuprot{i}' * model.Irot{i} * Xuprot{i};
    BC{p} = BC{p} + Xuprot{i}' * Brot{i} * Xuprot{i};

    j = i;
    while model.parent(j) ~= 0
        j = model.parent(j);
        if j > 0
            H{j,i} = S{j}' * F1;
            H{i,j} = H{j,i}';
            C{j,i} = S{j}' * F2;
            C{i,j} = (Sd{j}' * F1 + S{j}' * F3)';
        end
        F1 = Xup{j}' * F1;
        F2 = Xup{j}' * F2;
        F3 = Xup{j}' * F3;
    end
  end
end

for i = 1:model.NB
    for j = 1:model.NB
        if size(H{i,j},1) == 0
           H{i,j} = zeros(length(qd{i}), length(qd{j}));
           C{i,j} = zeros(length(qd{i}), length(qd{j}));
        end
    end
end
g = g';