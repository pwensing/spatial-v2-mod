function  [H, Hrot, C, Crot, g, grot ]  = C_rotor( model, q , qd, factorFunction)

if nargin == 3
    factorFunction = @(I,v)(factorFunctions(I,v, 1));
end

NJ = length(cell2mat(qd));


a_grav = get_gravity(model);


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
  
  B{i} = factorFunction(model.I{i}, v{i});
  Brot{i} = factorFunction(model.Irot{i}, vrot{i});
  
end

C = zeros(length(cell2mat(qd)));
H = zeros(length(cell2mat(qd)));
Hrot = zeros(length(cell2mat(qd)));
Crot = zeros(length(cell2mat(qd)));
grot = zeros(length(cell2mat(qd)),1);
g    = zeros(length(cell2mat(qd)),1);


for i = model.NB:-1:1
   J{i}  = cell(1,model.NB);
   Jd{i} = cell(1,model.NB);

   Jrot{i}  = cell(1,model.NB);
   Jrotd{i} = cell(1,model.NB);
   
   X = eye(6);
   Xr= eye(6);
   
   J{i}{i}  = X*S{i};
   Jd{i}{i} = X*Sd{i};
   
   Jrot{i}{i}  = Xr*Srot{i};
   Jrotd{i}{i} = Xr*Sdrot{i};
          
   j = i;       
   X  = X * Xup{j};
   Xr = Xr * Xuprot{j};
   j = model.parent(j);
   
   while j > 0
      J{i}{j} = X*S{j};
      Jd{i}{j} = X*Sd{j};
      
      Jrot{i}{j} = Xr*S{j};
      Jrotd{i}{j} = Xr*Sd{j};
      
      X = X * Xup{j};
      Xr = Xr * Xup{j};
      j = model.parent(j);
   end
   
   for j = 1:model.NB
       if isempty( J{i}{j} )
           J{i}{j} = zeros(6 ,length(qd{j}));
           Jd{i}{j} = zeros(6,length(qd{j}));
           Jrot{i}{j} = zeros(6 ,length(qd{j}));
           Jrotd{i}{j} = zeros(6,length(qd{j}));
       end
   end
   J{i} = cell2mat(J{i});
   Jd{i} = cell2mat(Jd{i});
   
   Jrot{i} = cell2mat(Jrot{i});
   Jrotd{i} = cell2mat(Jrotd{i});
end


for i = 1:model.NB
   H = H + J{i}' * model.I{i} * J{i};
   Hrot = Hrot + Jrot{i}' * model.Irot{i} * Jrot{i};

   Crot = Crot + Jrot{i}' *( Brot{i} * Jrot{i} + model.Irot{i} * Jrotd{i});
   C = C + J{i}' * model.I{i} * Jd{i} + J{i}' * B{i} * J{i}; 
   
   grot = grot - Jrot{i}' * model.Irot{i} * agrot{i} ;
   g    = g    - J{i}'    * model.I{i}    * ag{i}    ;
   
end