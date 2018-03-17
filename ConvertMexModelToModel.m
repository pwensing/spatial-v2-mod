function [model] = ConvertMexModelToModel(model, jtype, Xtree, I, Irot, Xrotor)
    model.parent =  max( model.parent-5,0);
    model.NB = model.NB-5;
    model.parent = model.parent(6:end);
    model.gc.body = model.gc.body-5;
    jtype{6} = 'fb';
    a=cell(model.NB,1);
    gr = model.gr;
    model.gr = {};
   
    for i = 1:model.NB
        model.I{i}     = I{i+5};
        model.jtype{i} = jtype{i+5};
        model.Xtree{i} = Xtree{i+5};
        model.Irot{i} = Irot{i+5};
        model.Xrotor{i} = Xrotor{i+5};
        model.gr{i}    = gr(i+5);
    end 
    model.gr{1}=1;
end
